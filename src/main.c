#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#define GL_SILENCE_DEPRECATION
#include "../dependencies/glew/2.2.01/include/GL/glew.h"
#include "../dependencies/glfw/3.4/include/GLFW/glfw3.h"
#include "../dependencies/glut/include/glut.h"
#include "fem.h"
#include "formation_mtrx.h"

#define LOAD 25000.  // Приложенная нагрузка
#define FLT_MAX 100.
#define KOEF_X 25.
#define KOEF_Y 14.
#define WINDTH_WINDOW 1200          // Высота окна
#define HEIGHT_WINDOW 800           // Ширина окна
#define STRESS_SCALE 1000.0         // Масштаб для визуализации напряжений
#define DEFORMATION_SCALE 10.0      // Масштаб для визуализации деформаций
#define GLUT_WINDOW_POSITION_X 100  // Позиция окна (по оси x)
#define GLUT_WINDOW_POSITION_Y 100  // Позиция окна (по оси y)

// Глобальные переменные для управления визуализацией
static int showDeformed = 1;           // Показывать деформированную модель
static int showStress = 1;             // Показывать напряжения
static int showValues = 0;             // Показывать числовые значения
static double zoom = 1.0;              // Масштаб отображения
static float animationProgress = 0.f;  // Шаг анимации
static int isAnimating = 0;            // Флаг анимации
static int animationDirection = 1;     // Направление анимации

void drawMashForSolve(int argc, char **argv);
void drawModel(void);
void display(void);
void init(void);
void updateAnimation(int __attribute__((unused)) value);
void keyboard(unsigned char key, int __attribute__((unused)) x,
              int __attribute__((unused)) y);
void renderText(float x, float y, const char *text);

static int nelem;               // кол-во треугольных элементов
static int nys;                 // число узлов К.Э модели
static double **car = NULL;     // массив координат узлов элемента
static int **jt03 = NULL;       // массив номеров узлов элемента
static double *u = NULL;        // массив перемещений узлов
static double **stress = NULL;  // массив напряжений

int main(int argc, char **argv) {
  int ndofysla = 2;  // кол-во степеней свободы одного узла
  double *dataCar;
  int *data_jt03;
  short fileErr = readFromFile("../nodes/node1.txt", &nys, &dataCar, &car,
                               &nelem, &data_jt03, &jt03);
  if (fileErr == 1) {
    free_memory(4, dataCar, car, data_jt03, jt03);
    exit(1);
  } else if (fileErr == 2) {
    free_memory(3, car, data_jt03, jt03);
    exit(1);
  } else if (fileErr == 3) {
    free_memory(3, dataCar, car, jt03);
    exit(1);
  }
  int ndof = nys * ndofysla;  // общее число степеней свободы
  // задаем материал
  double h = 1.0;
  double e = 2.1e5;
  double puas = 0.3;
  // глобальная матрица жесткости kglb[ndof][ndof]
  double *dataKGLB = (double *)calloc(ndof * ndof, sizeof(double));
  double **kglb = (double **)calloc(ndof, sizeof(double *));
  for (int i = 0; i < ndof; i++) {
    kglb[i] = dataKGLB + i * ndof;
  }
  if (kglb == NULL) {
    free_memory(5, kglb, dataCar, car, data_jt03, jt03);
    exit(1);
  }
  u = (double *)malloc(ndof * sizeof(double));  // массив перемещений узлов
  if (u == NULL) {
    free_memory(6, kglb, dataCar, car, data_jt03, jt03, u);
    exit(1);
  }
  double *r = (double *)malloc(ndof * sizeof(double));  // массив нагрузок
  // массив x (рабочий LDLT)
  double *x = (double *)malloc(ndof * sizeof(double));
  // расчет матрицы лок. жесткости и добавление ее в глоб. матрицу
  AssembleLocalStiffnessToGlobal(kglb, jt03, car, nelem, e, h, puas, ndofysla);
  int lenNodePres = 0, lenNodeZakrU = 0, lenNodeZakrV = 0;
  int *nodePres = NULL;   // массив нагруженных узлов
  int *nodeZakrU = NULL;  // массив закрепленных узлов по X
  int *nodeZakrV = NULL;  // массив закрепленных узлов по Y
  // нахождение закрепленных и нагруженных узлов
  FillConstrainedLoadedNodes(&nodePres, &lenNodePres, &nodeZakrU, &lenNodeZakrU,
                             &nodeZakrV, &lenNodeZakrV, car, nys);
  SetLoadVector(r, lenNodePres, nodePres, ndofysla, ndof,
                LOAD);  // задаем вектор нагрузок
  MakeConstrained(nodeZakrV, lenNodeZakrV, kglb,
                  ndofysla);  // задаем закрепления по V
  MakeConstrained(nodeZakrU, lenNodeZakrU, kglb,
                  ndofysla);  // задаем закрепления по U
  // решение СЛАУ методом разложения в LDLT
  bool ierr = solveLinearSystemLDLT(kglb, u, r, x, ndof);
  if (ierr) {  // ошибка разложения в LDLT или диаагонального решения
    free_memory(12, nodePres, nodeZakrU, nodeZakrV, u, r, x, dataKGLB, kglb,
                dataCar, car, data_jt03, jt03);
    exit(1);
  }
  // расчет деформаций, напряжений
  double *dataStrain = NULL;  // массив деформаций
  double **strain = NULL;
  makeDoubleMtrx(&dataStrain, &strain, 4, nelem);
  if (strain == NULL) {
    free_memory(13, strain, nodePres, nodeZakrU, nodeZakrV, u, r, x, dataKGLB,
                kglb, dataCar, car, data_jt03, jt03);
    exit(1);
  }
  double *dataStress = NULL;  // массив напряжений
  makeDoubleMtrx(&dataStress, &stress, 4, nelem);
  if (stress == NULL) {
    free_memory(15, stress, dataStrain, strain, nodePres, nodeZakrU, nodeZakrV,
                u, r, x, dataKGLB, kglb, dataCar, car, data_jt03, jt03);
    exit(1);
  }
  stressModel(ndofysla, nelem, jt03, car, e, puas, u, strain, stress);
  writeResult("result.txt", jt03, strain, stress, r, u, nelem, nys, ndof);
  drawMashForSolve(argc, argv);  // отрисовка модели, разбитой на КЭ
  // освобождение памяти из под матрицы
  free_memory(16, dataStress, stress, dataStrain, strain, nodePres, nodeZakrU,
              nodeZakrV, u, r, x, dataKGLB, kglb, dataCar, car, data_jt03,
              jt03);
}

void drawMashForSolve(int argc, char **argv) {
  // Инициализация GLUT
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(WINDTH_WINDOW, HEIGHT_WINDOW);
  glutInitWindowPosition(GLUT_WINDOW_POSITION_X, GLUT_WINDOW_POSITION_Y);
  glutCreateWindow("FEM Visualization");
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  init();
  glutMainLoop();
}

// функция отрисовки конечных элементов модели
void drawModel(void) {
  for (int i = 0; i < nelem; i++) {
    // Вычисление среднего напряжения для элемента
    int ielem = jt03[0][i];
    double avgStress = 0.0;
    for (int j = 0; j < 3; j++) {
      avgStress += stress[j][ielem - 1];
    }
    avgStress /= 3.0;
    // Нормализация напряжения для цветовой схемы
    double normalizedStress = fabs(avgStress) / STRESS_SCALE;
    normalizedStress = (normalizedStress > 1.0) ? 1.0 : normalizedStress;
    // Установка цвета в зависимости от режима отображения
    if (showStress) {
      glColor3f(normalizedStress, 0.0, 1.0 - normalizedStress);
    } else {
      glColor3f(0.0, 0.5, 0.0);
    }
    // Отрисовка треугольника
    glBegin(GL_TRIANGLES);
    for (int j = 0; j < 3; j++) {
      double x = car[0][jt03[j][i] - 1];
      double y = car[1][jt03[j][i] - 1];
      if (showDeformed) {
        if (x != 0. && y != 0) {  // опираясь на схему нагружения образца
          x += u[jt03[j][i] * 2 - 2] * DEFORMATION_SCALE * animationProgress;
          y += u[jt03[j][i] * 2 - 1] * DEFORMATION_SCALE * animationProgress;
        }
      }
      glVertex2f(KOEF_X + x * zoom / 2., KOEF_Y + y * zoom / 2.);
    }
    glEnd();
    // Отрисовка границ элементов
    glColor3f(0.7, 0.7, 0.7);
    int index = 0;
    for (int j = 0; j < 3; j++) {
      glBegin(GL_LINES);
      double x1 = car[0][jt03[index][i] - 1];
      double y1 = car[1][jt03[index][i] - 1];
      double x2 = car[0][jt03[(index + 1) % 3][i] - 1];
      double y2 = car[1][jt03[(index + 1) % 3][i] - 1];
      if (showDeformed) {
        if (x1 != 0. && y1 != 0) {
          x1 +=
              u[jt03[index][i] * 2 - 2] * DEFORMATION_SCALE * animationProgress;
          y1 +=
              u[jt03[index][i] * 2 - 1] * DEFORMATION_SCALE * animationProgress;
        }
        if (x2 != 0. && y2 != 0) {
          x2 += u[jt03[(index + 1) % 3][i] * 2 - 2] * DEFORMATION_SCALE *
                animationProgress;
          y2 += u[jt03[(index + 1) % 3][i] * 2 - 1] * DEFORMATION_SCALE *
                animationProgress;
        }
      }
      glVertex2f(KOEF_X + x1 * zoom / 2., KOEF_Y + y1 * zoom / 2.);
      glVertex2f(KOEF_X + x2 * zoom / 2., KOEF_Y + y2 * zoom / 2.);
      glEnd();
      index = (index + 1) % 3;
    }
    // Отображение числовых значений напряжений
    if (showValues) {
      char text[32];
      sprintf(text, "%.1f", avgStress);
      double centerX = 0.0, centerY = 0.0;
      for (int j = 0; j < 3; j++) {
        centerX += car[0][jt03[j][i] - 1];
        centerY += car[1][jt03[j][i] - 1];
      }
      centerX = KOEF_X + (centerX / 3.0) * zoom / 2.;
      centerY = KOEF_Y + (centerY / 3.0) * zoom / 2.;
      glColor3f(0.0, 0.0, 0.0);
      renderText(centerX, centerY, text);
    }
  }
}

void display(void) {
  glClear(GL_COLOR_BUFFER_BIT);
  drawModel();
  glFlush();
}

void init(void) {
  double minX = FLT_MAX, maxX = -FLT_MAX;
  double minY = FLT_MAX, maxY = -FLT_MAX;
  for (int i = 0; i < nys; i++) {
    if (car[0][i] < minX) minX = car[0][i];
    if (car[0][i] > maxX) maxX = car[0][i];
    if (car[1][i] < minY) minY = car[1][i];
    if (car[1][i] > maxY) maxY = car[1][i];
  }
  gluOrtho2D(minX, maxX, minY, maxY);
  glClearColor(1., 1., 1., 1.);  // Белый фон
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
}
void updateAnimation(int __attribute__((unused)) value) {
  if (isAnimating) {
    animationProgress += 0.01f * animationDirection;
    if (animationProgress >= 1.0f) {
      animationProgress = 1.0f;
      animationDirection = -1;
    } else if (animationProgress <= 0.0f) {
      animationProgress = 0.0f;
      animationDirection = 1;
    }
    glutPostRedisplay();
    glutTimerFunc(16, updateAnimation, 0);
  }
}

void keyboard(unsigned char key, int __attribute__((unused)) x,
              int __attribute__((unused)) y) {
  switch (key) {
    case 'd':  // Переключение деформированной/недеформированной модели
      showDeformed = !showDeformed;
      break;
    case 's':  // Переключение отображения напряжений
      showStress = !showStress;
      break;
    case 'v':  // Переключение отображения числовых значений
      showValues = !showValues;
      break;
    case '=':  // Увеличение масштаба (чтобы не зажимать shift)
      zoom *= 1.1;
      break;
    case '-':  // Уменьшение масштаба
      zoom /= 1.1;
      break;
    case 'a':  // Включение/выключение анимации
      isAnimating = !isAnimating;
      if (isAnimating) {
        animationProgress = 0.0f;
        animationDirection = 1;
        glutTimerFunc(16, updateAnimation, 0);
      }
      break;
  }
  glutPostRedisplay();
}

// Функция для отображения текста
void renderText(float x, float y, const char *text) {
  glRasterPos2f(x, y);
  for (const char *c = text; *c != '\0'; c++) {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, *c);
  }
}
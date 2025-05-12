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

#define LOAD 25000.
#define FLT_MAX 100.
#define KOEF_X 25.
#define KOEF_Y 14.
#define WINDTH_WINDOW 1200
#define HEIGHT_WINDOW 800
#define STRESS_SCALE 1000.0     // Масштаб для визуализации напряжений
#define DEFORMATION_SCALE 10.0  // Масштаб для визуализации деформаций

// Глобальные переменные для управления визуализацией
static int showDeformed = 1;  // Показывать деформированную модель
static int showStress = 1;    // Показывать напряжения
static int showValues = 0;    // Показывать числовые значения
static double zoom = 1.0;     // Масштаб отображения

void drawMashForSolve(int argc, char **argv);
void drawModel(void);
void display(void);
void init(void);
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
  double *dataGEST;
  // локальная матрица жесткости gest[6][6]
  double **gest = NULL;
  makeDoubleMtrx(&dataGEST, &gest, 6, 6);
  if (gest == NULL) {
    free_memory(5, gest, dataCar, car, data_jt03, jt03);
    exit(1);
  }
  // глобальная матрица жесткости kglb[ndof][ndof]
  double *dataKGLB = (double *)calloc(ndof * ndof, sizeof(double));
  double **kglb = (double **)calloc(ndof, sizeof(double *));
  for (int i = 0; i < ndof; i++) {
    kglb[i] = dataKGLB + i * ndof;
  }
  if (kglb == NULL) {
    free_memory(7, kglb, dataGEST, gest, dataCar, car, data_jt03, jt03);
    exit(1);
  }
  u = (double *)malloc(ndof * sizeof(double));  // массив перемещений узлов
  if (u == NULL) {
    free_memory(8, kglb, dataGEST, gest, dataCar, car, data_jt03, jt03, u);
    exit(1);
  }
  double *r = (double *)malloc(ndof * sizeof(double));  // массив нагрузок
  // массив x (рабочий LDLT)
  double *x = (double *)malloc(ndof * sizeof(double));
  // расчет матрицы лок. жесткости и добавление ее в глоб. матрицу
  AssembleLocalStiffnessToGlobal(gest, kglb, jt03, car, nelem, e, h, puas,
                                 ndofysla);
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
    free_memory(14, nodePres, nodeZakrU, nodeZakrV, u, r, x, dataKGLB, kglb,
                dataGEST, gest, dataCar, car, data_jt03, jt03);
    exit(1);
  }
  // расчет деформаций, напряжений
  double *dataStrain = NULL;  // массив деформаций
  double **strain = NULL;
  makeDoubleMtrx(&dataStrain, &strain, 4, nelem);
  if (strain == NULL) {
    free_memory(15, strain, nodePres, nodeZakrU, nodeZakrV, u, r, x, dataKGLB,
                kglb, dataGEST, gest, dataCar, car, data_jt03, jt03);
    exit(1);
  }
  double *dataStress = NULL;  // массив напряжений
  makeDoubleMtrx(&dataStress, &stress, 4, nelem);
  if (stress == NULL) {
    free_memory(17, stress, dataStrain, strain, nodePres, nodeZakrU, nodeZakrV,
                u, r, x, dataKGLB, kglb, dataGEST, gest, dataCar, car,
                data_jt03, jt03);
    exit(1);
  }
  stressModel(ndofysla, nelem, jt03, car, e, puas, u, strain, stress);
  writeResult("result.txt", jt03, strain, stress, r, u, nelem, nys, ndof);
  drawMashForSolve(argc, argv);  // отрисовка модели, разбитой на КЭ
  // освобождение памяти из под матрицы
  free_memory(18, dataStress, stress, dataStrain, strain, nodePres, nodeZakrU,
              nodeZakrV, u, r, x, dataKGLB, kglb, dataGEST, gest, dataCar, car,
              data_jt03, jt03);
}

void drawMashForSolve(int argc, char **argv) {
  // Инициализация GLUT
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(WINDTH_WINDOW, HEIGHT_WINDOW);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("FEM Visualization");
  // Регистрация callback-функций
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  // Отрисовка модели
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
          x += u[jt03[j][i] * 2 - 2] * DEFORMATION_SCALE;
          y += u[jt03[j][i] * 2 - 1] * DEFORMATION_SCALE;
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
          x1 += u[jt03[index][i] * 2 - 2] * DEFORMATION_SCALE;
          y1 += u[jt03[index][i] * 2 - 1] * DEFORMATION_SCALE;
        }
        if (x2 != 0. && y2 != 0) {
          x2 += u[jt03[(index + 1) % 3][i] * 2 - 2] * DEFORMATION_SCALE;
          y2 += u[jt03[(index + 1) % 3][i] * 2 - 1] * DEFORMATION_SCALE;
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
#include "formation_mtrx.h"

#include <stdio.h>
#include <stdlib.h>

#include "fem.h"

// расчет матрицы деформаций
double **formationDeformMtrx(double **deformMtrx, coord coord1, coord coord2,
                             coord coord3, double a2) {
  double a;
  a = (a2 > 0.) ? 1 / (a2) : 1.;
  deformMtrx[0][0] = a * (coord2.y - coord3.y);
  deformMtrx[0][1] = 0;
  deformMtrx[0][2] = a * (coord3.y - coord1.y);
  deformMtrx[0][3] = 0;
  deformMtrx[0][4] = a * (coord1.y - coord2.y);
  deformMtrx[0][5] = 0;
  deformMtrx[1][0] = 0;
  deformMtrx[1][1] = a * (coord3.x - coord2.x);
  deformMtrx[1][2] = 0;
  deformMtrx[1][3] = a * (coord1.x - coord3.x);
  deformMtrx[1][4] = 0;
  deformMtrx[1][5] = a * (coord2.x - coord1.x);
  deformMtrx[2][0] = a * (coord3.x - coord2.x);
  deformMtrx[2][1] = a * (coord2.y - coord3.y);
  deformMtrx[2][2] = a * (coord1.x - coord3.x);
  deformMtrx[2][3] = a * (coord3.y - coord1.y);
  deformMtrx[2][4] = a * (coord2.x - coord1.x);
  deformMtrx[2][5] = a * (coord1.y - coord2.y);
  return deformMtrx;
}

// расчет матрицы упругости
double **formationElastMtrx(double **elastMtrx, double e, double puas) {
  double ekoef = e / (1 - puas * puas);
  elastMtrx[0][0] = ekoef;
  elastMtrx[0][1] = ekoef * puas;
  elastMtrx[0][2] = 0;
  elastMtrx[1][0] = ekoef * puas;
  elastMtrx[1][1] = ekoef;
  elastMtrx[1][2] = 0;
  elastMtrx[2][0] = 0;
  elastMtrx[2][1] = 0;
  elastMtrx[2][2] = ekoef * (1 - puas) / 2;
  return elastMtrx;
}

// расчет матрицы напряжений
void stressPlanElem(coord coord1, coord coord2, coord coord3, double e,
                    double puas, double **deformMtrx, double **strsMatr) {
  // формирование матрицы упругости elastMtrx[3][3]
  double *dataElastMtrx = (double *)malloc(3 * 3 * sizeof(double));
  double **elastMtrx = NULL;
  makeDoubleMtrx(&dataElastMtrx, &elastMtrx, 3, 3);
  // заполнение матрицы упругости
  elastMtrx = formationElastMtrx(elastMtrx, e, puas);
  // заполнение матрицы деформаций deformMtrx[3][6]
  double a2 = coord2.x * coord3.y - coord3.x * coord2.y - coord1.x * coord3.y +
              coord1.y * coord3.x + coord1.x * coord2.y - coord1.y * coord2.x;
  deformMtrx = formationDeformMtrx(deformMtrx, coord1, coord2, coord3, a2);
  // заполнение матрицы напряжений
  // strsMatr[3][6]=elastMtrx[3][3]*deformMtrx[3][6]
  double sum;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 6; j++) {
      sum = 0.;
      for (int k = 0; k < 3; k++) {
        sum += elastMtrx[i][k] * deformMtrx[k][j];
      }
      strsMatr[i][j] = sum;
    }
  }
  free_memory(2, dataElastMtrx, elastMtrx);
}

// заполнение матрицы жесткости gest[6][6] трегольного конечного элемента
void planeElement(coord coord1, coord coord2, coord coord3, double e, double h,
                  double puas, double **gest) {
  double sum;
  // формирование матрица деформаций deformMtrx[3][6]
  double *dataDefMtrx = (double *)malloc(3 * 6 * sizeof(double));
  double **deformMtrx = NULL;
  makeDoubleMtrx(&dataDefMtrx, &deformMtrx, 3, 6);
  // заполнение матрицы деформаций
  double a2 = coord2.x * coord3.y - coord3.x * coord2.y - coord1.x * coord3.y +
              coord1.y * coord3.x + coord1.x * coord2.y - coord1.y * coord2.x;
  // формирование матрицы напряжений strsMatr[3][6]
  double *dataStrsMatr = (double *)malloc(3 * 6 * sizeof(double));
  double **strsMatr = NULL;
  makeDoubleMtrx(&dataStrsMatr, &strsMatr, 3, 6);
  // заполнение матрицы напряжений
  stressPlanElem(coord1, coord2, coord3, e, puas, deformMtrx, strsMatr);
  double vol = h * a2 * 0.5;
  // вычисление локальной матрицы жесткости
  // gest[6][6]=deformMtrx(trans)[6][3]*strsMatr[3][6];
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      sum = 0.;
      for (int k = 0; k < 3; k++) {
        sum += deformMtrx[k][i] * strsMatr[k][j];
      }
      gest[i][j] = sum * vol;
    }
  }
  // освобождаем память
  free_memory(4, dataDefMtrx, dataStrsMatr, deformMtrx, strsMatr);
}

// функция сборки глобальной матрицы жесткости kglb[ndof][ndof]
void assemblyGlobMatr(int ndofysla, nodeNumber node, double **gest,
                      double **kglb) {
  int il, jl, ig, jg;  // начальные позиции в лок. и глоб. матрицах
  int iblok, jblok;    // добавочные коэф. к позициям матриц
  int nys[3];
  nys[0] = node.iys1 - 1;
  nys[1] = node.iys2 - 1;
  nys[2] = node.iys3 - 1;
  for (int iy = 0; iy < 3; iy++) {
    for (int jy = 0; jy < 3; jy++) {
      il = iy * ndofysla;
      jl = jy * ndofysla;
      ig = nys[iy] * ndofysla;
      jg = nys[jy] * ndofysla;
      for (iblok = 0; iblok < ndofysla; iblok++) {
        for (jblok = 0; jblok < ndofysla; jblok++) {
          kglb[ig + iblok][jg + jblok] += gest[il + iblok][jl + jblok];
        }
      }
    }
  }
}

// функция заполнения массивов закрепленных узлов и массив нагруженных узлов
void FillConstrainedLoadedNodes(int **nodePres, int *lenNodePres,
                                int **nodeZakrU, int *lenNodeZakrU,
                                int **nodeZakrV, int *lenNodeZakrV,
                                double **car, int nys) {
  for (int i = 0; i < nys; i++) {
    if ((int)car[0][i] == 100) {
      (*lenNodePres)++;
      if (!(*lenNodePres)) {
        *nodePres = (int *)malloc(sizeof(int));
        (*nodePres)[*lenNodePres] = i;
      } else {
        *nodePres = (int *)realloc(*nodePres, *lenNodePres * sizeof(int));
        (*nodePres)[*lenNodePres - 1] = i;
      }
    }
    if (!(int)car[0][i]) {
      (*lenNodeZakrU)++;
      if (!(*lenNodeZakrU)) {
        *nodeZakrU = (int *)malloc(sizeof(int));
        (*nodeZakrU)[*lenNodeZakrU] = i;
      } else {
        *nodeZakrU = (int *)realloc(*nodeZakrU, *lenNodeZakrU * sizeof(int));
        (*nodeZakrU)[*lenNodeZakrU - 1] = i;
      }
    }
    if (!(int)car[1][i]) {
      (*lenNodeZakrV)++;
      if (!(*lenNodeZakrV)) {
        *nodeZakrV = (int *)malloc(sizeof(int));
        (*nodeZakrV)[*lenNodeZakrV] = i;
      } else {
        *nodeZakrV = (int *)realloc(*nodeZakrV, *lenNodeZakrV * sizeof(int));
        (*nodeZakrV)[*lenNodeZakrV - 1] = i;
      }
    }
  }
}

void MakeConstrained(int *nodeZakr, int lenNodeZakr, double **kglb,
                     int ndofysla) {
  for (int i = 0; i < lenNodeZakr; i++) {
    int kdof = nodeZakr[i] * ndofysla;
    kglb[kdof][kdof] += 1.e38;
  }
}

void SetLoadVector(double *r, int lenNodePres, int *nodePres, int ndofysla,
                   int ndof, float load) {
  for (int i = 0; i < ndof; i++) {
    r[i] = 0.;
  }
  for (int i = 0; i < lenNodePres; i++) {
    r[nodePres[i] * ndofysla] = load / lenNodePres;
  }
}

void makeDoubleMtrx(double **dataMtrx, double ***mtrx, int row, int col) {
  *dataMtrx = (double *)malloc(row * col * sizeof(double));
  if (*dataMtrx == NULL) {
    printf("Can't allocate memory for dataMtrx\n");
    exit(1);
  }
  *mtrx = (double **)malloc(row * sizeof(double *));
  if (*mtrx == NULL) {
    printf("Can't allocate memory for row pointers\n");
    free(*dataMtrx);
    exit(1);
  }
  for (int i = 0; i < row; i++) {
    (*mtrx)[i] = *dataMtrx + i * col;
  }
}

void makeIntegerMtrx(int **dataMtrx, int ***mtrx, int row, int col) {
  *dataMtrx = (int *)malloc(row * col * sizeof(int));
  if (*dataMtrx == NULL) {
    printf("Can't allocate memory for dataMtrx\n");
    exit(1);
  }
  *mtrx = (int **)malloc(row * sizeof(int *));
  if (*mtrx == NULL) {
    printf("Can't allocate memory for row pointers\n");
    free(*dataMtrx);
    exit(1);
  }
  for (int i = 0; i < row; i++) {
    (*mtrx)[i] = *dataMtrx + i * col;
  }
}

short readFromFile(char *filename, int *nys, double **dataCar, double ***car,
                   int *nelem, int **data_jt03, int ***jt03) {
  bool err = 0;
  FILE *file = fopen(filename, "r");
  if (file == NULL) {
    err = 1;
  } else {
    fscanf(file, "%d", nys);
    makeDoubleMtrx(dataCar, car, 3, *nys);  // массив координат узлов элемента
  }
  if (*car == NULL) {
    err = 2;
  }
  for (int i = 0; i < *nys; i++) {
    fscanf(file, "%lf%lf%lf", &(*car)[0][i], &(*car)[1][i], &(*car)[2][i]);
  }
  fscanf(file, "%d", nelem);
  makeIntegerMtrx(data_jt03, jt03, 3, *nelem);  // массив номеров узлов элемента
  if (*jt03 == NULL) {
    err = 3;
  }
  for (int i = 0; i < *nelem; i++) {
    fscanf(file, "%d%d%d", &(*jt03)[0][i], &(*jt03)[1][i], &(*jt03)[2][i]);
  }

  fclose(file);
  return err;
}

void writeResult(char *filename, int **jt03, double **strain, double **stress,
                 double *r, double *u, int nelem, int nys, int ndof) {
  FILE *file = fopen(filename, "w");
  fprintf(file, "Число элементов - %d\n", nelem);
  fprintf(file, "Число узлов - %d\n", nys);
  fprintf(file, "Число степеней свободы - %d\n", ndof);
  fprintf(file, "\n");
  fprintf(file, "Вектор нагрузок\n");
  int index = 1;
  for (int i = 0; i <= ndof; i += 2) {
    fprintf(file, "       ru%d       %12.4e       rv%d       %12.4e\n", index,
            r[i], index, r[i + 1]);
    index++;
  }
  fprintf(file, "\n");
  fprintf(file, "Результат расчета перемещений\n");
  index = 1;
  for (int i = 0; i <= ndof; i += 2) {
    fprintf(file, "       u%d       %12.4e       v%d       %12.4e\n", index,
            u[i], index, u[i + 1]);
    index++;
  }
  fprintf(file, "\n");
  fprintf(file, "Результат расчета деформаций, напряжений\n");
  for (int j = 0; j < nelem; j++) {
    int ielem = jt03[0][j];
    fprintf(file, "       %d       ", ielem);
    for (int i = 0; i < 4; i++) {
      fprintf(file, "       %12.4f", strain[i][ielem - 1]);
    }
    fprintf(file, "       |");
    for (int i = 0; i < 4; i++) {
      fprintf(file, "       %12.4f", stress[i][ielem - 1]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
}
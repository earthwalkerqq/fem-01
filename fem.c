#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "fem.h"
#include "formation_mtrx.h"

// расчет матрицы лок. жесткости и добавление ее в глоб. матрицу
void AssembleLocalStiffnessToGlobal(double **gest, double **kglb, int **jt03,
                                    double **car, int nelem, double e, double h,
                                    double puas, int ndofysla) {
#pragma omp parallel shared(nelem, jt03, car)
  {
#pragma omp for
    for (int ielem = 0; ielem < nelem; ielem++) {
      nodeNumber node = {
          jt03[0][ielem],
          jt03[1][ielem],
          jt03[2][ielem],
      };

      coord coord1 = {car[0][node.iys1 - 1], car[1][node.iys1 - 1]};
      coord coord2 = {car[0][node.iys2 - 1], car[1][node.iys2 - 1]};
      coord coord3 = {car[0][node.iys3 - 1], car[1][node.iys3 - 1]};
      planeElement(coord1, coord2, coord3, e, h, puas, gest);
      assemblyGlobMatr(ndofysla, node, gest, kglb);
    }
  }
}

// разложение матрицы глобальной матрицы жесткости kglb[ndof][ndof] в LDLT
bool matrLDLT(int ndof, double **kglb) {
  double diag;
  int i = 0;
  bool ierr = false;
  double sum = kglb[i][i];
  if (fabs(sum) < 1.0e-20) {
    printf("Разложение невозможно. Нулевой диигольный элемент");
    return true;
  }
#pragma omp parallel shared(kglb, sum, ndof)
  {
#pragma omp for
    for (int j = 1; j < ndof; j++) {
      kglb[j][0] /= sum;
      kglb[0][j] = 0.0;
    }
  }
  for (i = 1; i < ndof - 1; i++) {
    diag = kglb[i][i];
    for (int j = 0; j < i; j++) {
      diag -= kglb[j][j] * (kglb[i][j] * kglb[i][j]);
    }
    kglb[i][i] = diag;
    if (fabs(diag) < 1.e-20) {
      printf("Разложение невозможно. Нулевой диагональный элемент - %d ",
             i + 1);
      ierr = true;
      return ierr;
    }
#pragma omp parallel shared(kglb, i, ndof, diag) private(sum)
    {
#pragma omp for
      for (int k = i + 1; k < ndof; k++) {
        sum = kglb[k][i];
        for (int j = 0; j < i; j++) {
          sum -= kglb[j][j] * (kglb[k][j] * kglb[i][j]);
        }
        kglb[k][i] = sum / diag;
        kglb[i][k] = 0.0;
      }
    }
  }
  i = ndof - 1;
  diag = kglb[i][i];
  for (int j = 0; j < i; j++) {
    diag -= kglb[j][j] * (kglb[i][j] * kglb[i][j]);
  }
  kglb[i][i] = diag;
  return ierr;
}

// прямая подстановка. Решение L*X=R
void direktLDLT(int ndof, double **kglb, double *x,
                double *r) { // r - массив нагрузок; x - рабочий массив
  // разложение симметричной матрицы kglb в LDLT
  double sum;
  x[0] = r[0];
  for (int i = 1; i < ndof; i++) {
    sum = 0.0;
#pragma omp parallel shared(kglb, i) reduction(+ : sum)
    {
#pragma omp for
      for (int j = 0; j < i; j++) {
        sum += kglb[i][j] * x[j];
      }
    }
    x[i] = r[i] - sum;
  }
}

// Диагональное решение. Решение D*Y=X
bool diagLDLT(int ndof, double **kglb, double *x) {
  bool ierr = false;
  double diag;
#pragma omp parallel shared(kglb, ndof)
  {
#pragma omp for
    for (int i = 0; i < ndof; i++) {
      diag = kglb[i][i];
      if (fabs(diag) < 1.0e-20) {
        printf("Решение невозможно. Нулевой диагональный элемент - %d ", i);
        ierr = true;
      }
      x[i] /= diag;
    }
  }
  return ierr;
}

// обратная подстановка. Решение LT*U=X
void rechLDLT(int ndof, double **kglb, double *u, double *x) {
  // u - массив перемещений
  double sum;
  u[ndof - 1] = x[ndof - 1];
  for (int i = ndof - 2; i >= 0; i--) {
    sum = 0.0;
#pragma omp parallel shared(kglb, ndof, i) reduction(+ : sum)
    {
#pragma omp for
      for (int j = i + 1; j < ndof; j++) {
        sum += kglb[j][i] * u[j];
      }
    }
    u[i] = x[i] - sum;
  }
}

bool solveLinearSystemLDLT(double **kglb, double *u, double *r, double *x,
                           int ndof) {
  bool ierr = matrLDLT(ndof, kglb);
  if (!ierr) {
    direktLDLT(ndof, kglb, x, r);
    ierr = diagLDLT(ndof, kglb, x);
    if (!ierr) {
      rechLDLT(ndof, kglb, u, x);
    }
  }
  return ierr;
}

// расчет напряжений и деформаций в элементах модели
void stressModel(int ndofysla, int nys, int nelem, int **jt03, double **car,
                 double h, double e, double puas, double *u, double **strain,
                 double **stress) {
  // формирование матрицы деформаций deformMtrx[3][6]
  double *dataDeformMtrx = (double *)malloc(3 * 6 * sizeof(double));
  double **deformMtrx = (double **)malloc(3 * sizeof(double *));
  for (int i = 0; i < 3; i++) {
    deformMtrx[i] = dataDeformMtrx + i * 3;
  }
  // формирование матрицы напряжений strsMtrx[3][6]
  double *dataStrsMtrx = (double *)malloc(3 * 6 * sizeof(double));
  double **strsMtrx = (double **)malloc(3 * sizeof(double *));
  for (int i = 0; i < 3; i++) {
    strsMtrx[i] = dataStrsMtrx + i * 3;
  }
  // перемещения(uElement[6]), напряжения(eStress[3]) и деформации(eStrain[3])
  // одного элемента
  double *uElement = (double *)malloc(6 * sizeof(double));
  double *eStress = (double *)malloc(3 * sizeof(double));
  double *eStrain = (double *)malloc(3 * sizeof(double));
#pragma omp parallel shared(jt03, car, uElement, eStress, eStrain)
  {
#pragma omp for
    for (int ielem = 0; ielem < nelem; ielem++) {
      nodeNumber node = {jt03[0][ielem] - 1, jt03[1][ielem] - 1,
                         jt03[2][ielem] - 1};
      coord coord1 = {car[0][node.iys1], car[1][node.iys1]};
      coord coord2 = {car[0][node.iys2], car[1][node.iys2]};
      coord coord3 = {car[0][node.iys3], car[1][node.iys3]};
      stressPlanElem(coord1, coord2, coord3, h, e, puas, deformMtrx, strsMtrx);
      uElement[0] = u[node.iys1 * ndofysla];
      uElement[1] = u[node.iys1 * ndofysla + 1];
      uElement[2] = u[node.iys2 * ndofysla];
      uElement[3] = u[node.iys2 * ndofysla + 1];
      uElement[4] = u[node.iys3 * ndofysla];
      uElement[5] = u[node.iys3 * ndofysla + 1];
      double sum;
      for (int i = 0; i < 3; i++) {
        sum = 0.;
        for (int j = 0; j < 6; j++) {
          sum += strsMtrx[i][j] * uElement[j];
        }
        eStress[i] = sum;
        sum = 0.;
        for (int j = 0; j < 6; j++) {
          sum += deformMtrx[i][j] * uElement[j];
        }
        eStrain[i] = sum;
      }
      deformComp dC = {eStrain[0], eStrain[1], eStrain[2]};
      strain[0][ielem] = dC.ex;
      strain[1][ielem] = dC.ey;
      strain[2][ielem] = dC.exy;
      strain[3][ielem] =
          sqrt(2.0) / 3.0 *
          sqrt(dC.ex * dC.ey + (dC.ex - dC.ey) * (dC.ex - dC.ey) +
               dC.ey * dC.ey + 1.5 * dC.exy * dC.exy);
      stressComp sC = {eStress[0], eStress[1], eStress[2]};
      stress[0][ielem] = sC.sx;
      stress[1][ielem] = sC.sy;
      stress[2][ielem] = sC.sxy;
      stress[3][ielem] = sqrt(sC.sx * sC.sx - sC.sx * sC.sy + sC.sy * sC.sy +
                              3.0 * sC.sxy * sC.sxy);
    }
  }
  free_memory(7, dataDeformMtrx, dataStrsMtrx, strsMtrx, deformMtrx, uElement,
              eStress, eStrain);
}

void free_memory(int count, ...) {
  va_list args;
  va_start(args, count);
  for (int i = 0; i < count; i++) {
    void *pnt = va_arg(args, void *);
    free(pnt);
  }
  va_end(args);
}
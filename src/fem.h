#ifndef FEM_H
#define FEM_H
#include <stdbool.h>

typedef struct {
  double x;
  double y;
} coord;

typedef struct {
  int iys1;
  int iys2;
  int iys3;
} nodeNumber;

typedef struct {
  double ex;
  double ey;
  double exy;
} deformComp;

typedef struct {
  double sx;
  double sy;
  double sxy;
} stressComp;

bool matrLDLT(int ndof, double **kglb);
void direktLDLT(int ndof, double **kglb, double *x, double *r);
bool diagLDLT(int ndof, double **kglb, double *x);
void rechLDLT(int ndof, double **kglb, double *u, double *x);
bool solveLinearSystemLDLT(double **kglb, double *u, double *r, double *x,
                           int ndof);
void stressModel(int ndofysla, int nelem, int **jt03, double **car, double e,
                 double puas, double *u, double **strain, double **stress);
void AssembleLocalStiffnessToGlobal(double **kglb, int **jt03, double **car,
                                    int nelem, double e, double h, double puas,
                                    int ndofysla);
void free_memory(int, ...);

#endif
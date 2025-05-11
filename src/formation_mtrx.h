#ifndef FORMATION_MTRX
#define FORMATION_MTRX

#include "fem.h"

double **formationDeformMtrx(double **deformMtrx, coord coord1, coord coord2,
                             coord coord3, double a2);
double **formationElastMtrx(double **elastMtrx, double e, double puas);
void planeElement(coord coord1, coord coord2, coord coord3, double e, double h,
                  double puas, double **gest);
void assemblyGlobMatr(int ndofysla, nodeNumber node, double **gest,
                      double **kglb);
void stressPlanElem(coord coord1, coord coord2, coord coord3, double e,
                    double puas, double **deformMtrx, double **strsMatr);
void FillConstrainedLoadedNodes(int **nodePres, int *lenNodePres,
                                int **nodeZakrU, int *lenNodeZakrU,
                                int **nodeZakrV, int *lenNodeZakrV,
                                double **car, int nys);
void MakeConstrained(int *nodeZakr, int lenNodeZakr, double **kglb,
                     int ndofysla);
void SetLoadVector(double *r, int lenNodePres, int *nodePres, int ndofysla,
                   int ndof, float load);
void makeDoubleMtrx(double **dataMtrx, double ***mtrx, int row, int col);
void makeIntegerMtrx(int **dataMtrx, int ***mtrx, int row, int col);
short readFromFile(char *filename, int *nys, double **dataCar, double ***car,
                   int *nelem, int **data_jt03, int ***jt03);
void writeResult(char *filename, int **jt03, double **strain, double **stress,
                 double *r, double *u, int nelem, int nys, int ndof);

#endif
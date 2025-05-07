#ifndef TEST_FEM_H
#define TEST_FEM_H

#include <check.h>
#include <stdlib.h>

#include "formation_mtrx.h"
#include "fem.h"

coord coord1 = {0., 0.};
coord coord2 = {20., 0.};
coord coord3 = {20., 20.};
double e = 2.1e5;
double puas = 0.3;
double h = 1.;
int ndofysla = 2;
double car[3][4] = {
    {0., 1000., 1000., 0.},
    {0., 0., 200., 200.},
    {0., 0., 0., 0.}
};
int jt03[3][2] = {
    {1, 1},
    {2, 3},
    {3, 4}
};
int nelem = 2;
int ndof = 8;

Suite *test_formatationDeformMtrx(void);
Suite *test_formatationElastMtrx(void);
Suite *test_formatationGestMtrx(void);
Suite *test_formatationStressMtrx(void);
Suite *test_formatationGlobalMtrx(void);

#endif
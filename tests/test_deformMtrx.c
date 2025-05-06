#include "test_fem.h"

START_TEST(test_formatationDefomrmMtrx) {
    double *dataDeformMtrx = (double *)malloc(3 * 6 * sizeof(double));
    double **deformMxtr = (double **)malloc(3 * sizeof(double *));
    for (int i = 0; i < 3; i++) {
        deformMxtr[i] = dataDeformMtrx + i * 6;
    }
    coord coord1 = {0., 0.};
    coord coord2 = {20., 0.};
    coord coord3 = {20., 20.};
    double a2 = coord2.x * coord3.y - coord2.y * coord3.x +
                coord3.x * coord1.y - coord3.y * coord1.x +
                coord1.x * coord2.y - coord1.y * coord2.x;
    deformMxtr = formationDeformMtrx(deformMxtr, coord1, coord2, coord3, a2);
    double resMtrx[3][6] = {
        {a2 * (-20.), 0., a2 * 20., 0., 0., 0.,},
        {0., 0., 0., a2 * (-20.), 0., a2 * 20.},
        {0., a2 * (-20.), a2 * (-20.), a2 * 20., a2 * 20., 0.}
    };
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 6; col++) {
            ck_assert_double_eq(deformMxtr[row][col], resMtrx[row][col]);
        }
    }
    free_memory(2, dataDeformMtrx, deformMxtr);
} END_TEST

START_TEST(test_formatationElastMtrx) {
    double *dataElastMtrx = (double *)malloc(3 * 3 * sizeof(double));
    double **elastMtrx = (double **)malloc(3 * sizeof(double *));
    for (int i = 0; i < 3; i++) {
        elastMtrx[i] = dataElastMtrx + i * 3;
    }
    double e = 2.1e5;
    double puas = 0.3;
    elastMtrx = formationElastMtrx(elastMtrx, e, puas);
    // дописать тест!!
}
#include "test_fem.h"

START_TEST(test_deformMtrx_formatation) {
    coord coord1 = {0., 0.};
    coord coord2 = {20., 0.};
    coord coord3 = {20., 20.};

    double *dataDeformMtrx = (double *)malloc(3 * 6 * sizeof(double));
    double **deformMxtr = (double **)malloc(3 * sizeof(double *));
    for (int i = 0; i < 3; i++) {
        deformMxtr[i] = dataDeformMtrx + i * 6;
    }
    
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

Suite *test_formatationDeformMtrx(void) {
    Suite *s = suite_create("formatationDeformMtrx");
    TCase *tc = tcase_create(" test_deformMtrx ");

    suite_add_tcase(s, tc);
    tcase_add_test(tc, test_deformMtrx_formatation);

    suite_add_tcase(s, tc);
    return s;
}
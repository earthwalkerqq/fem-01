#include "test_fem.h"

START_TEST(test_planeElem_formatation) {
    coord coord1 = {0., 0.};
    coord coord2 = {20., 0.};
    coord coord3 = {20., 20.};
    double e = 2.1e5;
    double puas = 0.3;
    double h = 1.;

    double *dataGEST = (double *)malloc(6 * 6 * sizeof(double));
    double **gest = NULL;
    makeDoubleMtrx(&dataGEST, &gest, 6, 6);
    planeElement(coord1, coord2, coord3, e, h, puas, gest);
    double *dataDeformMtrx = (double *)malloc(3 * 6 * sizeof(double));
    double **deformMtrx = NULL;
    makeDoubleMtrx(&dataDeformMtrx, &deformMtrx, 3, 6);
    double *dataStrsMtrx = (double *)malloc(3 * 6 * sizeof(double));
    double **strsMtrx = NULL;
    makeDoubleMtrx(&dataStrsMtrx, &strsMtrx, 3, 6);
    stressPlanElem(coord1, coord2, coord3, e, puas, deformMtrx, strsMtrx);
    // double a2 = coord2.x * coord3.y - coord3.x * coord2.y + coord3.x * coord1.y -
    //             coord1.x * coord3.y + coord1.x * coord2.y - coord2.x * coord1.y;
    // double vol = h * a2 * 0.5;
    double resMtrx[6][6] = {
        {1.154e5, 0., -1.154e5, 3.462e4, 0., -3.462e4},
        {0., 4.038e4, 4.038e4, -4.038e44, -4.038e4, 0.},
        {-1.154e5, 4.038e4, 1.558e5, -7.5e4, -4.038e4, 3.462e4},
        {3.462e44, -4.038e4, -7.5e4, 1.558e5, 4.038e44, -1.154e5},
        {0., -4.038e4, -4.038e4, 4.038e4, 4.038e4, 0.},
        {-3.462e4, 0., 3.462e4, -1.154e5, 0., 1.154e5}
    };
    // double *dataResMtrx = (double *)malloc(6 * 6 * sizeof(double));
    // double **resMtrx = makeDoubleMtrx(dataResMtrx, 6, 6);
    // double sum;
    // for (int i = 0; i < 6; i++) {
    //     for (int j = 0; j < 6; j++) {
    //         sum = 0.0;
    //         for (int k = 0; k < 3; k++) {
    //             sum += strsMtrx[k][i] * deformMtrx[k][j];
    //         }
    //         resMtrx[i][j] = vol * sum;
    //     }
    // }
    for (int row = 0; row < 6; row++) {
        for (int col = 0; col < 6; col++) {
            ck_assert_double_eq(resMtrx[row][col], gest[row][col]);
        }
    }
    free_memory(8, dataDeformMtrx, dataGEST, dataStrsMtrx,
                deformMtrx, gest, strsMtrx);
} END_TEST

Suite *test_formatationPlaneElem(void) {
    Suite *s = suite_create("formatationPlaneElem");
    TCase *tc = tcase_create(" test_planeElem ");

    suite_add_tcase(s, tc);
    tcase_add_test(tc, test_planeElem_formatation);

    suite_add_tcase(s, tc);
    return s;
}
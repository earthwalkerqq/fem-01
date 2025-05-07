#include "test_fem.h"

START_TEST(test_formatationGlobalMtrx_1) {
    double *dataGEST = (double *)malloc(6 * 6 * sizeof(double));
    double **gest = makeDoubleMtrx(&dataGEST, 6, 6);
    double *dataKGLB = (double *)malloc(ndof * ndof * sizeof(double));
    double **kglb = makeDoubleMtrx(dataKGLB, ndof, ndof);
    for (int ielem = 0; ielem < nelem; ielem++) {
        nodeNumber node = {
            jt03[0][ielem],
            jt03[1][ielem],
            jt03[2][ielem]
        };
        planeElement(coord1, coord2, coord3, e, h, puas, gest);
        assemblyGlobMatr(ndofysla, node, gest, kglb);
    }
    double resMtrx[8][8] = {
        {}
    }
} END_TEST

START_TEST(test_formatationGlobalMtrx_2) {
    double *dataGEST = (double *)malloc(6 * 6 * sizeof(double));
    double **gest = makeDoubleMtrx(&dataGEST, 6, 6);
    double *dataKGLB = (double *)malloc(6 * 6 * sizeof(double));
    double **kglb = makeDoubleMtrx(&dataKGLB, 6, 6);
    nodeNumber node = {jt03[0][0], jt03[1][0], jt03[2][0]};
    planeElement(coord1, coord2, coord3, e, h, puas, gest);
    assemblyGlobMatr(ndofysla, node, gest, kglb);
    for (int row = 0; row < 6; row++) {
        for (int col = 0; col < 6; col++) {
            ck_assert_double_eq(gest[row][col], kglb[row][col]);
        }
    }
    free_memory(4, dataGEST, dataKGLB, gest, kglb)
} END_TEST

Suite *test_formatationGlobalMtrx(void) {
    Suite *s = suite_create("formatationGlobalMtrx");
    TCase *tc = tcase_create(" test_globalMtrx ");

    suite_add_tcase(s, tc);
    tcase_add_test(tc, test_formatationGlobalMtrx_1);
    tcase_add_test(tc, test_formatationGlobalMtrx_2);

    suite_add_tcase(s, tc);
    return s;
}
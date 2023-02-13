#include <stdio.h>
#include <memory.h>
#include <stdlib.h>

#include "lpm.h"
#include "utilities.h"

double *createRMatrix(int eulerflag, double angles[])
{
    /* Below notations are referred to wiki: https://en.wikipedia.org/wiki/Euler_angles */
    /* 1, 2, 3 represent the rotation angles angle1, angle2, gamma */
    /* x, y, z represent fixed frames */
    /* X, Y, Z mean the rotation matrices that give elemental rotations about axies x, y, z */
    /* e.g. X1 means a rotation about x axis by an angle which is angle1 */

    /* Examples of elemental rotation matrices, c1 is cos_alpha, s2 is sin_beta, etc. */
    // X1 = [1   0   0 ]  Y2 = [c2    0   s2]  Z3 = [c3  -s3   0]
    //      [0  c1  -s1]       [ 0    1    0]       [s3   c3   0]
    //      [0  s1   c1]       [-s2   0   c2]       [ 0    0   1]

    // In 2D, only Z3 is used, therefore R = I-I-Z3 = Z3
    // R = [c3  -s3   0]  or R = [c3  -s3]
    //     [s3   c3   0]         [s3   c3]
    //     [0     0   1]

    double angle1 = angles[0], angle2 = angles[1], angle3 = angles[2];
    static double R_matrix[NDIM*NDIM];

    if (eulerflag == 0)
    {
        /* R = X1-Y2-Z3; row major */
        // X1-Y2-Z3 = [c2c3          -c2s3          s2  ]
        //            [s1s2c3+c1s3   -s1s2s3+c1c3  -s1c2]
        //            [-c1s2c3+s1s3   c1s2s3+s1c3   c1c2]

        R_matrix[0] = cos(angle2) * cos(angle3);
        R_matrix[1] = -cos(angle2) * sin(angle3);
        R_matrix[2] = sin(angle2);
        R_matrix[3] = cos(angle1) * sin(angle3) + sin(angle1) * sin(angle2) * cos(angle3);
        R_matrix[4] = cos(angle1) * cos(angle3) - sin(angle1) * sin(angle2) * sin(angle3);
        R_matrix[5] = -sin(angle1) * cos(angle2);
        R_matrix[6] = sin(angle1) * sin(angle3) - cos(angle1) * sin(angle2) * cos(angle3);
        R_matrix[7] = sin(angle1) * cos(angle3) + cos(angle1) * sin(angle2) * sin(angle3);
        R_matrix[8] = cos(angle1) * cos(angle2);
    }
    else if (eulerflag == 1)
    {
        /* R = Z1-Y2-Z3, using Kocks convention, angle3 = pi - angle3' (angle3' is Kocks input); row major */
        // Z1-Y2-Z3 = [c1c2c3-s1s3  -c1c2s3-s1c3  s2c1]
        //            [s1c2c3+c1s3  -s1c2s3+c1c3  s2s1]
        //            [-s2c3         s2s3         c2  ]

        R_matrix[0] = -cos(angle1) * cos(angle2) * cos(angle3) - sin(angle1) * sin(angle3);
        R_matrix[1] = cos(angle3) * sin(angle1) - sin(angle3) * cos(angle2) * cos(angle1);
        R_matrix[2] = sin(angle2) * cos(angle1);
        R_matrix[3] = sin(angle3) * cos(angle1) - cos(angle3) * cos(angle2) * sin(angle1);
        R_matrix[4] = -cos(angle1) * cos(angle3) - sin(angle1) * cos(angle2) * sin(angle3);
        R_matrix[5] = sin(angle2) * sin(angle1);
        R_matrix[6] = cos(angle3) * sin(angle2);
        R_matrix[7] = sin(angle3) * sin(angle2);
        R_matrix[8] = cos(angle2);
    }
    else if (eulerflag == 2)
    {
        /* R = Z1-X2-Z3, using Bunge convention; row major */
        // Z1-X2-Z3 = [c1c3-s1c2s3  -c1s3-s1c2c3  s1s2]
        //            [c1c2s3+c3s1  c1c2c3-s1s3  -c1s2]
        //            [s2s3          s2c3         c2  ]

        R_matrix[0] = -sin(angle1) * cos(angle2) * sin(angle3) + cos(angle1) * cos(angle3);
        R_matrix[1] = -sin(angle3) * cos(angle1) - cos(angle3) * cos(angle2) * sin(angle1);
        R_matrix[2] = sin(angle1) * sin(angle2);

        R_matrix[3] = cos(angle3) * sin(angle1) + sin(angle3) * cos(angle2) * cos(angle1);
        R_matrix[4] = -sin(angle1) * sin(angle3) + cos(angle1) * cos(angle2) * cos(angle3);
        R_matrix[5] = -cos(angle1) * sin(angle2);

        R_matrix[6] = sin(angle2) * sin(angle3);
        R_matrix[7] = sin(angle2) * cos(angle3);
        R_matrix[8] = cos(angle2);
    }

    return R_matrix;
}

std::vector<std::array<double, NDIM>> createCuboidSC3D(double box[], UnitCell cell, double R_matrix[])
{
    /* settings for matrix-vector product, BLAS */
    CBLAS_LAYOUT layout = CblasRowMajor;
    CBLAS_TRANSPOSE trans = CblasNoTrans;

    int lda = 3, incx = 1, incy = 1;
    double blasAlpha = 1.0, blasBeta = 0.0;

    /* initialize the particle xyz_t (a larger system) */
    double a = sqrt(pow(box[1] - box[0], 2) + pow(box[3] - box[2], 2) + pow(box[5] - box[4], 2));
    double box_t[6] = {-a, a, -a, a, -a, a};

    /* model parameters */
    double hx = 2 * cell.radius;
    double hy = hx;
    double hz = hx;
    int particles_first_row = 1 + (int)floor((box_t[1] - box_t[0]) / hx);
    int rows = 1 + (int)floor((box_t[3] - box_t[2]) / hy);
    int layers = 1 + (int)floor((box_t[5] - box_t[4]) / hz);

    std::vector<std::array<double, NDIM>> xyz_t;
    double p[NDIM] = {0}, p_new[NDIM] = {0};
    for (int k = 1; k <= layers; k++)
    {
        p[2] = box_t[4] + hz * (k - 1); // z
        for (int j = 1; j <= rows; j++)
        {
            p[1] = box_t[2] + hy * (j - 1); // y
            for (int i = 1; i <= particles_first_row; i++)
            {
                p[0] = box[0] + hx * (i - 1); // x
                cblas_dgemv(layout, trans, NDIM, NDIM, blasAlpha, R_matrix, lda, p, incx, blasBeta, p_new, incy);

                /* test if the rotated system is within the specified domain */
                if (p_new[0] >= box[0] && p_new[0] <= box[1] &&
                    p_new[1] >= box[2] && p_new[1] <= box[3] &&
                    p_new[2] >= box[4] && p_new[2] <= box[5])
                {
                    std::array<double, NDIM> p_arr{p_new[0], p_new[1], p_new[2]};
                    xyz_t.push_back(p_arr);
                }
            }
        }
    }

    return xyz_t;
}
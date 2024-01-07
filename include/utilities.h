#pragma once
#ifndef UTILITIES_H
#define UTILITIES_H

#include "unit_cell.h"
#include "lpm.h"

MatrixXd createRMatrix(int eulerflag, double angles[])
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
    MatrixXd R_matrix(NDIM, NDIM);

    if (eulerflag == 0)
    {
        /* R = X1-Y2-Z3; row major */
        // X1-Y2-Z3 = [c2c3          -c2s3          s2  ]
        //            [s1s2c3+c1s3   -s1s2s3+c1c3  -s1c2]
        //            [-c1s2c3+s1s3   c1s2s3+s1c3   c1c2]

        R_matrix(0, 0) = cos(angle2) * cos(angle3);
        R_matrix(0, 1) = -cos(angle2) * sin(angle3);
        R_matrix(0, 2) = sin(angle2);
        R_matrix(1, 0) = cos(angle1) * sin(angle3) + sin(angle1) * sin(angle2) * cos(angle3);
        R_matrix(1, 1) = cos(angle1) * cos(angle3) - sin(angle1) * sin(angle2) * sin(angle3);
        R_matrix(1, 2) = -sin(angle1) * cos(angle2);
        R_matrix(2, 0) = sin(angle1) * sin(angle3) - cos(angle1) * sin(angle2) * cos(angle3);
        R_matrix(2, 1) = sin(angle1) * cos(angle3) + cos(angle1) * sin(angle2) * sin(angle3);
        R_matrix(2, 2) = cos(angle1) * cos(angle2);
    }
    else if (eulerflag == 1)
    {
        /* R = Z1-Y2-Z3, using Kocks convention, angle3 = pi - angle3' (angle3' is Kocks input); row major */
        // Z1-Y2-Z3 = [c1c2c3-s1s3  -c1c2s3-s1c3  s2c1]
        //            [s1c2c3+c1s3  -s1c2s3+c1c3  s2s1]
        //            [-s2c3         s2s3         c2  ]

        R_matrix(0, 0) = -cos(angle1) * cos(angle2) * cos(angle3) - sin(angle1) * sin(angle3);
        R_matrix(0, 1) = cos(angle3) * sin(angle1) - sin(angle3) * cos(angle2) * cos(angle1);
        R_matrix(0, 2) = sin(angle2) * cos(angle1);
        R_matrix(1, 0) = sin(angle3) * cos(angle1) - cos(angle3) * cos(angle2) * sin(angle1);
        R_matrix(1, 1) = -cos(angle1) * cos(angle3) - sin(angle1) * cos(angle2) * sin(angle3);
        R_matrix(1, 2) = sin(angle2) * sin(angle1);
        R_matrix(2, 0) = cos(angle3) * sin(angle2);
        R_matrix(2, 1) = sin(angle3) * sin(angle2);
        R_matrix(2, 2) = cos(angle2);
    }
    else if (eulerflag == 2)
    {
        /* R = Z1-X2-Z3, using Bunge convention; row major */
        // Z1-X2-Z3 = [c1c3-s1c2s3  -c1s3-s1c2c3  s1s2]
        //            [c1c2s3+c3s1  c1c2c3-s1s3  -c1s2]
        //            [s2s3          s2c3         c2  ]

        R_matrix(0, 0) = -sin(angle1) * cos(angle2) * sin(angle3) + cos(angle1) * cos(angle3);
        R_matrix(0, 1) = -sin(angle3) * cos(angle1) - cos(angle3) * cos(angle2) * sin(angle1);
        R_matrix(0, 2) = sin(angle1) * sin(angle2);
        R_matrix(1, 0) = cos(angle3) * sin(angle1) + sin(angle3) * cos(angle2) * cos(angle1);
        R_matrix(1, 1) = -sin(angle1) * sin(angle3) + cos(angle1) * cos(angle2) * cos(angle3);
        R_matrix(1, 2) = -cos(angle1) * sin(angle2);
        R_matrix(2, 0) = sin(angle2) * sin(angle3);
        R_matrix(2, 1) = sin(angle2) * cos(angle3);
        R_matrix(2, 2) = cos(angle2);
    }

    return R_matrix;
}

std::vector<std::array<double, NDIM>> createPlateSQ2D(std::array<double, 2 * NDIM> box, UnitCell cell, MatrixXd R_matrix)
{
    /* 2D square lattice (double layer neighbor) */
    double a = sqrt(pow(box[1] - box[0], 2) + pow(box[3] - box[2], 2));
    double box_t[6] = {-a, a, -a, a, -a, a};

    /* model parameters */
    double hx = 2 * cell.radius;
    double hy = hx;
    int particles_first_row = floor((box_t[1] - box_t[0]) / hx);
    int rows = floor((box_t[3] - box_t[2]) / hy);

    std::vector<std::array<double, NDIM>> xyz_t; /* Coordinate vector */
    VectorXd p(NDIM), p_new(NDIM);
    for (int j = 1; j <= rows; j++)
    {
        p(1) = box_t[2] + hy * (j - 1) + cell.radius;
        for (int i = 1; i <= particles_first_row; i++)
        {
            p(0) = box_t[0] + hx * (i - 1) + cell.radius;
            p_new = R_matrix * p;

            if (p_new(0) >= box[0] && p_new(0) <= box[1] &&
                p_new(1) >= box[2] && p_new(1) <= box[3])
            {
                std::array<double, NDIM> p_arr{p_new(0), p_new(1), 0};
                xyz_t.push_back(p_arr);
            }
        }
    }

    return xyz_t;
}

std::vector<std::array<double, NDIM>> createPlateHEX2D(std::array<double, 2 * NDIM> box, UnitCell cell, MatrixXd R_matrix)
{
    /* 2D hexagon lattice (double layer neighbor) */
    double a = sqrt(pow(box[1] - box[0], 2) + pow(box[3] - box[2], 2));
    double box_t[6] = {-a, a, -a, a, -a, a};
    /* model parameters */
    double hx = 2 * cell.radius;
    double hy = hx * sqrt(3.0) / 2.0;
    int particles_first_row = 1 + floor((box_t[1] - box_t[0]) / hx);
    int rows = 1 + floor((box_t[3] - box_t[2]) / hy);

    std::vector<std::array<double, NDIM>> xyz_t; /* Coordinate vector */
    VectorXd p(NDIM), p_new(NDIM);
    std::vector<std::array<double, NDIM>> xyz;
    for (int j = 1; j <= rows; j++)
    {
        p(1) = box_t[2] + hy * (j - 1);
        if (j % 2 == 1)
            for (int i = 1; i <= particles_first_row; i++)
            {
                p(0) = box_t[0] + hx * (i - 1);
                xyz.push_back(std::array<double, NDIM>{p(0), p(1), 0});
            }
        else
            for (int i = 1; i <= particles_first_row - 1; i++)
            {
                p(0) = box_t[0] + hx * (2 * i - 1) / 2.0;
                xyz.push_back(std::array<double, NDIM>{p(0), p(1), 0});
            }
    }

    /* rotate the system */
    for (std::array<double, NDIM> pt : xyz)
    {
        p(0) = pt[0];
        p(1) = pt[1];
        p(2) = pt[2];

        p_new = R_matrix * p;

        /* test if the rotated system is within the specified domain */
        if (p_new(0) >= box[0] && p_new(0) <= box[1] &&
            p_new(1) >= box[2] && p_new(1) <= box[3])
        {
            std::array<double, NDIM> p_arr{p_new(0), p_new(1), p_new(2)};
            xyz_t.push_back(p_arr);
        }
    }

    return xyz_t;
}

std::vector<std::array<double, NDIM>> createCuboidSC3D(std::array<double, 2 * NDIM> box, UnitCell cell, MatrixXd R_matrix)
{

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
    VectorXd p(NDIM), p_new(NDIM);
    for (int k = 1; k <= layers; k++)
    {
        p(2) = box_t[4] + hz * (k - 1); // z
        for (int j = 1; j <= rows; j++)
        {
            p(1) = box_t[2] + hy * (j - 1); // y
            for (int i = 1; i <= particles_first_row; i++)
            {
                p(0) = box[0] + hx * (i - 1); // x
                p_new = R_matrix * p;

                /* test if the rotated system is within the specified domain */
                if (p_new(0) >= box[0] && p_new(0) <= box[1] &&
                    p_new(1) >= box[2] && p_new(1) <= box[3] &&
                    p_new(2) >= box[4] && p_new(2) <= box[5])
                {
                    std::array<double, NDIM> p_arr{p_new(0), p_new(1), p_new(2)};
                    xyz_t.push_back(p_arr);
                }
            }
        }
    }

    return xyz_t;
}

std::vector<std::array<double, NDIM>> createCuboidFCC3D(std::array<double, 2 * NDIM> box, UnitCell cell, MatrixXd R_matrix)
{
    /* initialize the particle xyz_t (a larger system) */
    double a = sqrt(pow(box[1] - box[0], 2) + pow(box[3] - box[2], 2) + pow(box[5] - box[4], 2));
    double box_t[6] = {-a, a, -a, a, -a, a};

    /* model parameters */
    double hx = 4 * cell.radius / sqrt(2.0);
    double hy = hx;
    double hz = hx;
    int particles_first_row = 1 + (int)floor((box_t[1] - box_t[0]) / hx);
    int rows = 1 + (int)floor((box_t[3] - box_t[2]) / (hy / 2.0));
    int layers = 1 + (int)floor((box_t[5] - box_t[4]) / (hz / 2.0));

    std::vector<std::array<double, NDIM>> xyz;
    VectorXd p(NDIM), p_new(NDIM);

    for (int k = 1; k <= layers; k++)
    {
        p(2) = box_t[4] + hz / 2.0 * (k - 1);
        if (k % 2 == 1)
        {
            for (int j = 1; j <= rows; j++)
            {
                p(1) = box_t[2] + hy / 2.0 * (j - 1);
                if (j % 2 == 1)
                {
                    for (int i = 1; i <= particles_first_row; i++)
                    {
                        p(0) = box_t[0] + hx * (i - 1);
                        xyz.push_back(std::array<double, NDIM>{p(0), p(1), p(2)});
                    }
                }
                else
                {
                    for (int i = 1; i <= particles_first_row - 1; i++)
                    {
                        p(0) = box_t[0] + hx * (i - 1) + hx / 2.0;
                        xyz.push_back(std::array<double, NDIM>{p(0), p(1), p(2)});
                    }
                }
            }
        }
        else
        {
            for (int j = 1; j <= rows; j++)
            {
                p(1) = box_t[2] + hy / 2.0 * (j - 1);
                if (j % 2 == 1)
                {
                    for (int i = 1; i <= particles_first_row - 1; i++)
                    {
                        p(0) = box_t[0] + hx * (i - 1) + hx / 2.0;
                        xyz.push_back(std::array<double, NDIM>{p(0), p(1), p(2)});
                    }
                }
                else
                {
                    for (int i = 1; i <= particles_first_row; i++)
                    {
                        p(0) = box_t[0] + hx * (i - 1);
                        xyz.push_back(std::array<double, NDIM>{p(0), p(1), p(2)});
                    }
                }
            }
        }
    }

    std::vector<std::array<double, NDIM>> xyz_t;
    for (std::array<double, NDIM> pt : xyz)
    {
        p(0) = pt[0];
        p(1) = pt[1];
        p(2) = pt[2];

        p_new = R_matrix * p;

        /* test if the rotated system is within the specified domain */
        if (p_new(0) >= box[0] && p_new(0) <= box[1] &&
            p_new(1) >= box[2] && p_new(1) <= box[3] &&
            p_new(2) >= box[4] && p_new(2) <= box[5])
        {
            std::array<double, NDIM> p_arr{p_new(0), p_new(1), p_new(2)};
            xyz_t.push_back(p_arr);
        }
    }

    return xyz_t;
}

#endif

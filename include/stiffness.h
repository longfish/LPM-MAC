#pragma once
#ifndef STIFFNESS_H
#define STIFFNESS_H

#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include "lpm.h"
#include "unit_cell.h"
#include "bond.h"
#include "derivative.h"

template <int nlayer>
class Particle;

template <int nlayer>
class Stiffness
{
    StiffnessMode mode; // use finite difference or analytical approach to compute local stiffness

public:
    Eigen::SparseMatrix<double> K_matrix;
    Eigen::VectorXd residual;

    void reset(std::vector<Particle<nlayer> *> &pt_sys);
    void calcStiffness(std::vector<Particle<nlayer> *> &pt_sys);
    void updateStiffnessDispBC(std::vector<Particle<nlayer> *> &pt_sys);

    std::array<std::array<double, NDIM>, NDIM> localStiffness(Particle<nlayer> *pi, Particle<nlayer> *pj);
    std::array<std::array<double, NDIM>, NDIM> localStiffnessFD(Particle<nlayer> *pi, Particle<nlayer> *pj);
    std::array<std::array<double, NDIM>, NDIM> localStiffnessANA(Particle<nlayer> *pi, Particle<nlayer> *pj);

    Stiffness(std::vector<Particle<nlayer> *> &pt_sys, StiffnessMode p_mode)
        : K_matrix(pt_sys.size() * pt_sys[0]->cell.dim, pt_sys.size() * pt_sys[0]->cell.dim),
          residual(pt_sys.size() * pt_sys[0]->cell.dim)
    { // Given a particle system, construct a stiffness matrix and solver
        mode = p_mode;
    }
};

template <int nlayer>
void Stiffness<nlayer>::reset(std::vector<Particle<nlayer> *> &pt_sys)
{
    // reset the stiffness matrix to be zero
    K_matrix.setZero();
}

template <int nlayer>
std::array<std::array<double, NDIM>, NDIM> Stiffness<nlayer>::localStiffness(Particle<nlayer> *pi, Particle<nlayer> *pj)
{
    if (mode == StiffnessMode::Analytical)
        return localStiffnessANA(pi, pj);
    return localStiffnessFD(pi, pj);
}

std::array<std::array<double, NDIM>, NDIM> sumArray(std::array<std::array<double, NDIM>, NDIM> &a, std::array<std::array<double, NDIM>, NDIM> &b)
{
    std::array<std::array<double, NDIM>, NDIM> local;
    for (int i = 0; i < NDIM; i++)
        for (int j = 0; j < NDIM; j++)
            local[i][j] = a[i][j] + b[i][j];
    return local;
}

template <int nlayer>
std::array<std::array<double, NDIM>, NDIM> Stiffness<nlayer>::localStiffnessANA(Particle<nlayer> *pi, Particle<nlayer> *pj)
{
    bool flag1{false}, flag2{false};
    std::array<std::array<double, NDIM>, NDIM> K_local{0}, K_temp1{0}, K_temp2{0};
    if (pi->hasAFEMneighbor(pj, 0))
    {
        K_temp1 = fdu2dxyz1(pi, pj);
        flag1 = true;
    }
    if (pi->hasAFEMneighbor(pj, 1))
    {
        K_temp2 = fdu2dxyz2(pi, pj);
        flag2 = true;
    }

    if (pi->id == pj->id)
    {
        K_local = fdu2dxyz(pi);
        // if (pi->id == 40)
        // {
        //     printf("%f, %f, %f\n%f, %f, %f\n%f, %f, %f\n\n", K_local[0][0], K_local[0][1], K_local[0][2],
        //            K_local[1][0], K_local[1][1], K_local[1][2],
        //            K_local[2][0], K_local[2][1], K_local[2][2]);
        // }
    }
    else if (flag1 && flag2)
        K_local = sumArray(K_temp1, K_temp2);
    else if (flag1)
        K_local = K_temp1;
    else if (flag2)
        K_local = K_temp2;

    return K_local;
}

template <int nlayer>
std::array<std::array<double, NDIM>, NDIM> Stiffness<nlayer>::localStiffnessFD(Particle<nlayer> *pi, Particle<nlayer> *pj)
{
    std::vector<double> K_ij;
    std::array<double, NDIM> xyz_temp = pj->xyz;

    // common conns particles
    // std::vector<Particle<nlayer> *> common_conns;
    // set_intersection(pi->conns.begin(), pi->conns.end(), pj->conns.begin(), pj->conns.end(), back_inserter(common_conns));

    pi->updateParticleForce(); // update pi particle internal forces

    std::array<double, NDIM> Pin_temp = pi->Pin;
    // printf("id: %d, %f, %f, %f\n", pi->id, Pin_temp[0], Pin_temp[1], Pin_temp[2]);

    for (int r = 0; r < pj->cell.dim; ++r)
    {
        xyz_temp[r] += EPS * pj->cell.radius; // forward-difference
        pj->moveTo(xyz_temp);                 // move to a new position

        // update bforce of all common conns (need to update all nonlocal geometry and state variables before bforce calculation)
        for (Particle<nlayer> *pjj : pj->conns)
            pjj->updateBondsGeometry(); // update all bond information, e.g., dL, dL_total

        for (Particle<nlayer> *pjj : pj->conns)
            pjj->updateParticleStateVariables();

        for (Particle<nlayer> *pjj : pj->conns)
            pjj->updateBondsForce(); // update all bond forces

        pi->updateParticleForce(); // update pi particle internal forces

        for (int s = 0; s < pi->cell.dim; s++)
        {
            double K_value = (pi->Pin[s] - Pin_temp[s]) / EPS / (pi->cell.radius);
            K_ij.push_back(K_value); // for each K_ij, index order is 11, 12, 13, 21, ..., 32, 33
            // if(pi->id == 20)
            // printf("%f, ", K_value);
        }

        // resume the particle's original state
        xyz_temp[r] -= EPS * pj->cell.radius; // move back the particle position
        pj->moveTo(xyz_temp);
        pj->resumeParticleState();
    } // K_ij has finished

    std::array<std::array<double, NDIM>, NDIM> K_local{0};
    for (int r = 0; r < pi->cell.dim; r++)
        for (int s = 0; s < pi->cell.dim; s++)
            K_local[r][s] = K_ij[pi->cell.dim * s + r];
    return K_local;
}

// template <int nlayer>
// void Stiffness<nlayer>::calcStiffness2D(std::vector<Particle<nlayer> *> &pt_sys)
// {
//     // std::fill(K_global.begin(), K_global.end(), 0.0);

// #pragma omp parallel for if (mode == StiffnessMode::Analytical)
//     for (const auto &pi_iterator : pt_sys | indexed(0))
//     {
//         Particle<nlayer> *pi = pi_iterator.value();

//         for (const auto &pj_iterator : pi->conns | indexed(0))
//         {
//             int idx_j = (int)pj_iterator.index();
//             Particle<nlayer> *pj = pj_iterator.value();
//             std::array<std::array<double, NDIM>, NDIM> K_local = localStiffness(pi, pj);

//             // if (pi->id == 0)
//             // {
//             //     printf("j: %d, %f, %f, %f\n%f, %f, %f\n%f, %f, %f\n\n", pj->id, K_local[0][0], K_local[0][1], K_local[0][2],
//             //            K_local[1][0], K_local[1][1], K_local[1][2],
//             //            K_local[2][0], K_local[2][1], K_local[2][2]);
//             // }

//             if (pi->id == pj->id)
//             {
//                 K_global[K_pointer[pi->id]] += K_local[0][0];
//                 K_global[K_pointer[pi->id] + 1] += K_local[0][1];
//                 K_global[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq)] += K_local[1][1];

//                 JK[K_pointer[pi->id]] = pi->cell.dim * (pj->id + 1) - 1;
//                 JK[K_pointer[pi->id] + 1] = pi->cell.dim * (pj->id + 1);
//                 JK[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq)] = pi->cell.dim * (pj->id + 1);
//             }
//             else if (pj->id > pi->id)
//             {
//                 int num1 = pi->nconn_largeq - (pi->nconn - idx_j); // index difference between i and j, in i's conn list
//                 // if (pi->id == 30)
//                 //     printf("%d, ", num1);
//                 K_global[K_pointer[pi->id] + pi->cell.dim * num1] += 0.5 * K_local[0][0];
//                 K_global[K_pointer[pi->id] + pi->cell.dim * num1 + 1] += 0.5 * K_local[0][1];
//                 K_global[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 1] += 0.5 * K_local[1][0];
//                 K_global[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1] += 0.5 * K_local[1][1];

//                 JK[K_pointer[pi->id] + pi->cell.dim * num1] = pi->cell.dim * (pj->id + 1) - 1;
//                 JK[K_pointer[pi->id] + pi->cell.dim * num1 + 1] = pi->cell.dim * (pj->id + 1);
//                 JK[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 1] = pi->cell.dim * (pj->id + 1) - 1;
//                 JK[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1] = pi->cell.dim * (pj->id + 1);
//             }
//             else
//             {
//                 auto pt_i = std::find(pj->conns.begin(), pj->conns.end(), pi);
//                 int num2 = pj->nconn_largeq - (int)std::distance(pt_i, pj->conns.end());
//                 // if (pi->id == 40)
//                 //      printf("%d, ", num2);
//                 K_global[K_pointer[pj->id] + pi->cell.dim * num2] += 0.5 * K_local[0][0];
//                 K_global[K_pointer[pj->id] + pi->cell.dim * num2 + 1] += 0.5 * K_local[1][0];
//                 K_global[K_pointer[pj->id] + pi->cell.dim * (pj->nconn_largeq) + pi->cell.dim * num2 - 1] += 0.5 * K_local[0][1];
//                 K_global[K_pointer[pj->id] + pi->cell.dim * (pj->nconn_largeq) + pi->cell.dim * num2] += 0.5 * K_local[1][1];
//             }

//             // if (pi->id == 40)
//             //     printf("%f, ", K_global[K_pointer[pj->id] + pi->cell.dim * 1 + 1]);
//         }

//         IK[pi->cell.dim * pi->id] = K_pointer[pi->id] + 1;
//         IK[pi->cell.dim * pi->id + 1] = K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + 1;
//         // if (pi->id == 40)
//         //     printf("%d, ", IK[pi->cell.dim * pi->id]);
//     }
//     IK[pt_sys[0]->cell.dim * pt_sys.size()] = K_pointer[pt_sys.size()] + 1;
// }

template <int nlayer>
void Stiffness<nlayer>::calcStiffness(std::vector<Particle<nlayer> *> &pt_sys)
{
    // Estimate the sparsity pattern
    // For simplicity, this example assumes you know the total number of non-zero entries or have a good estimate.
    // A more accurate approach would involve analyzing the connectivity of your system to determine the exact sparsity pattern.
    const int estimatedNonZerosPerRow = pt_sys.size() * pt_sys[0]->cell.dim * pt_sys[0]->cell.nneighbors_AFEM; /* your estimate here */
    ;
    K_matrix.resize(pt_sys.size() * pt_sys[0]->cell.dim, pt_sys.size() * pt_sys[0]->cell.dim);
    K_matrix.reserve(Eigen::VectorXi::Constant(pt_sys.size() * pt_sys[0]->cell.dim, estimatedNonZerosPerRow));

    for (const auto &pi_iterator : pt_sys | indexed(0))
    {
        Particle<nlayer> *pi = pi_iterator.value();

        for (const auto &pj_iterator : pi->conns | indexed(0))
        {
            int idx_j = (int)pj_iterator.index();
            Particle<nlayer> *pj = pj_iterator.value();

            std::array<std::array<double, NDIM>, NDIM> K_local = localStiffness(pi, pj);
            int insertRow = pi->id * pi->cell.dim;
            int insertCol = pj->id * pj->cell.dim;

            for (int i = 0; i < pi->cell.dim; ++i)
            {
                for (int j = 0; j < pj->cell.dim; ++j)
                {
                    K_matrix.insert(insertRow + i, insertCol + j) = K_local[i][j];
                }
            }
        }
    }

    K_matrix.makeCompressed(); // Optional, depending on subsequent operations
}

template <int nlayer>
void Stiffness<nlayer>::updateStiffnessDispBC(std::vector<Particle<nlayer> *> &pt_sys)
{
    // extract the diagonal vector of stiffness matrix
    Eigen::VectorXd diag = K_matrix.diagonal();
    double norm_diag = diag.norm();

    // update the stiffness matrix
    for (Particle<nlayer> *pi : pt_sys)
    {
        for (int k = 0; k < pi->cell.dim; k++)
        {
            if (pi->disp_constraint[k] == 1 || abs(pi->damage_visual - 1.0) < EPS) // if the particle is under disp bc or totally damaged
            {
                // zero out the entire row and column for the current particle
                int row = pi->id * pi->cell.dim + k;

                for (int col = 0; col < K_matrix.outerSize(); ++col)
                {
                    K_matrix.coeffRef(row, col) = 0.0; // Zero out the row
                    K_matrix.coeffRef(col, row) = 0.0; // Zero out the column
                }

                K_matrix.coeffRef(row, row) = norm_diag;
            }
        }
    }

    K_matrix.makeCompressed();
}

#endif
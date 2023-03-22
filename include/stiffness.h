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
    MKL_INT *IK;
    int *K_pointer; // start index for each particle in the global stiffness matrix
    double *residual;
    std::vector<MKL_INT> JK;
    std::vector<double> K_global;

    void initialize(std::vector<Particle<nlayer> *> &pt_sys);
    void reset(std::vector<Particle<nlayer> *> &pt_sys);
    void calcStiffness3D(std::vector<Particle<nlayer> *> &pt_sys);
    void calcStiffness2D(std::vector<Particle<nlayer> *> &pt_sys);
    void updateStiffnessDispBC(std::vector<Particle<nlayer> *> &pt_sys);

    std::array<std::array<double, NDIM>, NDIM> localStiffness(Particle<nlayer> *pi, Particle<nlayer> *pj);
    std::array<std::array<double, NDIM>, NDIM> localStiffnessFD(Particle<nlayer> *pi, Particle<nlayer> *pj);
    std::array<std::array<double, NDIM>, NDIM> localStiffnessANA(Particle<nlayer> *pi, Particle<nlayer> *pj);

    Stiffness(std::vector<Particle<nlayer> *> &pt_sys, StiffnessMode p_mode)
    { // Given a particle system, construct a stiffness matrix and solver
        mode = p_mode;
        initialize(pt_sys);
    }

    ~Stiffness()
    {
        delete[] IK;
        delete[] residual;
        delete[] K_pointer;
    }
};

template <int nlayer>
void Stiffness<nlayer>::initialize(std::vector<Particle<nlayer> *> &pt_sys)
{
    K_pointer = new int[pt_sys.size() + 1];
    K_pointer[pt_sys[0]->id] = 0;
    for (auto pt : pt_sys)
    {
        if (pt->cell.dim == 2)
            K_pointer[pt->id + 1] = K_pointer[pt->id] + (pt->cell.dim) * (pt->cell.dim) * (pt->nconn_largeq) - 1;
        else
            K_pointer[pt->id + 1] = K_pointer[pt->id] + (pt->cell.dim) * (pt->cell.dim) * (pt->nconn_largeq) - 3;
    }

    IK = new MKL_INT[pt_sys[0]->cell.dim * pt_sys.size() + 1];
    JK = std::vector<MKL_INT>(K_pointer[pt_sys.size()]);
    residual = new double[pt_sys[0]->cell.dim * pt_sys.size()];
    K_global = std::vector<double>(K_pointer[pt_sys.size()]);
}

template <int nlayer>
void Stiffness<nlayer>::reset(std::vector<Particle<nlayer> *> &pt_sys)
{
    for (auto pt : pt_sys)
    {
        if (pt->cell.dim == 2)
            K_pointer[pt->id + 1] = K_pointer[pt->id] + (pt->cell.dim) * (pt->cell.dim) * (pt->nconn_largeq) - 1;
        else
            K_pointer[pt->id + 1] = K_pointer[pt->id] + (pt->cell.dim) * (pt->cell.dim) * (pt->nconn_largeq) - 3;
    }

    JK = std::vector<MKL_INT>(K_pointer[pt_sys.size()]);
    K_global = std::vector<double>(K_pointer[pt_sys.size()]);
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

        // update bforce of all common conns
        for (Particle<nlayer> *pjj : pj->conns) // pj->conns or common_conns
        {
            pjj->updateBondsGeometry(); // update all bond information, e.g., dL, ddL
            pjj->updateBondsForce();    // update all bond forces
        }

        pi->updateParticleForce(); // update pi particle internal forces

        for (int s = 0; s < pi->cell.dim; s++)
        {
            double K_value = (pi->Pin[s] - Pin_temp[s]) / EPS / (pi->cell.radius);
            K_ij.push_back(K_value); // for each K_ij, index order is 11, 12, 13, 21, ..., 32, 33
            // if(pi->id == 20)
            // printf("%f, ", K_value);
        }

        xyz_temp[r] -= EPS * pj->cell.radius; // move back the particle position
        pj->moveTo(xyz_temp);
        for (Particle<nlayer> *pjj : pj->conns) // pj->conns or common_conns
            pjj->resumeParticle();
    } // K_ij has finished

    std::array<std::array<double, NDIM>, NDIM> K_local{0};
    for (int r = 0; r < pi->cell.dim; r++)
        for (int s = 0; s < pi->cell.dim; s++)
            K_local[r][s] = K_ij[pi->cell.dim * s + r];
    return K_local;
}

template <int nlayer>
void Stiffness<nlayer>::calcStiffness3D(std::vector<Particle<nlayer> *> &pt_sys)
{
    // std::fill(K_global.begin(), K_global.end(), 0.0);

#pragma omp parallel for if (mode == StiffnessMode::Analytical)
    for (const auto &pi_iterator : pt_sys | indexed(0))
    {
        Particle<nlayer> *pi = pi_iterator.value();

        for (const auto &pj_iterator : pi->conns | indexed(0))
        {
            int idx_j = (int)pj_iterator.index();
            Particle<nlayer> *pj = pj_iterator.value();
            std::array<std::array<double, NDIM>, NDIM> K_local = localStiffness(pi, pj);

            // if (pi->id == 30)
            // {
            //     printf("j: %d, %f, %f, %f\n%f, %f, %f\n%f, %f, %f\n\n", pj->id, K_local[0][0], K_local[0][1], K_local[0][2],
            //            K_local[1][0], K_local[1][1], K_local[1][2],
            //            K_local[2][0], K_local[2][1], K_local[2][2]);
            // }

            if (pi->id == pj->id)
            {
                K_global[K_pointer[pi->id]] += K_local[0][0];
                K_global[K_pointer[pi->id] + 1] += K_local[0][1];
                K_global[K_pointer[pi->id] + 2] += K_local[0][2];
                K_global[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq)] += K_local[1][1];
                K_global[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + 1] += K_local[1][2];
                K_global[K_pointer[pi->id] + 2 * pi->cell.dim * (pi->nconn_largeq) - 1] += K_local[2][2];

                JK[K_pointer[pi->id]] = pi->cell.dim * (pj->id + 1) - 2;
                JK[K_pointer[pi->id] + 1] = pi->cell.dim * (pj->id + 1) - 1;
                JK[K_pointer[pi->id] + 2] = pi->cell.dim * (pj->id + 1);
                JK[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq)] = pi->cell.dim * (pj->id + 1) - 1;
                JK[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + 1] = pi->cell.dim * (pj->id + 1);
                JK[K_pointer[pi->id] + 2 * pi->cell.dim * (pi->nconn_largeq) - 1] = pi->cell.dim * (pj->id + 1);
            }
            else if (pj->id > pi->id)
            {
                int num1 = pi->nconn_largeq - (pi->nconn - idx_j); // index difference between i and j, in i's conn list
                // if (pi->id == 30)
                //     printf("%d, ", num1);
                K_global[K_pointer[pi->id] + pi->cell.dim * num1] += 0.5 * K_local[0][0];
                K_global[K_pointer[pi->id] + pi->cell.dim * num1 + 1] += 0.5 * K_local[0][1];
                K_global[K_pointer[pi->id] + pi->cell.dim * num1 + 2] += 0.5 * K_local[0][2];
                K_global[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 1] += 0.5 * K_local[1][0];
                K_global[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1] += 0.5 * K_local[1][1];
                K_global[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 + 1] += 0.5 * K_local[1][2];
                K_global[K_pointer[pi->id] + 2 * pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 3] += 0.5 * K_local[2][0];
                K_global[K_pointer[pi->id] + 2 * pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 2] += 0.5 * K_local[2][1];
                K_global[K_pointer[pi->id] + 2 * pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 1] += 0.5 * K_local[2][2];

                JK[K_pointer[pi->id] + pi->cell.dim * num1] = pi->cell.dim * (pj->id + 1) - 2;
                JK[K_pointer[pi->id] + pi->cell.dim * num1 + 1] = pi->cell.dim * (pj->id + 1) - 1;
                JK[K_pointer[pi->id] + pi->cell.dim * num1 + 2] = pi->cell.dim * (pj->id + 1);
                JK[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 1] = pi->cell.dim * (pj->id + 1) - 2;
                JK[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1] = pi->cell.dim * (pj->id + 1) - 1;
                JK[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 + 1] = pi->cell.dim * (pj->id + 1);
                JK[K_pointer[pi->id] + 2 * pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 3] = pi->cell.dim * (pj->id + 1) - 2;
                JK[K_pointer[pi->id] + 2 * pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 2] = pi->cell.dim * (pj->id + 1) - 1;
                JK[K_pointer[pi->id] + 2 * pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 1] = pi->cell.dim * (pj->id + 1);
            }
            else
            {
                auto pt_i = std::find(pj->conns.begin(), pj->conns.end(), pi);
                int num2 = pj->nconn_largeq - (int)std::distance(pt_i, pj->conns.end());
                // if (pi->id == 40)
                //      printf("%d, ", num2);
                K_global[K_pointer[pj->id] + pi->cell.dim * num2] += 0.5 * K_local[0][0];
                K_global[K_pointer[pj->id] + pi->cell.dim * num2 + 1] += 0.5 * K_local[1][0];
                K_global[K_pointer[pj->id] + pi->cell.dim * num2 + 2] += 0.5 * K_local[2][0];
                K_global[K_pointer[pj->id] + pi->cell.dim * (pj->nconn_largeq) + pi->cell.dim * num2 - 1] += 0.5 * K_local[0][1];
                K_global[K_pointer[pj->id] + pi->cell.dim * (pj->nconn_largeq) + pi->cell.dim * num2] += 0.5 * K_local[1][1];
                K_global[K_pointer[pj->id] + pi->cell.dim * (pj->nconn_largeq) + pi->cell.dim * num2 + 1] += 0.5 * K_local[2][1];
                K_global[K_pointer[pj->id] + 2 * pi->cell.dim * (pj->nconn_largeq) + pi->cell.dim * num2 - 3] += 0.5 * K_local[0][2];
                K_global[K_pointer[pj->id] + 2 * pi->cell.dim * (pj->nconn_largeq) + pi->cell.dim * num2 - 2] += 0.5 * K_local[1][2];
                K_global[K_pointer[pj->id] + 2 * pi->cell.dim * (pj->nconn_largeq) + pi->cell.dim * num2 - 1] += 0.5 * K_local[2][2];
            }

            // if (pi->id == 40)
            //     printf("%f, ", K_global[K_pointer[pj->id] + pi->cell.dim * 1 + 1]);
        }

        IK[pi->cell.dim * pi->id] = K_pointer[pi->id] + 1;
        IK[pi->cell.dim * pi->id + 1] = K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + 1;
        IK[pi->cell.dim * pi->id + 2] = K_pointer[pi->id] + 2 * pi->cell.dim * (pi->nconn_largeq);
        // if (pi->id == 40)
        //     printf("%d, ", IK[pi->cell.dim * pi->id]);
    }
    IK[pt_sys[0]->cell.dim * pt_sys.size()] = K_pointer[pt_sys.size()] + 1;
}

template <int nlayer>
void Stiffness<nlayer>::calcStiffness2D(std::vector<Particle<nlayer> *> &pt_sys)
{
    // std::fill(K_global.begin(), K_global.end(), 0.0);

#pragma omp parallel for if (mode == StiffnessMode::Analytical)
    for (const auto &pi_iterator : pt_sys | indexed(0))
    {
        Particle<nlayer> *pi = pi_iterator.value();

        for (const auto &pj_iterator : pi->conns | indexed(0))
        {
            int idx_j = (int)pj_iterator.index();
            Particle<nlayer> *pj = pj_iterator.value();
            std::array<std::array<double, NDIM>, NDIM> K_local = localStiffness(pi, pj);

            // if (pi->id == 0)
            // {
            //     printf("j: %d, %f, %f, %f\n%f, %f, %f\n%f, %f, %f\n\n", pj->id, K_local[0][0], K_local[0][1], K_local[0][2],
            //            K_local[1][0], K_local[1][1], K_local[1][2],
            //            K_local[2][0], K_local[2][1], K_local[2][2]);
            // }

            if (pi->id == pj->id)
            {
                K_global[K_pointer[pi->id]] += K_local[0][0];
                K_global[K_pointer[pi->id] + 1] += K_local[0][1];
                K_global[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq)] += K_local[1][1];

                JK[K_pointer[pi->id]] = pi->cell.dim * (pj->id + 1) - 1;
                JK[K_pointer[pi->id] + 1] = pi->cell.dim * (pj->id + 1);
                JK[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq)] = pi->cell.dim * (pj->id + 1);
            }
            else if (pj->id > pi->id)
            {
                int num1 = pi->nconn_largeq - (pi->nconn - idx_j); // index difference between i and j, in i's conn list
                // if (pi->id == 30)
                //     printf("%d, ", num1);
                K_global[K_pointer[pi->id] + pi->cell.dim * num1] += 0.5 * K_local[0][0];
                K_global[K_pointer[pi->id] + pi->cell.dim * num1 + 1] += 0.5 * K_local[0][1];
                K_global[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 1] += 0.5 * K_local[1][0];
                K_global[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1] += 0.5 * K_local[1][1];

                JK[K_pointer[pi->id] + pi->cell.dim * num1] = pi->cell.dim * (pj->id + 1) - 1;
                JK[K_pointer[pi->id] + pi->cell.dim * num1 + 1] = pi->cell.dim * (pj->id + 1);
                JK[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 1] = pi->cell.dim * (pj->id + 1) - 1;
                JK[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1] = pi->cell.dim * (pj->id + 1);
            }
            else
            {
                auto pt_i = std::find(pj->conns.begin(), pj->conns.end(), pi);
                int num2 = pj->nconn_largeq - (int)std::distance(pt_i, pj->conns.end());
                // if (pi->id == 40)
                //      printf("%d, ", num2);
                K_global[K_pointer[pj->id] + pi->cell.dim * num2] += 0.5 * K_local[0][0];
                K_global[K_pointer[pj->id] + pi->cell.dim * num2 + 1] += 0.5 * K_local[1][0];
                K_global[K_pointer[pj->id] + pi->cell.dim * (pj->nconn_largeq) + pi->cell.dim * num2 - 1] += 0.5 * K_local[0][1];
                K_global[K_pointer[pj->id] + pi->cell.dim * (pj->nconn_largeq) + pi->cell.dim * num2] += 0.5 * K_local[1][1];
            }

            // if (pi->id == 40)
            //     printf("%f, ", K_global[K_pointer[pj->id] + pi->cell.dim * 1 + 1]);
        }

        IK[pi->cell.dim * pi->id] = K_pointer[pi->id] + 1;
        IK[pi->cell.dim * pi->id + 1] = K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + 1;
        // if (pi->id == 40)
        //     printf("%d, ", IK[pi->cell.dim * pi->id]);
    }
    IK[pt_sys[0]->cell.dim * pt_sys.size()] = K_pointer[pt_sys.size()] + 1;
}

template <int nlayer>
void Stiffness<nlayer>::updateStiffnessDispBC(std::vector<Particle<nlayer> *> &pt_sys)
{
    double *diag = new double[pt_sys[0]->cell.dim * pt_sys.size()]; /* diagonal vector of stiffness matrix */

    /* extract the diagonal vector of stiffness matrix */
    for (Particle<nlayer> *pt : pt_sys)
    {
        diag[(pt->id) * (pt->cell.dim)] = K_global[K_pointer[pt->id]];
        diag[(pt->id) * (pt->cell.dim) + 1] = K_global[K_pointer[pt->id] + (pt->cell.dim) * (pt->nconn_largeq)];
        if (pt->cell.dim == 3)
            diag[(pt->id) * (pt->cell.dim) + 2] = K_global[K_pointer[pt->id] + 2 * (pt->cell.dim) * (pt->nconn_largeq) - 1];
    }

    /* compute the norm of the diagonal */
    double norm_diag = cblas_dnrm2(pt_sys[0]->cell.dim * pt_sys.size(), diag, 1); /* Euclidean norm (L2 norm) */
    delete[] diag;

    /* update the stiffness matrix */
    for (Particle<nlayer> *pi : pt_sys)
    {
        for (int k = 0; k < pt_sys[0]->cell.dim; k++)
        {
            if (pi->disp_constraint[k] == 1 || abs(pi->damage_visual - 1) < EPS) // if the particle is under disp bc or totally damaged
            {
                // printf("disp id: %d, dis: %d\n", pi->id, k);
                for (const auto &pj_iterator : pi->conns | indexed(0))
                {
                    int idx_j = (int)pj_iterator.index();
                    Particle<nlayer> *pj = pj_iterator.value();

                    if (pj->id > pi->id)
                        continue;
                    else if (pj->id == pi->id)
                    {
                        if (k == 0)
                        {
                            int start_idx = (pi->cell.dim == 2) ? K_pointer[pi->id] : (K_pointer[pi->id] + 1);
                            for (int kk = start_idx; kk <= K_pointer[pi->id] + (pi->cell.dim) * (pi->nconn_largeq) - 1; kk++)
                                K_global[kk] = 0.0;
                            K_global[K_pointer[pi->id]] = norm_diag;
                        }
                        else if (k == 1)
                        {
                            K_global[K_pointer[pi->id] + 1] = 0.0;
                            int start_idx = (pi->cell.dim == 2) ? (K_pointer[pi->id] + (pi->cell.dim) * (pi->nconn_largeq)) : (K_pointer[pi->id] + (pi->cell.dim) * (pi->nconn_largeq) + 1);
                            int end_idx = (pi->cell.dim == 2) ? (K_pointer[(pi->id) + 1] - 1) : K_pointer[pi->id] + 2 * (pi->cell.dim) * (pi->nconn_largeq) - 2;
                            for (int kk = start_idx; kk <= end_idx; kk++)
                                K_global[kk] = 0.0;
                            K_global[K_pointer[pi->id] + (pi->cell.dim) * (pi->nconn_largeq)] = norm_diag;
                        }
                        else // k==2
                        {
                            K_global[K_pointer[pi->id] + 2] = 0.0;
                            K_global[K_pointer[pi->id] + (pi->cell.dim) * (pi->nconn_largeq) + 1] = 0.0;
                            for (int kk = K_pointer[pi->id] + 2 * (pi->cell.dim) * (pi->nconn_largeq) - 1; kk <= K_pointer[(pi->id) + 1] - 1; kk++)
                                K_global[kk] = 0.0;
                            K_global[K_pointer[pi->id] + 2 * (pi->cell.dim) * (pi->nconn_largeq) - 1] = norm_diag;
                        }
                    }
                    else // pj->id < pi->id
                    {
                        auto pt_i = std::find(pj->conns.begin(), pj->conns.end(), pi);
                        int num2 = pj->nconn_largeq - (int)std::distance(pt_i, pj->conns.end());
                        K_global[K_pointer[pj->id] + pj->cell.dim * num2 + k] = 0;
                        K_global[K_pointer[pj->id] + pj->cell.dim * (pj->nconn_largeq) + pj->cell.dim * num2 - 1 + k] = 0;
                        if (pi->cell.dim == 3)
                            K_global[K_pointer[pj->id] + 2 * pj->cell.dim * (pj->nconn_largeq) + pj->cell.dim * num2 - 3 + k] = 0;
                    }
                }
                residual[(pi->cell.dim) * (pi->id) + k] = 0.0;
            }
        }
    }
}

#endif
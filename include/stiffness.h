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

template <int nlayer>
class Particle;

template <int nlayer>
class Stiffness
{
    const double eps = 1e-6;
    int mode{0}; // 0 to use finite difference; 1 to use analytical approach to compute local stiffness

public:
    MKL_INT *IK, *JK;
    int *K_pointer; // start index for each particle in the global stiffness matrix
    double *K_global;

    void initialize(std::vector<Particle<nlayer> *> &ptsystem);
    void initializeStiffness3D(std::vector<Particle<nlayer> *> &ptsystem);
    void initializeStiffness2D(std::vector<Particle<nlayer> *> &ptsystem);
    void updateStiffnessDispBC(std::vector<Particle<nlayer> *> &ptsystem);

    std::array<std::array<double, NDIM>, NDIM> localStiffness(Particle<nlayer> *pi, Particle<nlayer> *pj);
    std::array<std::array<double, NDIM>, NDIM> localStiffnessFD(Particle<nlayer> *pi, Particle<nlayer> *pj);
    std::array<std::array<double, NDIM>, NDIM> localStiffnessANA(Particle<nlayer> *pi, Particle<nlayer> *pj);

    Stiffness(std::vector<Particle<nlayer> *> &ptsystem, int p_mode)
    { // Given a particle system, construct a stiffness matrix and solver

        mode = p_mode;
        initialize(ptsystem);
        if (ptsystem[0]->cell.dim == 2)
            initializeStiffness2D(ptsystem);
        else
            initializeStiffness3D(ptsystem);
    }

    ~Stiffness()
    {
        delete[] IK;
        delete[] JK;
        delete[] K_global;
        delete[] K_pointer;
    }
};

template <int nlayer>
void Stiffness<nlayer>::initialize(std::vector<Particle<nlayer> *> &ptsystem)
{
    K_pointer = new int[ptsystem.size() + 1];
    K_pointer[0] = 0;
    for (auto pt : ptsystem)
    {
        if (pt->cell.dim == 2)
            K_pointer[pt->id + 1] = K_pointer[pt->id] + pt->cell.dim * pt->cell.dim * pt->nconn_largeq - 1;
        else
            K_pointer[pt->id + 1] = K_pointer[pt->id] + pt->cell.dim * pt->cell.dim * pt->nconn_largeq - 3;
    }

    JK = new MKL_INT[K_pointer[ptsystem.size()]];
    IK = new MKL_INT[ptsystem[0]->cell.dim * ptsystem.size() + 1];
    K_global = new double[K_pointer[ptsystem.size()]];
}

template <int nlayer>
std::array<std::array<double, NDIM>, NDIM> Stiffness<nlayer>::localStiffness(Particle<nlayer> *pi, Particle<nlayer> *pj)
{
    if (mode == 1)
        return localStiffnessANA(pi, pj);
    return localStiffnessFD(pi, pj); // mode == 0
}

template <int nlayer>
std::array<std::array<double, NDIM>, NDIM> Stiffness<nlayer>::localStiffnessANA(Particle<nlayer> *pi, Particle<nlayer> *pj)
{
    std::array<std::array<double, NDIM>, NDIM> K_local;
    return K_local;
}

template <int nlayer>
std::array<std::array<double, NDIM>, NDIM> Stiffness<nlayer>::localStiffnessFD(Particle<nlayer> *pi, Particle<nlayer> *pj)
{
    std::vector<double> K_ij;
    std::array<double, NDIM> xyz_temp = pj->xyz;

    // common conns particles
    std::vector<Particle<nlayer> *> common_conns;
    set_intersection(pi->conns.begin(), pi->conns.end(), pj->conns.begin(), pj->conns.end(), back_inserter(common_conns));

    pi->updateParticleForce(); // update pi particle internal forces

    std::array<double, NDIM> Pin_temp = pi->Pin;
    // printf("id: %d, %f, %f, %f\n", pi->id, Pin_temp[0], Pin_temp[1], Pin_temp[2]);

    for (int r = 0; r < pj->cell.dim; ++r)
    {
        xyz_temp[r] += eps * pj->cell.radius; // forward-difference
        pj->moveTo(xyz_temp);                 // move to a new position

        // update bforce of all common conns
        for (Particle<nlayer> *pjj : common_conns) // pj->conns)
        {
            pjj->updateBondsGeometry(); // update all bond information, e.g., dL, ddL
            pjj->updateBondsForce();    // update all bond forces
        }

        pi->updateParticleForce(); // update pi particle internal forces

        for (int s = 0; s < pi->cell.dim; s++)
        {
            double K_value = (pi->Pin[s] - Pin_temp[s]) / eps / (pi->cell.radius);
            K_ij.push_back(K_value); // for each K_ij, index order is 11, 12, 13, 21, ..., 32, 33
            // if(pi->id == 20)
            // printf("%f, ", K_value);
        }

        xyz_temp[r] -= eps * pj->cell.radius; // move back the particle position
        pj->moveTo(xyz_temp);
        for (Particle<nlayer> *pjj : common_conns) // pj->conns)
            pjj->resumeParticle();
    } // K_ij has finished

    std::array<std::array<double, NDIM>, NDIM> K_local;
    for (int r = 0; r < NDIM; r++)
        for (int s = 0; s < NDIM; s++)
            K_local[r][s] = K_ij[NDIM * s + r];
    return K_local;
}

template <int nlayer>
void Stiffness<nlayer>::initializeStiffness3D(std::vector<Particle<nlayer> *> &ptsystem)
{
#pragma omp parallel if (mode)
    for (const auto &pi_iterator : ptsystem | indexed(0))
    {
        Particle<nlayer> *pi = pi_iterator.value();

        for (const auto &pj_iterator : pi->conns | indexed(0))
        {
            int idx_j = (int)pj_iterator.index();
            Particle<nlayer> *pj = pj_iterator.value();
            std::array<std::array<double, NDIM>, NDIM> K_local = localStiffness(pi, pj);

            // if (pi->id == 30)
            // {
            //     printf("%f, %f, %f\n", K_local[0][0], K_local[1][1], K_local[2][2]);
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
    IK[ptsystem[0]->cell.dim * ptsystem.size()] = K_pointer[ptsystem.size()] + 1;
}

template <int nlayer>
void Stiffness<nlayer>::initializeStiffness2D(std::vector<Particle<nlayer> *> &ptsystem)
{
}

template <int nlayer>
void Stiffness<nlayer>::updateStiffnessDispBC(std::vector<Particle<nlayer> *> &ptsystem)
{
    double *diag = new double[ptsystem[0]->cell.dim * ptsystem.size()]; /* diagonal vector of stiffness matrix */

    /* extract the diagonal vector of stiffness matrix */
    for (Particle<nlayer> *pt : ptsystem)
    {
        diag[(pt->id) * (pt->cell.dim)] = K_global[K_pointer[pt->id]];
        diag[(pt->id) * (pt->cell.dim) + 1] = K_global[K_pointer[pt->id] + (pt->cell.dim) * (pt->nconn_largeq)];
        if (pt->cell.dim == 3)
            diag[(pt->id) * (pt->cell.dim) + 2] = K_global[K_pointer[pt->id] + 2 * (pt->cell.dim) * (pt->nconn_largeq) - 1];
    }

    /* compute the norm of the diagonal */
    double norm_diag = cblas_dnrm2(ptsystem[0]->cell.dim * ptsystem.size(), diag, 1); /* Euclidean norm (L2 norm) */

    /* update the stiffness matrix */
    for (Particle<nlayer> *pi : ptsystem)
    {
        for (int k = 0; k < ptsystem[0]->cell.dim; k++)
        {
            if (pi->disp_constraint[k] == 1 || pi->frozen == 1)
            {
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
                            for (int kk = K_pointer[pi->id] + 1; kk <= K_pointer[pi->id] + (pi->cell.dim) * (pi->nconn_largeq) - 1; kk++)
                                K_global[kk] = 0.0;
                            K_global[K_pointer[pi->id]] = norm_diag;
                        }
                        else if (k == 1)
                        {
                            K_global[K_pointer[pi->id] + 1] = 0.0;
                            for (int kk = K_pointer[pi->id] + (pi->cell.dim) * (pi->nconn_largeq) + 1; kk <= K_pointer[pi->id] + 2 * (pi->cell.dim) * (pi->nconn_largeq) - 2; kk++)
                                K_global[kk] = 0.0;
                            K_global[K_pointer[pi->id] + (pi->cell.dim) * (pi->nconn_largeq)] = norm_diag;
                        }
                        else // k==2
                        {
                            K_global[K_pointer[pi->id] + 2] = 0.0;
                            K_global[K_pointer[pi->id] + (pi->cell.dim) * (pi->nconn_largeq) + 1] = 0.0;
                            for (int kk = K_pointer[pi->id] + 2 * (pi->cell.dim) * (pi->nconn_largeq) - 1; kk <= K_pointer[pi->id + 1] - 1; kk++)
                                K_global[kk] = 0.0;
                            K_global[K_pointer[pi->id] + 2 * (pi->cell.dim) * (pi->nconn_largeq) - 1] = norm_diag;
                        }
                    }
                    else
                    {
                        int num1 = pi->nconn_largeq - (pi->nconn - idx_j); // index difference between i and j, in i's conn list
                        K_global[K_pointer[pi->id] + pi->cell.dim * num1 + k] = 0;
                        K_global[K_pointer[pi->id] + pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 1 + k] = 0;
                        K_global[K_pointer[pi->id] + 2 * pi->cell.dim * (pi->nconn_largeq) + pi->cell.dim * num1 - 3 + k] = 0;
                    }
                }
            }
            residual[(pi->cell.dim) * (pi->id) + k] = 0.0;
        }
    }
}
#endif
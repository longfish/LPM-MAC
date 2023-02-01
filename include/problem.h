#pragma once
#ifndef PROBLEM_H
#define PROBLEM_H

#include <vector>
#include <array>
#include <algorithm>

#include "lpm.h"
#include "unit_cell.h"
#include "bond.h"

template <int nlayer>
class Bond;

template <int nlayer>
class Particle;

template <int nlayer>
class LPMProblem
{
    // Given a particle system, construct a stiffness matrix, RHS and solver
public:
    int nparticle; // number of particles
    MKL_INT *IK, *JK;
    double *K_global, *residual, *reaction_force;
    std::vector<Particle<nlayer> *> ptsystem; // system of particles

    LPMProblem(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell);

    void createParticles(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell);
    void createBonds();
    void createConnections();
};

template <int nlayer>
LPMProblem<nlayer>::LPMProblem(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell)
{
    createParticles(p_xyz, p_cell);
    createBonds();
    createConnections();
    nparticle = ptsystem.size();
}

template <int nlayer>
void LPMProblem<nlayer>::createParticles(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell)
{
    for (auto xyz : p_xyz)
    {
        Particle<nlayer> *pt = new Particle<nlayer>(xyz[0], xyz[1], xyz[2], p_cell);
        ptsystem.push_back(pt);
    }
}

template <int nlayer>
void LPMProblem<nlayer>::createBonds()
{
#pragma omp parallel for
    for (Particle<nlayer> *p1 : ptsystem)
    {
        for (Particle<nlayer> *p2 : ptsystem)
        {
            double distance = p1->distanceTo(*p2);
            if ((distance < 1.01 * p1->cell.neighbor2_cutoff) && (p1->id != p2->id))
            {
                int layer = 1;
                if (distance < 1.01 * p1->cell.neighbor1_cutoff)
                    layer = 0;

                Bond<nlayer> *bd = new Bond<nlayer>(p1, p2, layer, distance);
                p1->bond_layers[0].push_back(bd);
                p1->neighbors.push_back(p2);
            }
        }
        p1->nb = p1->neighbors.size();
    }
}

template <int nlayer>
void LPMProblem<nlayer>::createConnections()
{
#pragma omp parallel for 
    for (Particle<nlayer> *p1 : ptsystem)
    {
        for (int i = 0; i < nlayer; i++)
        {
            // loop forward bond particles
            for (Bond<nlayer> *bd_fw : p1->bond_layers[i])
            {
                p1->conns.push_back(bd_fw->p2);

                // loop backward bond particles
                for (Bond<nlayer> *bd_bw : bd_fw->p2->bond_layers[i])
                    p1->conns.push_back(bd_bw->p2);
            }
        }
        std::sort(p1->conns.begin(), p1->conns.end(), [](Particle<nlayer> *a, Particle<nlayer> *b)
                  { return a->id < b->id; });
        p1->conns.erase(unique(p1->conns.begin(), p1->conns.end()), p1->conns.end());
        p1->nconn = p1->conns.size();
    }
}

#endif
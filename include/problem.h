#pragma once
#ifndef PROBLEM_H
#define PROBLEM_H

#include <vector>
#include <array>

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
    std::vector<Particle<nlayer>> ptsystem; // system of particles

    LPMProblem(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell)
    {
        createParticles(p_xyz, p_cell);
        createBonds();
        createConnections();
        nparticle = ptsystem.size();
    }

    void createParticles(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell);
    void createBonds();
    void createConnections();
};

#endif
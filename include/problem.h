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
    int nparticle{0};                                    // number of particles
    MKL_INT *IK, *JK;
    double *K_global, *residual, *reaction_force;
    std::vector<Particle<nlayer>> psystem;                 // system of particles

    LPMProblem(){}
    LPMProblem(std::vector<std::array<double, NDIM>> p_xyz, UnitCell p_cell)
    {
        assemblyParticles(p_xyz, p_cell);
        nparticle = psystem.size();
    }

    void assemblyParticles(std::vector<std::array<double, NDIM>> p_xyz, UnitCell p_cell);
};

#endif
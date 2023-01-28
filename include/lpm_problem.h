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
    int m_nparticle{0};                                    // number of particles
    MKL_INT *IK, *JK;
    double *K_global, *residual, *reaction_force;
    std::vector<Particle<nlayer>> m_psystem;                 // system of particles

    Problem(const std::vector<Particle<nlayer>> &psystem)
    {
        m_psystem = psystem;
        m_nparticle = psystem.size()
    }
};

#endif
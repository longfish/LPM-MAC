#pragma once
#ifndef SOLVER_H
#define SOLVER_H

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
class Solver
{
public:
    MKL_INT *IK, *JK;
    double *K_global, *residual, *reaction_force;

    void updateStiffness3DFiniteDifference();
    void updateStiffness2DFiniteDifference();

    Solver(std::vector<Particle<nlayer> *> &ptsystem)
    { // Given a particle system, construct a stiffness matrix and solver
        if (ptsystem[0]->cell.dim == 2)
            updateStiffness2DFiniteDifference();
        else
            updateStiffness3DFiniteDifference();
    }

    ~Solver()
    {
        delete[] IK;
        delete[] JK;
        delete[] K_global;
        delete[] residual;
        delete[] reaction_force;
    }
};

#endif
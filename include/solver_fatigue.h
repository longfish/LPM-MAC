#pragma once
#ifndef SOLVER_FATIGUE_H
#define SOLVER_FATIGUE_H

#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include "lpm.h"
#include "load_step.h"
#include "stiffness.h"
#include "unit_cell.h"
#include "assembly.h"
#include "solver.h"

template <int nlayer>
class SolverFatigue : public Solver<nlayer>
{
public:
    Stiffness<nlayer> stiffness;

    SolverFatigue(Assembly<nlayer> &ass, const StiffnessMode &p_stiff_mode, const SolverMode &p_sol_mode, const std::string &p_dumpFile)
        : Solver<nlayer>{ass, p_stiff_mode, p_sol_mode, p_dumpFile} {}

};


#endif
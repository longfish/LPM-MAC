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
    SolverFatigue(Assembly<nlayer> &p_ass, const StiffnessMode &p_stiff_mode, const SolverMode &p_sol_mode, const std::string &p_dumpFile, const int &p_niter, const double &p_tol)
        : Solver<nlayer>{p_ass, p_stiff_mode, p_sol_mode, p_dumpFile, p_niter, p_tol} {}

    int getCycleJumpN();
    int updateFatigueDamage();
    bool solveProblemStep(LoadStep<nlayer> &load);
    void solveProblemOneCycle(std::vector<LoadStep<nlayer>> &load);
    void solveProblemCyclic(std::vector<std::vector<LoadStep<nlayer>>> &cycle_loads);
    void solveProblemStatic(std::vector<LoadStep<nlayer>> &load, int &start_index);
};

template <int nlayer>
int SolverFatigue<nlayer>::getCycleJumpN()
{
    std::vector<int> vec_njump;
    for (Particle<nlayer> *pt : this->ass.pt_sys)
        vec_njump.push_back(pt->calcNCycleJump());

    return *std::min_element(vec_njump.begin(), vec_njump.end()); // return minimum cycle jump number
}

template <int nlayer>
int SolverFatigue<nlayer>::updateFatigueDamage()
{
    this->ass.updateStateVar(); // compute the damage rate
    int dN = getCycleJumpN();   // determine the cycle-jump number
    std::cout << dN << std::endl;

    // set the ncycle_jump of the particle
    for (Particle<nlayer> *pt : this->ass.pt_sys)
    {
        pt->ncycle_jump = dN;               // set the ncycle_jump
        pt->updateParticleStateVariables(); // update the current damage
        pt->ncycle_jump = 0;                // set the ncycle_jump to be zero
    }

    return dN;
}

template <int nlayer>
void SolverFatigue<nlayer>::solveProblemStatic(std::vector<LoadStep<nlayer>> &load, int &start_index)
{
    int n_step = load.size();
    bool is_converged{true}; // flag to determine whether need to cut the loading into half
    for (int i = 0; i < n_step; i++)
    {
    restart:
        printf("Loading step-%d, iteration starts:\n", i + start_index);
        double t1 = omp_get_wtime();

        is_converged = solveProblemStep(load[i]);
        if (!is_converged)
        {
            this->ass.resetStateVar(true); // reset the state_var to last converged ones
            load[i].loadCutHalf();
            load.insert(load.begin() + i, load[i]);
            ++n_step;
            printf("Step-%d not converging\n\n", i + start_index);
            goto restart;
        }

        std::cout << 8613 << ',' << this->ass.pt_sys[8613]->Pin[0] << ',' << this->ass.pt_sys[8613]->Pin[1] << ',' << this->ass.pt_sys[8613]->Pin[2] << std::endl;

        double t2 = omp_get_wtime();
        printf("Loading step %d has finished, spent %f seconds\n\nData output ...\n\n", i + start_index, t2 - t1);

        this->ass.writeDump(this->dumpFile, i);
    }

    start_index += n_step;
}

template <int nlayer>
void SolverFatigue<nlayer>::solveProblemCyclic(std::vector<std::vector<LoadStep<nlayer>>> &cycle_loads)
{
    int N{0};
    int n_jump{1};
    do
    {
        solveProblemOneCycle(cycle_loads[N]);
        // n_jump = updateFatigueDamage();
        // if (n_jump < 1)
        //     n_jump = 1;
        // else if (n_jump > 2 * 1e5)
        //     n_jump = 2 * 1e5;
        N += n_jump;

        this->ass.updateForceState();
        this->ass.writeDump(this->dumpFile, N);

    } while (N < cycle_loads.size());
}

template <int nlayer>
void SolverFatigue<nlayer>::solveProblemOneCycle(std::vector<LoadStep<nlayer>> &load)
{
    int n_step = load.size();
    this->ass.clearStateVar(false); // clear state variables, and keep the particle damage value unchanged

    printf("Cycle starts:\n");
    for (int i = 0; i < n_step; ++i)
    {
        printf("Substep-%d, iteration starts:\n", i + 1);
        double t1 = omp_get_wtime();
        bool is_converged = solveProblemStep(load[i], i + 1);

        this->ass.updateStateVar(); // update the state variables in current cycle

        double t2 = omp_get_wtime();
        printf("Loading step %d has finished, spent %f seconds\n\n", i + 1, t2 - t1);
    }
}

template <int nlayer>
bool SolverFatigue<nlayer>::solveProblemStep(LoadStep<nlayer> &load_step)
{
    // keep the damage unchanged, update the deformation field

    bool new_broken{false};
    int n_newton{0};

    this->updateForceBC(load_step);
    this->updateDisplacementBC(load_step);

    do
    {
        // update the stiffness matrix using current state variables (bdamage)
        this->stiffness.reset(this->ass.pt_sys);
        double t11 = omp_get_wtime();
        if (this->ass.pt_sys[0]->cell.dim == 2)
            this->stiffness.calcStiffness2D(this->ass.pt_sys);
        else
            this->stiffness.calcStiffness3D(this->ass.pt_sys);
        this->stiffness.updateStiffnessDispBC(this->ass.pt_sys);
        double t12 = omp_get_wtime();
        printf("Stiffness matrix calculation costs %f seconds\n", t12 - t11);

        // balance the system using current bond configuration and state variables
        n_newton = this->NewtonIteration();
        if (n_newton >= this->max_NR_iter)
            break; // abnormal return

        new_broken = this->ass.updateBrokenBonds();
        if (new_broken)
            printf("Updating broken bonds\n");
        this->ass.updateGeometry();
        this->ass.updateForceState();
    } while (new_broken);

    if (n_newton < this->max_NR_iter)
        return true; // normal return
    else
        return false; // abnormal return
}

#endif
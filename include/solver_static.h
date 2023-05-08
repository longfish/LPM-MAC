#pragma once
#ifndef SOLVER_STATIC_H
#define SOLVER_STATIC_H

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
class SolverStatic : public Solver<nlayer>
{
    const int max_iter = 30;      /* maximum Newton iteration number */
    const double tol_iter = 1e-5; /* newton iteration tolerance */

public:
    SolverStatic(Assembly<nlayer> &p_ass, const StiffnessMode &p_stiff_mode, const SolverMode &p_sol_mode, const std::string &p_dumpFile)
        : Solver<nlayer>{p_ass, p_stiff_mode, p_sol_mode, p_dumpFile} {}

    bool solveProblemStep(LoadStep<nlayer> &load, int &dump_step);
    void solveProblem(std::vector<LoadStep<nlayer>> &load);

    int NewtonIteration(); // return number of Newton iterations
};

template <int nlayer>
void SolverStatic<nlayer>::solveProblem(std::vector<LoadStep<nlayer>> &load)
{
    // ass.writeDump(dumpFile, 0);

    int n_step = load.size(), dump_step{0};
    bool is_converged{true}; // flag to determine whether need to cut the loading into half
    for (int i = 0; i < n_step; i++)
    {
    restart:
        printf("Loading step-%d, iteration starts:\n", i + 1);
        double t1 = omp_get_wtime();

        is_converged = solveProblemStep(load[i], dump_step);
        if (!is_converged)
        {
            this->ass.resetStateVar(true); // reset the state_var to last converged ones
            load[i].loadCutHalf();
            load.insert(load.begin() + i, load[i]);
            ++n_step;
            printf("Step-%d not converging\n\n", i + 1);
            goto restart;
        }

        // ass.updateStateVar();
        this->ass.storeStateVar(); // store converged state variables

        double t2 = omp_get_wtime();
        printf("Loading step %d has finished, spent %f seconds\n\nData output ...\n\n", i + 1, t2 - t1);

        // ass.writeDump(dumpFile, i);
    }
}

// couple damage and bond stretch
template <int nlayer>
bool SolverStatic<nlayer>::solveProblemStep(LoadStep<nlayer> &load_step, int &dump_step)
{
    bool new_damaged{false}, new_broken{false};
    int m{0}, n_newton{0};

    this->updateForceBC(load_step);
    this->updateDisplacementBC(load_step);

    while (m == 0 || new_damaged || new_broken)
    {
        this->ass.writeDump(this->dumpFile, dump_step++);

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
        n_newton = NewtonIteration();
        if (n_newton >= max_iter)
            break;

        if (m % 2 == 0)
        {
            new_damaged = this->ass.updateStateVar();
            if (new_damaged)
                printf("Updating damage\n");
            this->ass.updateGeometry();
            this->ass.updateForceState();
            ++m;
            if (new_damaged)
                continue;
        }
        if (m % 2 == 1)
        {
            new_broken = this->ass.updateBrokenBonds();
            if (new_broken)
                printf("Updating broken bonds\n");
            this->ass.updateGeometry();
            this->ass.updateForceState();
            ++m;
        }
    }

    if (n_newton < max_iter)
        return true; // normal return
    else
        return false; // abnormal return
}

template <int nlayer>
int SolverStatic<nlayer>::NewtonIteration()
{
    // update the force state and residual
    this->ass.updateGeometry();
    this->ass.updateForceState();
    this->updateRR();

    // compute the Euclidean norm (L2 norm)
    double norm_residual = cblas_dnrm2(this->problem_size, this->stiffness.residual, 1);
    double norm_reaction_force = cblas_dnrm2(this->reaction_force.size(), this->reaction_force.data(), 1);
    double tol_multiplier = MAX(norm_residual, norm_reaction_force);
    char tempChar1[] = "residual", tempChar2[] = "reaction";
    printf("|  Norm of residual is %.5e, norm of reaction is %.5e, tolerance criterion is based on ", norm_residual, norm_reaction_force);
    if (norm_residual > norm_reaction_force)
        printf("%s force\n", tempChar1);
    else
        printf("%s force\n", tempChar2);

    int ni{0};
    while (norm_residual > tol_iter * tol_multiplier)
    {
        if (++ni > max_iter)
            return max_iter; // abnormal return

        printf("|  |  Iteration-%d: ", ni);
        this->solveLinearSystem(); // solve for the incremental displacement

        this->ass.updateGeometry();
        this->ass.updateForceState();
        this->updateRR(); /* update the RHS risidual force vector */
        norm_residual = cblas_dnrm2(this->problem_size, this->stiffness.residual, 1);
        printf("|  |  Norm of residual is %.3e, residual ratio is %.3e\n", norm_residual, norm_residual / tol_multiplier);
    }

    return ni; // normal return, return number of iterations
}

#endif
#pragma once
#ifndef SOLVER_STATIC_NODAMAGE_H
#define SOLVER_STATIC_NODAMAGE_H

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
class SolverStaticNodamage : public Solver<nlayer>
{
public:
    SolverStaticNodamage(Assembly<nlayer> &p_ass, const StiffnessMode &p_stiff_mode, const SolverMode &p_sol_mode, const std::string &p_dumpFile, const int &p_niter, const double &p_tol)
        : Solver<nlayer>{p_ass, p_stiff_mode, p_sol_mode, p_dumpFile, p_niter, p_tol} {}

    bool solveProblemStep(LoadStep<nlayer> &load, int &dump_step);
    void solveProblem(std::vector<LoadStep<nlayer>> &load, int &start_index);
};

template <int nlayer>
void SolverStaticNodamage<nlayer>::solveProblem(std::vector<LoadStep<nlayer>> &load, int &start_index)
{

    int n_step = load.size(), dump_step{0};
    bool is_converged{true}; // flag to determine whether need to cut the loading into half
    for (int i = start_index; i < start_index + n_step; i++)
    {
    restart:
        printf("Loading step-%d, iteration starts:\n", i);
        double t1 = omp_get_wtime();

        is_converged = solveProblemStep(load[i], dump_step);
        if (!is_converged)
        {
            load[i].loadCutHalf();
            load.insert(load.begin() + i, load[i]);
            ++n_step;
            printf("Step-%d not converging\n\n", i);
            goto restart;
        }

        std::cout << 8613 << ',' << this->ass.pt_sys[8613]->Pin[0] << ',' << this->ass.pt_sys[8613]->Pin[1] << ',' << this->ass.pt_sys[8613]->Pin[2] << std::endl;

        double t2 = omp_get_wtime();
        printf("Loading step %d has finished, spent %f seconds\n\nData output ...\n\n", i, t2 - t1);

        this->ass.writeDump(this->dumpFile, i);
    }

    start_index += n_step;
}

// couple damage and bond stretch
template <int nlayer>
bool SolverStaticNodamage<nlayer>::solveProblemStep(LoadStep<nlayer> &load_step, int &dump_step)
{
    int n_newton{0};

    this->updateForceBC(load_step);
    this->updateDisplacementBC(load_step);

    // update the stiffness matrix
    this->stiffness.reset(this->ass.pt_sys);
    double t11 = omp_get_wtime();
    if (this->ass.pt_sys[0]->cell.dim == 2)
        this->stiffness.calcStiffness2D(this->ass.pt_sys);
    else
        this->stiffness.calcStiffness3D(this->ass.pt_sys);
    this->stiffness.updateStiffnessDispBC(this->ass.pt_sys);
    double t12 = omp_get_wtime();
    printf("Stiffness matrix calculation costs %f seconds\n", t12 - t11);

    // balance the system using current bond configuration
    n_newton = this->NewtonIteration();
    if (n_newton >= this->max_NR_iter)
        return false; // abnormal return

    this->ass.updateGeometry();
    this->ass.updateForceState();
    return true; // normal return
}

#endif
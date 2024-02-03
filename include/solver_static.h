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
public:
    int undamaged_pt_type{0};

    SolverStatic(const int &p_undamaged, Assembly<nlayer> &p_ass, const StiffnessMode &p_stiff_mode, const SolverMode &p_sol_mode, const std::string &p_dumpFile, const int &p_niter, const double &p_tol)
        : Solver<nlayer>{p_ass, p_stiff_mode, p_sol_mode, p_dumpFile, p_niter, p_tol}, undamaged_pt_type{p_undamaged} {}

    bool updateStaticDamage();
    bool solveProblemStep(LoadStep<nlayer> &load, double &dt);
    void solveProblem(std::vector<LoadStep<nlayer>> &load, int &start_index);
};

template <int nlayer>
bool SolverStatic<nlayer>::updateStaticDamage()
{
    // update nonlocal damage dot
    for (Particle<nlayer> *p1 : this->ass.pt_sys)
    {
        if (p1->type != undamaged_pt_type)
        {
            double A = 0;
            p1->Ddot_nonlocal = 0; // nonlocal damage dot
            for (Particle<nlayer> *p2 : p1->neighbors_nonlocal)
            {
                if (p2->type != undamaged_pt_type)
                {
                    double dis = p1->distanceTo(p2);
                    double V_m = p2->cell.particle_volume * p2->nb / p2->cell.nneighbors;
                    p1->Ddot_nonlocal += p2->Ddot_local * func_phi(dis, p1->nonlocal_L) * V_m;
                    A += func_phi(dis, p1->nonlocal_L) * V_m;
                }
            }
            p1->Ddot_nonlocal /= A;
        }
    }

    // update local-wise damage
    bool any_damaged{false};
    for (Particle<nlayer> *pt : this->ass.pt_sys)
        if (pt->type != undamaged_pt_type)
            any_damaged = pt->updateParticleStaticDamage() || any_damaged;

    return any_damaged;
}

template <int nlayer>
void SolverStatic<nlayer>::solveProblem(std::vector<LoadStep<nlayer>> &load, int &start_index)
{
    // ass.writeDump(dumpFile, 0);

    int n_step = load.size(), dump_step{0};
    bool is_converged{true}; // flag to determine whether need to cut the loading into half
    double dt = 1;
    for (int i = 0; i < n_step; i++)
    {
    restart:
        printf("Loading step-%d, iteration starts:\n", i + start_index);
        double t1 = omp_get_wtime();

        is_converged = solveProblemStep(load[i], dt);
        if (!is_converged)
        {
            this->ass.resetStateVar(true); // reset the state_var to last converged ones
            load[i].loadCutHalf();
            load.insert(load.begin() + i, load[i]);
            ++n_step;
            printf("Step-%d not converging\n\n", i + start_index);
            goto restart;
        }

        this->ass.writeDump(this->dumpFile, i + start_index);

        double t2 = omp_get_wtime();
        printf("Loading step %d has finished, spent %f seconds\n\nData output ...\n\n", i + start_index, t2 - t1);
    }
    start_index += n_step;
}

// couple damage and bond stretch
template <int nlayer>
bool SolverStatic<nlayer>::solveProblemStep(LoadStep<nlayer> &load_step, double &dt)
{
    this->updateForceBC(load_step);
    this->updateDisplacementBC(load_step);

    int n_newton{0};
    bool new_damaged{false};

    do
    {
        // update the stiffness matrix using current state variables (bdamage)
        this->stiffness.reset(this->ass.pt_sys);
        double t11 = omp_get_wtime();
        this->stiffness.calcStiffness(this->ass.pt_sys);
        //this->stiffness.updateStiffnessDispBC(this->ass.pt_sys);
        double t12 = omp_get_wtime();
        printf("Stiffness matrix calculation costs %f seconds\n", t12 - t11);

        // balance the system using current bond configuration and state variables
        n_newton = this->NewtonIteration();
        if (n_newton >= this->max_NR_iter)
            break;

        this->ass.updateStateVar();
        new_damaged = updateStaticDamage();
        this->ass.storeStateVar(); // store converged state variables, last_var = var

        if (new_damaged)
        {
            printf("Updating damage\n");
            this->ass.updateGeometry();
            this->ass.updateForceState();
            // this->ass.resetStateVar(false); // reset var = last_var
        }

        // this->ass.writeDump(this->dumpFile, 0);

    } while (new_damaged);

    // std::cout << this->ass.pt_sys[5087]->damage
    //           << ',' << this->ass.pt_sys[5087]->state_var[0]
    //           << ',' << this->ass.pt_sys[5087]->state_var[1]
    //           << ',' << this->ass.pt_sys[5087]->state_var[2] << std::endl;
    // std::cout << this->ass.pt_sys[8614]->state_var[0] << ',' << this->ass.pt_sys[8614]->state_var[0] << std::endl;
    // std::cout << m << ',' << this->ass.pt_sys[5030]->damage << ',' << this->ass.pt_sys[5030]->state_var[0] << ',' << this->ass.pt_sys[5030]->state_var[1] << std::endl;

    if (n_newton < this->max_NR_iter)
        return true; // normal return
    else
        return false; // abnormal return
}

#endif
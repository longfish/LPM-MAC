#pragma once
#ifndef SOLVER_FATIGUE_H
#define SOLVER_FATIGUE_H

#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include "lpm.h"
#include "load_step.h"
#include "load_step_fatigue.h"
#include "stiffness.h"
#include "unit_cell.h"
#include "assembly.h"
#include "solver.h"

template <int nlayer>
class SolverFatigue : public Solver<nlayer>
{
public:
    int undamaged_pt_type{0};
    TimeMapMode t_mode;                // fatigue time mapping mode
    double tau{0};                     // fatigue time mapping parameter
    std::vector<double> load_spectrum; // fatigue loading data

    SolverFatigue(const int &p_undamaged, Assembly<nlayer> &p_ass, const StiffnessMode &p_stiff_mode, const SolverMode &p_sol_mode, const TimeMapMode &p_t_mode, const double &p_tau, const std::string &p_dumpFile, const int &p_niter, const double &p_tol)
        : undamaged_pt_type{p_undamaged}, Solver<nlayer>{p_ass, p_stiff_mode, p_sol_mode, p_dumpFile, p_niter, p_tol}, t_mode(p_t_mode), tau(p_tau) {}

    std::vector<LoadStep<nlayer>> generateLoadCycle(FatigueLoadType ltype, int N, int n_interval, std::vector<std::vector<Particle<nlayer> *>> pt_groups);

    void readLoad(const std::string &loadFile);
    bool updateFatigueDamage(double dNdt);
    bool solveProblemStep(LoadStep<nlayer> &load);
    void solveProblemOneCycle(std::vector<LoadStep<nlayer>> &load, double dNdt);
    void solveProblemCyclic(FatigueLoadType loadType, int n_interval, std::vector<std::vector<Particle<nlayer> *>> pt_groups);
    void solveProblemStatic(std::vector<LoadStep<nlayer>> &load, int &start_index);
};

template <int nlayer>
void SolverFatigue<nlayer>::readLoad(const std::string &loadFile)
{
    // first line is f_min
    // cyclic loading starts from the 2nd line -> f_max
    // pattern is f_min, f_max, f_min, f_max , ..., f_max, f_min

    FILE *fpt;
    fpt = fopen(loadFile.c_str(), "r+"); /* read-only */

    if (fpt == NULL)
    {
        printf("\'%s\' does not exist!\n", loadFile.c_str());
        exit(1);
    }

    double data;
    while (fscanf(fpt, "%lf", &(data)) > 0)
    {
        load_spectrum.push_back(data);
    }

    fclose(fpt);
}

int func_dNdt_linear(const double &tau)
{
    return (int)(1 / tau);
}

int func_dNdt_exponental(const double &tau, const double &t)
{
    return (int)(1 / tau * exp(t / tau));
}

template <int nlayer>
bool SolverFatigue<nlayer>::updateFatigueDamage(double dNdt)
{
    // update nonlocal damage dot
    for (Particle<nlayer> *p1 : this->ass.pt_sys)
    {
        // p1->Ddot_nonlocal = p1->Ddot_local;
        if (p1->type != undamaged_pt_type)
        {
            double A = 0;
            p1->Ddot_nonlocal = 0; // nonlocal damage dot

            // compute nonlocal damage rate based on Gaussian
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

            // compute nonlocal damage rate using connection averaging
            // for (Particle<nlayer> *p2 : p1->conns)
            // {
            //     double V_m = p2->cell.particle_volume * p2->nb / p2->cell.nneighbors;
            //     p1->Ddot_nonlocal += p2->Ddot_local * V_m;
            //     A += V_m;
            // }

            // compute nonlocal damage rate using neighbor averaging
            // for (Particle<nlayer> *p2 : p1->neighbors)
            // {
            //     double V_m = p2->cell.particle_volume * p2->nb / p2->cell.nneighbors;
            //     p1->Ddot_nonlocal += p2->Ddot_local * V_m;
            //     A += V_m;
            // }
            p1->Ddot_nonlocal /= A;
        }
    }

    bool any_damaged{false};
    for (Particle<nlayer> *pt : this->ass.pt_sys)
    {
        if (pt->type != undamaged_pt_type)
            any_damaged = pt->updateParticleFatigueDamage(dNdt) || any_damaged; // update the fatigue damage
    }

    return any_damaged;
}

template <int nlayer>
std::vector<LoadStep<nlayer>> SolverFatigue<nlayer>::generateLoadCycle(FatigueLoadType ltype, int N, int n_interval, std::vector<std::vector<Particle<nlayer> *>> pt_groups)
{
    // generate the load step for cycle-N
    double f_0{load_spectrum[2 * N]}, f_1{load_spectrum[2 * N + 1]}, f_2{load_spectrum[2 * N + 2]};
    std::vector<LoadStep<nlayer>> load_cycle;

    if (ltype == FatigueLoadType::LoadUniaxialDisp)
    {
        // loading part
        for (int i = 0; i < n_interval; ++i)
            load_cycle.push_back(LoadUniaxialDisp<nlayer>((f_1 - f_0) / n_interval, pt_groups[0], pt_groups[1], pt_groups[2]));

        // unloading part
        for (int i = 0; i < n_interval; ++i)
            load_cycle.push_back(LoadUniaxialDisp<nlayer>((f_2 - f_1) / n_interval, pt_groups[0], pt_groups[1], pt_groups[2]));
    }

    if (ltype == FatigueLoadType::LoadUniaxialForce)
    {
        // loading part
        for (int i = 0; i < n_interval; ++i)
            load_cycle.push_back(LoadUniaxialForce<nlayer>((f_1 - f_0) / n_interval, pt_groups[0], pt_groups[1], pt_groups[2]));

        // unloading part
        for (int i = 0; i < n_interval; ++i)
            load_cycle.push_back(LoadUniaxialForce<nlayer>((f_2 - f_1) / n_interval, pt_groups[0], pt_groups[1], pt_groups[2]));
    }

    if (ltype == FatigueLoadType::LoadDogBoneDisp)
    {
        // loading part
        for (int i = 0; i < n_interval; ++i)
            load_cycle.push_back(LoadDogBoneDisp<nlayer>(0, (f_1 - f_0) / n_interval, 0, pt_groups[0], pt_groups[1]));

        // unloading part
        for (int i = 0; i < n_interval; ++i)
            load_cycle.push_back(LoadDogBoneDisp<nlayer>(0, (f_2 - f_1) / n_interval, 0, pt_groups[0], pt_groups[1]));
    }

    if (ltype == FatigueLoadType::LoadCTForce)
    {
        // loading part
        for (int i = 0; i < n_interval; ++i)
            load_cycle.push_back(LoadCTForce<nlayer>((f_1 - f_0) / n_interval, pt_groups[0], pt_groups[1], pt_groups[2]));

        // unloading part
        for (int i = 0; i < n_interval; ++i)
            load_cycle.push_back(LoadCTForce<nlayer>((f_2 - f_1) / n_interval, pt_groups[0], pt_groups[1], pt_groups[2]));
    }

    return load_cycle;
}

template <int nlayer>
void SolverFatigue<nlayer>::solveProblemCyclic(FatigueLoadType loadType, int n_interval, std::vector<std::vector<Particle<nlayer> *>> pt_groups)
{
    // solve fatigue problem by dividing n intervals for a single loading step

    double n_cycles = 0.5 * (load_spectrum.size() - 1);
    double t{0}, dNdt{1}; // fake simulation time and derivative
    int N{0};             // real cycle number

    do
    {
        printf("Cycle-%d starts:\n\n", N + 1);

        // generate the load step for cycle-N
        std::vector<LoadStep<nlayer>> load_cycle = generateLoadCycle(loadType, N, n_interval, pt_groups);

        dNdt = (t_mode == TimeMapMode::Linear) ? func_dNdt_linear(tau) : func_dNdt_exponental(tau, t++);
        solveProblemOneCycle(load_cycle, dNdt);

        // compute average strain energy
        double ave_u{0};
        int n_pt = 0;
        for (Particle<nlayer> *pt : this->ass.pt_sys)
        {
            if (pt->type != undamaged_pt_type)
            {
                n_pt++;
                ave_u += pt->state_var[0];
            }
        }
        std::cout << "Average strain energy is: " << ave_u / n_pt << std::endl;
        // std::cout << this->ass.pt_sys[685]->damage << ',' << this->ass.pt_sys[685]->state_var[0] << std::endl;

        N += (int)dNdt;
        this->ass.writeDump(this->dumpFile, N);
    } while (N < n_cycles);
}

template <int nlayer>
void SolverFatigue<nlayer>::solveProblemOneCycle(std::vector<LoadStep<nlayer>> &load, double dNdt)
{
    bool new_damaged{false};

    int n_step = load.size();
    for (int i = 0; i < n_step; ++i)
    {
        printf("Substep-%d, iteration starts:\n", i);
        double t1 = omp_get_wtime();
        bool is_converged = solveProblemStep(load[i]);
        this->ass.updateStateVar(); // update the state variables in current cycle
        double t2 = omp_get_wtime();
        // this->ass.writeDump(this->dumpFile, 0);
        printf("Loading step %d has finished, spent %f seconds\n\n", i, t2 - t1);
    }

    new_damaged = updateFatigueDamage(dNdt); // delta_t is 1, so dN = dNdt * 1
    this->ass.storeStateVar();               // store converged state variables
    this->ass.updateGeometry();
    this->ass.updateForceState();
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

        this->ass.updateStateVar();
        this->ass.storeStateVar(); // store converged state variables

        double t2 = omp_get_wtime();
        printf("Loading step %d has finished, spent %f seconds\n\nData output ...\n\n", i + start_index, t2 - t1);

        this->ass.writeDump(this->dumpFile, i);
    }

    start_index += n_step;
}

template <int nlayer>
bool SolverFatigue<nlayer>::solveProblemStep(LoadStep<nlayer> &load_step)
{
    // keep the damage unchanged, update the deformation field
    this->updateForceBC(load_step);
    this->updateDisplacementBC(load_step);

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
    int n_newton = this->NewtonIteration();
    if (n_newton < this->max_NR_iter)
        return true; // normal return
    else
        return false; // abnormal return
}

#endif
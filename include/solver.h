#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include "lpm.h"
#include "load_step.h"
#include "stiffness.h"
#include "unit_cell.h"
#include "assembly.h"

template <int nlayer>
class Solver
{

public:
    int problem_size, max_NR_iter;
    double tol_NR_iter;
    std::string dumpFile;
    SolverMode sol_mode;
    Eigen::VectorXd disp;
    Eigen::VectorXd reaction_force;

    Stiffness<nlayer> stiffness;
    Assembly<nlayer> ass;

    void updateDisplacementBC(LoadStep<nlayer> &load_step);
    void updateForceBC(LoadStep<nlayer> &load_step);
    void updateRR(); // update residual force and reaction force
    void LPM_CG();

    void solveLinearSystem();
    int NewtonIteration(); // return number of Newton iterations

    Solver(Assembly<nlayer> &p_ass, const StiffnessMode &p_stiff_mode, const SolverMode &p_sol_mode, const std::string &p_dumpFile, const int &p_niter, const double &p_tol)
        : sol_mode{p_sol_mode},
          problem_size((p_ass.pt_sys[0]->cell.dim) * p_ass.pt_sys.size()),
          disp((p_ass.pt_sys[0]->cell.dim) * p_ass.pt_sys.size()),
          stiffness(p_ass.pt_sys, p_stiff_mode),
          dumpFile(p_dumpFile),
          ass(p_ass),
          max_NR_iter(p_niter),
          tol_NR_iter(p_tol),
          reaction_force(0)
    { // stiff_mode = 0 means finite difference; 1 means analytical
    }
};

template <int nlayer>
int Solver<nlayer>::NewtonIteration()
{
    // update the force state and residual
    ass.updateGeometry();
    ass.updateForceState();
    updateRR();

    // compute the Euclidean norm (L2 norm)
    double norm_residual = stiffness.residual.norm();
    double norm_reaction_force = reaction_force.norm();
    double tol_multiplier = std::max(norm_residual, norm_reaction_force);
    char tempChar1[] = "residual", tempChar2[] = "reaction";
    printf("|  Norm of residual is %.5e, norm of reaction is %.5e, tolerance criterion is based on ", norm_residual, norm_reaction_force);
    if (norm_residual > norm_reaction_force)
        printf("%s force\n", tempChar1);
    else
        printf("%s force\n", tempChar2);

    int ni{0};
    while (norm_residual > tol_NR_iter * tol_multiplier)
    {
        if (++ni > max_NR_iter)
            return max_NR_iter; // abnormal return

        printf("|  |  Iteration-%d: ", ni);
        solveLinearSystem(); // solve for the incremental displacement

        ass.updateGeometry();
        ass.updateForceState();
        updateRR(); /* update the RHS risidual force vector */
        norm_residual = stiffness.residual.norm();
        printf("|  |  Norm of residual is %.3e, residual ratio is %.3e\n", norm_residual, norm_residual / tol_multiplier);
    }

    return ni; // normal return, return number of iterations
}

template <int nlayer>
void Solver<nlayer>::solveLinearSystem()
{
    if (sol_mode == SolverMode::CG)
        LPM_CG();

    /* update the position */
    for (auto pt : ass.pt_sys)
    {
        std::array<double, NDIM> dxyz{0, 0, 0};
        for (int j = 0; j < pt->cell.dim; j++)
            dxyz[j] += disp((pt->cell.dim) * (pt->id) + j);

        pt->moveBy(dxyz);
    }
}

template <int nlayer>
void Solver<nlayer>::updateRR()
{
    std::vector<double> temp_force;
    for (Particle<nlayer> *pt : ass.pt_sys)
    {
        for (int k = 0; k < pt->cell.dim; k++)
        {
            /* residual force (exclude the DoF that being applied to displacement BC) */
            stiffness.residual((pt->cell.dim) * (pt->id) + k) = (1 - pt->disp_constraint[k]) * (pt->Pex[k] - pt->Pin[k]);

            /* reaction force (DoF being applied to displacement BC) */
            if (pt->disp_constraint[k] == 1)
                temp_force.push_back(pt->Pin[k]);
        }
    }

    reaction_force.resize(temp_force.size());
    for (size_t i = 0; i < temp_force.size(); ++i)
        reaction_force(i) = temp_force[i];
}

template <int nlayer>
void Solver<nlayer>::updateDisplacementBC(LoadStep<nlayer> &load_step)
{
    for (DispBC<nlayer> bc : load_step.dispBCs)
    {
        for (Particle<nlayer> *pt : bc.group)
        {
            if (bc.flag == 'x')
            {
                std::array<double, NDIM> dxyz{bc.step, 0, 0};
                if (bc.load_mode == LoadMode::Relative)
                    pt->moveBy(dxyz);
                else if (bc.load_mode == LoadMode::Absolute)
                    pt->moveTo(dxyz);
                pt->disp_constraint[0] = 1;
            }
            if (bc.flag == 'y')
            {
                std::array<double, NDIM> dxyz{0, bc.step, 0};
                if (bc.load_mode == LoadMode::Relative)
                    pt->moveBy(dxyz);
                else if (bc.load_mode == LoadMode::Absolute)
                    pt->moveTo(dxyz);
                pt->disp_constraint[1] = 1;
            }
            if (bc.flag == 'z')
            {
                std::array<double, NDIM> dxyz{0, 0, bc.step};
                if (bc.load_mode == LoadMode::Relative)
                    pt->moveBy(dxyz);
                else if (bc.load_mode == LoadMode::Absolute)
                    pt->moveTo(dxyz);
                pt->disp_constraint[2] = 1;
            }
        }
    }
}

template <int nlayer>
void Solver<nlayer>::updateForceBC(LoadStep<nlayer> &load_step)
{
    for (ForceBC<nlayer> bc : load_step.forceBCs)
    {
        int num_forceBC = (int)bc.group.size();
        for (Particle<nlayer> *pt : bc.group)
        {
            if (bc.load_mode == LoadMode::Relative)
            {
                pt->Pex[0] += bc.fx / num_forceBC;
                pt->Pex[1] += bc.fy / num_forceBC;
                pt->Pex[2] += bc.fz / num_forceBC;
            }
            else if (bc.load_mode == LoadMode::Absolute)
            {
                pt->Pex[0] = bc.fx / num_forceBC;
                pt->Pex[1] = bc.fy / num_forceBC;
                pt->Pex[2] = bc.fz / num_forceBC;
            }
        }
    }
}

template <int nlayer>
void Solver<nlayer>::LPM_CG()
{
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;

    cg.compute(stiffness.K_matrix);

    if (cg.info() != Eigen::Success)
    {
        // Decomposition failed
        std::cout << "Matrix decomposition failed." << std::endl;
        return;
    }

    disp = cg.solve(stiffness.residual);

    if (cg.info() != Eigen::Success)
    {
        // Solving failed
        std::cout << "Solving failed." << std::endl;
        return;
    }

    std::cout << "The system has been solved after " << cg.iterations() << " iterations." << std::endl;
    std::cout << "The estimated error is " << cg.error() << std::endl;
}

#endif
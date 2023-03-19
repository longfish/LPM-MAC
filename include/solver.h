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
    const int max_iter = 20;      /* maximum global iteration number */
    const double tol_iter = 1e-5; /* newton iteration tolerance */

    int problem_size;
    double *disp;
    std::string dumpFile;
    std::vector<double> reaction_force;
    SolverMode sol_mode;

public:
    Stiffness<nlayer> stiffness;

    void updateDisplacementBC(Assembly<nlayer> &ass, LoadStep<nlayer> &load_step);
    void updateForceBC(Assembly<nlayer> &ass, LoadStep<nlayer> &load_step);
    void updateRR(Assembly<nlayer> &ass); // update residual force and reaction force
    void LPM_PARDISO();
    void LPM_CG();

    int NewtonIteration(Assembly<nlayer> &ass); // return number of Newton iterations

    Solver(Assembly<nlayer> &ass, const StiffnessMode &p_stiff_mode, const SolverMode &p_sol_mode, const std::string &p_dumpFile)
        : sol_mode{p_sol_mode}, stiffness(ass.pt_sys, p_stiff_mode), dumpFile(p_dumpFile)
    { // stiff_mode = 0 means finite difference; 1 means analytical
        problem_size = (ass.pt_sys[0]->cell.dim) * ass.pt_sys.size();
        disp = new double[problem_size];
        dumpFile = p_dumpFile;
    }

    ~Solver()
    {
        delete[] disp;
    }

    void solveProblem(Assembly<nlayer> &ass, std::vector<LoadStep<nlayer>> &load);
    void solveStepwise(Assembly<nlayer> &ass);
};

template <int nlayer>
void Solver<nlayer>::solveStepwise(Assembly<nlayer> &ass)
{
    if (sol_mode == SolverMode::PARDISO)
        LPM_PARDISO();
    else
        LPM_CG();
    
    /* update the position */
    for (auto pt : ass.pt_sys)
    {
        std::array<double, NDIM> dxyz{0, 0, 0};
        for (int j = 0; j < pt->cell.dim; j++)
            dxyz[j] += disp[(pt->cell.dim) * (pt->id) + j];

        pt->moveBy(dxyz);
    }
}

template <int nlayer>
void Solver<nlayer>::solveProblem(Assembly<nlayer> &ass, std::vector<LoadStep<nlayer>> &load)
{
    ass.writeDump(dumpFile, 0);

    for (int i = 0; i < load.size(); i++)
    {
        printf("Loading step-%d: ", i + 1);

        double t1 = omp_get_wtime();
        if (ass.pt_sys[0]->cell.dim == 2)
            stiffness.initializeStiffness2D(ass.pt_sys);
        else
            stiffness.initializeStiffness3D(ass.pt_sys);
        double t2 = omp_get_wtime();
        printf("Stiffness matrix calculation costs %f seconds\n", t2 - t1);

        updateDisplacementBC(ass, load[i]);
        updateForceBC(ass, load[i]);
        ass.updateForceState();
        updateRR(ass);

        int ni = NewtonIteration(ass);
        ass.writeDump(dumpFile, i);

        printf("Loading step %d has finished in %d iterations\n\nData output ...\n", i + 1, ni);
    }
}

template <int nlayer>
int Solver<nlayer>::NewtonIteration(Assembly<nlayer> &ass)
{
    // compute the Euclidean norm (L2 norm)
    double norm_residual = cblas_dnrm2(problem_size, stiffness.residual, 1);
    double norm_reaction_force = cblas_dnrm2(reaction_force.size(), reaction_force.data(), 1);
    double tol_multiplier = MAX(norm_residual, norm_reaction_force);
    char tempChar1[] = "residual", tempChar2[] = "reaction";
    printf("Norm of residual is %.5e, norm of reaction is %.5e, tolerance criterion is based on ", norm_residual, norm_reaction_force);
    if (norm_residual > norm_reaction_force)
        printf("%s force\n", tempChar1);
    else
        printf("%s force\n", tempChar2);

    int ni{0};
    while (norm_residual > tol_iter * tol_multiplier && ni < max_iter)
    {
        printf("    Iteration-%d: ", ++ni);
        // double t1 = omp_get_wtime();
        // if (ass.pt_sys[0]->cell.dim == 2)
        //     stiffness.initializeStiffness2D(ass.pt_sys);
        // else
        //     stiffness.initializeStiffness3D(ass.pt_sys);
        // double t2 = omp_get_wtime();
        // printf("Stiffness matrix calculation costs %f seconds\n", t2 - t1);

        stiffness.updateStiffnessDispBC(ass.pt_sys);
        solveStepwise(ass); // solve for the incremental displacement
        ass.updateForceState();
        updateRR(ass); /* update the RHS risidual force vector */
        norm_residual = cblas_dnrm2(problem_size, stiffness.residual, 1);
        printf("    Norm of residual is %.3e, residual ratio is %.3e\n", norm_residual, norm_residual / tol_multiplier);
    }

    return ni;
}

template <int nlayer>
void Solver<nlayer>::updateRR(Assembly<nlayer> &ass)
{
    for (Particle<nlayer> *pt : ass.pt_sys)
    {
        for (int k = 0; k < pt->cell.dim; k++)
        {
            /* residual force (exclude the DoF that being applied to displacement BC) */
            stiffness.residual[(pt->cell.dim) * (pt->id) + k] = (1 - pt->disp_constraint[k]) * (pt->Pex[k] - pt->Pin[k]);

            /* reaction force (DoF being applied to displacement BC) */
            if (pt->disp_constraint[k] == 1)
                reaction_force.push_back(pt->Pin[k]);
        }
    }
}

template <int nlayer>
void Solver<nlayer>::updateDisplacementBC(Assembly<nlayer> &ass, LoadStep<nlayer> &load_step)
{
    for (DispBC<nlayer> bc : load_step.dispBCs)
    {
        for (Particle<nlayer> *pt : bc.group)
        {
            if (bc.flag == 'x')
            {
                std::array<double, NDIM> dxyz{bc.step, 0, 0};
                pt->moveBy(dxyz);
                pt->disp_constraint[0] = 1;
            }
            if (bc.flag == 'y')
            {
                std::array<double, NDIM> dxyz{0, bc.step, 0};
                pt->moveBy(dxyz);
                pt->disp_constraint[1] = 1;
            }
            if (bc.flag == 'z')
            {
                std::array<double, NDIM> dxyz{0, 0, bc.step};
                pt->moveBy(dxyz);
                pt->disp_constraint[2] = 1;
            }
        }
    }
}

template <int nlayer>
void Solver<nlayer>::updateForceBC(Assembly<nlayer> &ass, LoadStep<nlayer> &load_step)
{
    for (ForceBC<nlayer> bc : load_step.forceBCs)
    {
        int num_forceBC = (int)bc.group.size();
        for (Particle<nlayer> *pt : bc.group)
        {
            pt->Pex[0] += bc.fx / num_forceBC;
            pt->Pex[1] += bc.fy / num_forceBC;
            pt->Pex[2] += bc.fz / num_forceBC;
        }
    }
}

template <int nlayer>
void Solver<nlayer>::LPM_PARDISO()
{
    MKL_INT n, idum, maxfct, mnum, mtype, phase, error, error1, msglvl, nrhs, iter;
    MKL_INT iparm[64];
    double ddum;

    void *pt[64]; /* Internal solver memory pointer */

    for (int i = 0; i < 64; i++)
        iparm[i] = 0;
    iparm[0] = 1;  /* No solver default */
    iparm[1] = 3;  /* The parallel (OpenMP) version of the nested dissection algorithm is used */
    iparm[3] = 0;  /* No iterative-direct algorithm */
    iparm[4] = 0;  /* No user fill-in reducing permutation */
    iparm[5] = 0;  /* Write solution into x */
    iparm[6] = 0;  /* Not in use */
    iparm[7] = 0;  /* Max numbers of iterative refinement steps */
    iparm[8] = 0;  /* Not in use */
    iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
    // iparm[11] = 0;        /* Not in use */
    iparm[12] = 0; /* Maximum weighted matching algorithm is switched-off (default for
                          symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
                   // iparm[13] = 0;        /* Output: Number of perturbed pivots */
    // iparm[14] = 0;        /* Not in use */
    // iparm[15] = 0;        /* Not in use */
    // iparm[16] = 0;        /* Not in use */
    // iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    // iparm[18] = -1;       /* Output: Mflops for LU factorization */
    // iparm[19] = 0;        /* Output: Numbers of CG Iterations */

    n = problem_size; /* Data number */
    maxfct = 1;       /* Maximum number of numerical factorizations */
    mnum = 1;         /* Which factorization to use */
    msglvl = 0;       /* 0, no print statistical info; 1, print statistical info */
    error = 0;        /* Initialize error flag */

    mtype = 2; /* Real symmetric positive definite, -2: real+symmetric+indefinite */
    nrhs = 1;  /* Number of right hand sides */
    iter = 1;  /* Iteration number */

    for (int i = 0; i < 64; i++)
        pt[i] = 0; /* Initiliaze the internal solver memory pointer */

    /* Reordering and Symbolic Factorization. This step also allocates all memory that is  */
    /* necessary for the factorization */
    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, stiffness.K_global, stiffness.IK, stiffness.JK, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
    {
        printf("\nERROR during symbolic factorization: " IFORMAT, error);
        exit(1);
    }
    // printf("PARDISO: Size of factors(MB): %f", MAX(iparm[14], iparm[15] + iparm[16]) / 1000.0);
    printf("    PARDISO: Size of factors(MB): %f", iparm[16] / 1000.0);

    /* Numerical factorization */
    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, stiffness.K_global, stiffness.IK, stiffness.JK, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
    {
        printf("\nERROR during numerical factorization: " IFORMAT, error);
        exit(2);
    }

    /* Back substitution and iterative refinement */
    phase = 33;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, stiffness.K_global, stiffness.IK, stiffness.JK, &idum, &nrhs, iparm, &msglvl, stiffness.residual, disp, &error);
    if (error != 0)
    {
        printf("\nERROR during solution: " IFORMAT, error);
        exit(3);
    }
    printf("\n    Solve completed at iteration: " IFORMAT "\n", iter);

    /* Termination and release of memory */
    phase = -1; /* Release internal memory. */
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, stiffness.IK, stiffness.JK, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error1);
    if (error1 != 0)
    {
        printf("\nERROR on release stage: " IFORMAT, error1);
        exit(4);
    }
}

template <int nlayer>
void Solver<nlayer>::LPM_CG()
{
    MKL_INT n, rci_request, itercount, mkl_disable_fast_mm;
    MKL_INT ipar[128];
    double dpar[128], *tmp;

    /* matrix descriptor */
    struct matrix_descr descrA;
    sparse_matrix_t csrA;
    sparse_operation_t transA = SPARSE_OPERATION_NON_TRANSPOSE;
    descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;
    mkl_disable_fast_mm = 1; /* avoid memory leaks */

    /* initial setting */
    n = problem_size; /* Data number */
    tmp = new double[4 * n];
    mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ONE, n, n, stiffness.IK, stiffness.IK + 1, stiffness.JK, stiffness.K_global);

    /* initial guess for the displacement vector */
    for (int i = 0; i < n; i++)
        disp[i] = 0;

    /* initialize the solver */
    dcg_init(&n, disp, stiffness.residual, &rci_request, ipar, dpar, tmp);
    if (rci_request != 0)
        goto failure;

    /* modify the initialized solver parameters */
    ipar[4] = problem_size;
    ipar[8] = 1; /* default value is 0, does not perform the residual stopping test; otherwise, perform the test */
    ipar[9] = 0; /* default value is 1, perform user defined stopping test; otherwise, does not perform the test */
    // ipar[10] = 1; /* use the preconditioned version of the CG method */
    dpar[0] = 1e-12;  /* specifies the relative tolerance, the default value is 1e-6 */
    dpar[1] = 1e-12; /* specifies the absolute tolerance, the default value is 0.0 */

    /* check the correctness and consistency of the newly set parameters */
    dcg_check(&n, disp, stiffness.residual, &rci_request, ipar, dpar, tmp);
    if (rci_request != 0)
        goto failure;

    /* compute the solution by RCI (residual)CG solver */
    /* reverse Communications starts here */
rci:
    dcg(&n, disp, stiffness.residual, &rci_request, ipar, dpar, tmp);
    /* if rci_request=0, then the solution was found according to the requested  */
    /* stopping tests. in this case, this means that it was found after 100      */
    /* iterations. */
    if (rci_request == 0)
        goto getsln;

    /* If rci_request=1, then compute the vector K_global*TMP[0]                  */
    /* and put the result in vector TMP[n]                                       */
    if (rci_request == 1)
    {
        mkl_sparse_d_mv(transA, 1.0, csrA, descrA, tmp, 0.0, &tmp[n]);
        goto rci;
    }

    /* If rci_request=anything else, then dcg subroutine failed                  */
    /* to compute the solution vector: solution[n]                               */
    goto failure;
    /* Reverse Communication ends here                                           */
    /* Get the current iteration number into itercount                           */
getsln:
    dcg_get(&n, disp, stiffness.residual, &rci_request, ipar, dpar, tmp, &itercount);
    printf("    The system has been solved after " IFORMAT " iterations\n", itercount);
    goto success;

failure:
    printf("    The computation FAILED as the solver has returned the ERROR code " IFORMAT "\n", rci_request);

success:
    /* free memory */
    delete[] tmp;
}

#endif
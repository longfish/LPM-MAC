#include "lpm.h"
#include "particle.h"
#include "assembly.h"
#include "utilities.h"
#include "stiffness.h"
#include "load_step.h"
#include "solver_static.h"

void run()
{
    double start = omp_get_wtime(); // record the CPU time, begin

    const int n_layer = 2; // number of neighbor layers (currently only support 2 layers of neighbors)
    double radius = 0.6;   // particle radius
    UnitCell cell(LatticeType::SimpleCubic3D, radius);

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{-0.0, 10.0, -0.4, 10.0, -0.4, 30.0};

    std::vector<std::array<double, NDIM>> sc_xyz = createCuboidSC3D(box, cell, R_matrix);
    Assembly<n_layer> pt_ass{sc_xyz, box, cell, ParticleType::Elastic}; // elastic bond
    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    double E0 = 69e3, mu0 = 0.3; // Young's modulus and Poisson's ratio
    std::vector<Particle<n_layer> *> top_group, bottom_group, internal_group;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        // assign boundary and internal particles
        if (p1->xyz[2] > box[5] - 1.5 * radius)
        {
            top_group.push_back(p1); // top
            p1->type = 1;
        }
        if (p1->xyz[2] < box[4] + 2 * radius)
        {
            bottom_group.push_back(p1); // bottom
            p1->type = 2;
        }
        if (p1->nb == cell.nneighbors)
            internal_group.push_back(p1); // particles with full neighbor list

        // assign material properties - need to cast to elastic particle
        ParticleElastic<n_layer> *elpt = dynamic_cast<ParticleElastic<n_layer> *>(p1);
        elpt->setParticleProperty(E0, mu0);
    }

    // simulation settings
    int n_steps = 1;          // number of loading steps
    double step_size = -1e-3; // step size for force or displacement loading

    std::vector<LoadStep<n_layer>> load; // load settings for multiple steps
    for (int i = 0; i < n_steps; i++)
    {
        LoadStep<n_layer> step; 

        // boundary conditions
        step.dispBCs.push_back(DispBC<n_layer>(top_group, 'x', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_group, 'y', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_group, 'z', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(bottom_group, 'z', step_size));
        load.push_back(step);
    }

    pt_ass.updateGeometry();
    pt_ass.updateForceState();

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    int max_iter = 1;                                                                                                         /* maximum Newton iteration number */
    double tol_iter = 1e-5;                                                                                                    /* newton iteration tolerance */
    SolverStatic<n_layer> solv{pt_ass, StiffnessMode::Analytical, SolverMode::PARDISO, "result_position.dump", max_iter, tol_iter}; // stiffness mode and solution mode
    solv.solveProblem(load);

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
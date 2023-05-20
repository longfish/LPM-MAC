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
    double radius = 0.2;   // particle radius
    UnitCell cell(LatticeType::Hexagon2D, radius);

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{0.0, 40.0, 0.0, 40.0, 0.0, 8.0};                                                                // thickness is used for force calculation
    Assembly<n_layer> pt_ass{"../geometry/geo1_CT_2DHEX.dump", "../geometry/geo1_CT_2DHEX.bond", cell, ParticleType::ElasticDamage}; // read coordinate from local files

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    double E0 = 205e3, mu0 = 0.29;                  // Steel, Young's modulus and Poisson's ratio, MPa
    double alpha = 0.96, beta = 100, kappa0 = 6e-3; // brittle damage parameters
    double comp_tensile_ratio = 0.9;

    std::vector<Particle<n_layer> *> top_group, bottom_group, mid_group;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        // assign boundary and internal particles
        if (p1->id == 5030)
        {
            mid_group.push_back(p1); // mid, to fix
            p1->type = 1;
        }
        if (p1->id == 8614)
        {
            top_group.push_back(p1); // top
            p1->type = 2;
        }
        if (p1->id == 1273)
        {
            bottom_group.push_back(p1); // bottom
            p1->type = 3;
        }

        // assign material properties - need to cast to elastic particle
        ParticleElasticDamage<n_layer> *elpt = dynamic_cast<ParticleElasticDamage<n_layer> *>(p1);
        elpt->setParticleProperty(E0, mu0, kappa0, alpha, beta, comp_tensile_ratio);
    }

    // simulation settings
    int n_steps = 80;         // number of loading steps
    double step_size = -2e-4; // step size for displacement loading
    // double step_size = -600; // step size for force loading

    std::vector<LoadStep<n_layer>> load; // load settings for multiple steps

    for (int i = 0; i < n_steps; i++)
    {
        LoadStep<n_layer> step;

        // boundary conditions
        step.dispBCs.push_back(DispBC<n_layer>(mid_group, LoadMode::Relative, 'x', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(mid_group, LoadMode::Relative, 'y', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_group, LoadMode::Relative, 'y', -step_size));
        step.dispBCs.push_back(DispBC<n_layer>(bottom_group, LoadMode::Relative, 'y', step_size));

        // below is used for force control loading
        // step.dispBCs.push_back(DispBC<n_layer>(top_group, 'x', 0.0));
        // step.dispBCs.push_back(DispBC<n_layer>(bottom_group, 'x', 0.0));
        // step.forceBCs.push_back(ForceBC<n_layer>(top_group, 0.0, -step_size, 0.0));
        // step.forceBCs.push_back(ForceBC<n_layer>(bottom_group, 0.0, step_size, 0.0));
        load.push_back(step);
    }
    load[0].dispBCs[2].step *= 50; // increase the elastic loading step size
    load[0].dispBCs[3].step *= 50;

    pt_ass.updateGeometry();
    pt_ass.updateForceState();
    pt_ass.updateStateVar();

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    int max_iter = 30, start_index = 0;                                                                                          /* maximum Newton iteration number */
    double tol_iter = 1e-5;                                                                                                      /* newton iteration tolerance */
    SolverStatic<n_layer> solv{pt_ass, StiffnessMode::Analytical, SolverMode::CG, "CT_2DHex_position.dump", max_iter, tol_iter}; // stiffness mode and solution mode
    solv.solveProblem(load, start_index);

    // output top loading point's reaction force
    // std::cout << 8613 << ',' << pt_ass.pt_sys[8613]->Pin[0] << ',' << pt_ass.pt_sys[8613]->Pin[1] << ',' << pt_ass.pt_sys[8613]->Pin[2] << std::endl;
    // std::cout << 8614 << ',' << pt_ass.pt_sys[8614]->Pin[0] << ',' << pt_ass.pt_sys[8614]->Pin[1] << ',' << pt_ass.pt_sys[8614]->Pin[2] << std::endl;
    // std::cout << 8615 << ',' << pt_ass.pt_sys[8615]->Pin[0] << ',' << pt_ass.pt_sys[8615]->Pin[1] << ',' << pt_ass.pt_sys[8615]->Pin[2] << std::endl;

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
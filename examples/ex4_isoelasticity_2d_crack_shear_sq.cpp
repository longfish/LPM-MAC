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
    double radius = 5e-3;  // particle radius
    UnitCell cell(LatticeType::Square2D, radius);

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    Eigen::MatrixXd R_matrix = createRMatrix(eulerflag, angles);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{0.0, 1.0, 0.0, 1.0, 0.0, 8.0};                                                                      // thickness is used for force calculation
    Assembly<n_layer> pt_ass{"../geometry/geo1_2DSQ_Shear.dump", "../geometry/geo1_2DSQ_Shear.bond", cell, ParticleType::ElasticDamage}; // read coordinate from local files

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    bool is_plane_stress = true;
    double E0 = 210e3, mu0 = 0.3;                                // polymer, Young's modulus and Poisson's ratio, MPa
    double k0 = 50, k1 = 5, c_t_ratio = 10, damage_thres = 0.99; // brittle damage parameters

    // simulation settings
    int n_steps = 100;                           // number of loading steps
    double cutoff_ratio = 1.5;                   // nonlocal cutoff ratio
    double step_size = 1.5e-4;                   // step size for displacement loading
    double nonlocal_L = 2e-2;                    // nonlocal length scale
    int undamaged_pt_type = 4;                   // particles that dont update damage
    int max_iter = 30, start_index = 0;          // maximum Newton iteration number
    double tol_iter = 1e-5;                      // newton iteration tolerance
    std::string dumpFile{"Shear_2d_crack_sq.dump"}; // output file name

    std::vector<Particle<n_layer> *> top_group, bottom_group;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        // assign boundary and internal particles
        if (p1->xyz[1] > 0.99)
        {
            top_group.push_back(p1); // top
            p1->type = 1;
        }
        if (p1->xyz[1] < 0.01)
        {
            bottom_group.push_back(p1); // bottom
            p1->type = 2;
        }
        // if (p1->xyz[0] >= 29.5 - 9.5 / 2 && p1->xyz[0] <= 29.5 + 9.5 / 2 && p1->xyz[1] >= 40 - 9.2)
        //     p1->type = 4;
        // if (p1->xyz[0] >= 29.5 - 9.5 / 2 && p1->xyz[0] <= 29.5 + 9.5 / 2 && p1->xyz[1] <= 9.2)
        //     p1->type = 4; // particles that not update damage

        // assign material properties - need to cast to elastic particle
        ParticleElasticDamage<n_layer> *elpt = dynamic_cast<ParticleElasticDamage<n_layer> *>(p1);
        elpt->setParticleProperty(nonlocal_L, is_plane_stress, E0, mu0, k0, k1, c_t_ratio, damage_thres);
    }

    std::vector<LoadStep<n_layer>> load; // load settings for multiple steps
    for (int i = 0; i < n_steps; i++)
    {
        LoadStep<n_layer> step;

        // boundary conditions
        step.dispBCs.push_back(DispBC<n_layer>(bottom_group, LoadMode::Relative, 'x', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(bottom_group, LoadMode::Relative, 'y', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_group, LoadMode::Relative, 'x', step_size));
        step.dispBCs.push_back(DispBC<n_layer>(top_group, LoadMode::Relative, 'y', 0));

        // below is used for force control loading
        // step.dispBCs.push_back(DispBC<n_layer>(top_group, 'x', 0.0));
        // step.dispBCs.push_back(DispBC<n_layer>(bottom_group, 'x', 0.0));
        // step.forceBCs.push_back(ForceBC<n_layer>(top_group, 0.0, -step_size, 0.0));
        // step.forceBCs.push_back(ForceBC<n_layer>(bottom_group, 0.0, step_size, 0.0));
        load.push_back(step);
    }
    // load[0].dispBCs[2].step *= 20; // increase the elastic loading step size
    // load[0].dispBCs[3].step *= 20;

    pt_ass.searchNonlocalNeighbors(cutoff_ratio);
    pt_ass.updateGeometry();
    pt_ass.updateForceState();
    pt_ass.updateStateVar();

    SolverStatic<n_layer> solv{undamaged_pt_type, pt_ass, StiffnessMode::Analytical, SolverMode::CG, dumpFile, max_iter, tol_iter}; // stiffness mode and solution mode

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    solv.solveProblem(load, start_index);

    // output top loading point's reaction force
    // std::cout << 8613 << ',' << pt_ass.pt_sys[8613]->Pin[0] << ',' << pt_ass.pt_sys[8613]->Pin[1] << ',' << pt_ass.pt_sys[8613]->Pin[2] << std::endl;
    // std::cout << 8614 << ',' << pt_ass.pt_sys[8614]->Pin[0] << ',' << pt_ass.pt_sys[8614]->Pin[1] << ',' << pt_ass.pt_sys[8614]->Pin[2] << std::endl;
    // std::cout << 8615 << ',' << pt_ass.pt_sys[8615]->Pin[0] << ',' << pt_ass.pt_sys[8615]->Pin[1] << ',' << pt_ass.pt_sys[8615]->Pin[2] << std::endl;

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
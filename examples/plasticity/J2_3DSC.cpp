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
    double radius = 0.25;  // particle radius
    UnitCell cell(LatticeType::SimpleCubic3D, radius);

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{-0.0, 10.0, -0.0, 10.0, -0.0, 10.0};

    std::vector<std::array<double, NDIM>> sc_xyz = createCuboidSC3D(box, cell, R_matrix);
    Assembly<n_layer> pt_ass{sc_xyz, box, cell, ParticleType::J2Plasticity}; // J2 plasticity particle
    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    bool is_plane_stress = false;
    double E0 = 146e3, mu0 = 0.3;                      // Young's modulus and Poisson's ratio
    double sigmay = 200, J2_xi = 0.0, J2_H = 38.714e3; // initial yield stress, isotropic hardening, hardening modulus
    double A = 10;                                     // plasticity damage parameter
    double critical_bstrain = 1;                       // critical bond strain

    // simulation settings
    int n_steps = 91;                      // number of loading steps
    double step_size = -2000;              // step size for force loading
    double nonlocal_L = 0;                 // nonlocal length scale
    int max_iter = 30, start_index = 0;    // maximum Newton iteration number
    double tol_iter = 1e-5;                // newton iteration tolerance
    int undamaged_pt_type = -1;            // undamaged particle type
    std::string dumpFile{"sc_j2_3d.dump"}; // output file name

    std::vector<Particle<n_layer> *> top_group, top_edge_y0, top_edge_x0, bottom_group, internal_group;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        // assign boundary and internal particles
        if (p1->xyz[2] > box[5] - 2 * radius)
        {
            top_group.push_back(p1); // top
            p1->type = 1;
            if (p1->xyz[1] < 0.2)
            {
                top_edge_y0.push_back(p1); // edge, y = 0
                p1->type = 3;
            }
            if (p1->xyz[0] < 0.2)
            {
                top_edge_x0.push_back(p1); // edge, x = 0
                p1->type = 4;
            }
        }
        if (p1->xyz[2] < box[4] + 2 * radius)
        {
            bottom_group.push_back(p1); // bottom
            p1->type = 2;
        }
        if (p1->nb == cell.nneighbors)
            internal_group.push_back(p1); // particles with full neighbor list

        // assign material properties - need to cast to elastic particle
        ParticleJ2Plasticity<n_layer> *elpt = dynamic_cast<ParticleJ2Plasticity<n_layer> *>(p1);
        elpt->setParticleProperty(nonlocal_L, is_plane_stress, E0, mu0, sigmay, J2_xi, J2_H, A, critical_bstrain);
    }

    std::vector<LoadStep<n_layer>> load; // load settings for multiple steps
    for (int i = 0; i < 14; i++)
    {
        LoadStep<n_layer> step;

        // boundary conditions
        step.dispBCs.push_back(DispBC<n_layer>(top_group, LoadMode::Relative, 'z', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_edge_x0, LoadMode::Relative, 'x', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_edge_y0, LoadMode::Relative, 'y', 0.0));
        step.forceBCs.push_back(ForceBC<n_layer>(bottom_group, LoadMode::Relative, 0.0, 0.0, step_size));
        load.push_back(step);
    }
    for (int i = 14; i < 48; i++)
    {
        LoadStep<n_layer> step;

        // boundary conditions
        step.dispBCs.push_back(DispBC<n_layer>(top_group, LoadMode::Relative, 'z', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_edge_x0, LoadMode::Relative, 'x', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_edge_y0, LoadMode::Relative, 'y', 0.0));
        step.forceBCs.push_back(ForceBC<n_layer>(bottom_group, LoadMode::Relative, 0.0, 0.0, -step_size));
        load.push_back(step);
    }
    for (int i = 48; i < n_steps; i++)
    {
        LoadStep<n_layer> step;

        // boundary conditions
        step.dispBCs.push_back(DispBC<n_layer>(top_group, LoadMode::Relative, 'z', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_edge_x0, LoadMode::Relative, 'x', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_edge_y0, LoadMode::Relative, 'y', 0.0));
        step.forceBCs.push_back(ForceBC<n_layer>(bottom_group, LoadMode::Relative, 0.0, 0.0, step_size));
        load.push_back(step);
    }

    // pt_ass.searchNonlocalNeighbors(cutoff_ratio);
    pt_ass.updateGeometry();
    pt_ass.updateForceState();

    SolverStatic<n_layer> solv{undamaged_pt_type, pt_ass, StiffnessMode::Analytical, SolverMode::CG, dumpFile, max_iter, tol_iter}; // stiffness mode and solution mode

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    // write down global matrices
    solv.solveProblem(load, start_index);

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
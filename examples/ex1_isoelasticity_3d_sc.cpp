#include "lpm.h"
#include "particle.h"
#include "assembly.h"
#include "utilities.h"
#include "stiffness.h"
#include "load_step.h"
#include "solver.h"

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
    std::array<double, 2 * NDIM> box{-0.0, 10.0, -0.4, 10.0, -0.4, 30.0};

    std::vector<std::array<double, NDIM>> sc_xyz = createCuboidSC3D(box, cell, R_matrix);
    Assembly<n_layer> pt_ass{sc_xyz, box, cell, BondType::Elastic}; // elastic bond with brittle damage law
    // pt_ass.writeBond("../geometry/test.bond");
    // pt_ass.writeConfigurationDump("../geometry/test.dump");
    // Assembly<n_layer> pt_ass{"../geometry/test.dump", "../geometry/test.bond", cell, BondType::Elastic}; // read coordinate from local files
    // Assembly<n_layer> pt_ass{"../geometry/test.dump", cell, BondType::Elastic}; // read coordinate from local files

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    double E0 = 69e3, mu0 = 0.3;      // Young's modulus and Poisson's ratio
    double critical_bstrain = 1.0e-2; // critical bond strain value at which bond will break
    int nbreak = 20;                  // limit the broken number of bonds in a single iteration, should be an even number

    std::vector<Particle<n_layer> *> top_group, bottom_group, internal_group;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        // assign boundary and internal particles
        if (p1->xyz[2] > box[5] - 1.5 * radius)
            top_group.push_back(p1); // top
        if (p1->xyz[2] < box[4] + 1.5 * radius)
            bottom_group.push_back(p1); // bottom
        if (p1->nb == cell.nneighbors)
            internal_group.push_back(p1); // particles with full neighbor list

        // assign material properties
        for (int i = 0; i < n_layer; ++i)
        {
            for (auto bd : p1->bond_layers[i])
            {
                // cast to elastic bond (or other type of bonds)
                BondElastic<n_layer> *elbd = dynamic_cast<BondElastic<n_layer> *>(bd);
                elbd->setBondProperty(E0, mu0, critical_bstrain, nbreak);
            }
        }
    }

    // simulation settings
    int n_steps = 1;          // number of loading steps
    double step_size = -1e-3; // step size for force or displacement loading

    std::vector<LoadStep<n_layer>> load; // load settings for multiple steps
    for (int i = 0; i < n_steps; i++)
    {
        LoadStep<n_layer> step{1}; // 1 means tension loading, -1 means compression loading

        // boundary conditions
        step.dispBCs.push_back(DispBC<n_layer>(top_group, 'x', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_group, 'y', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_group, 'z', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(bottom_group, 'z', step_size));
        load.push_back(step);
    }

    pt_ass.updateForceState();

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    Solver<n_layer> solv{pt_ass, StiffnessMode::Analytical, SolverMode::CG, "result_position.dump"}; // stiffness mode and solution mode
    solv.solveProblem(pt_ass, load);

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
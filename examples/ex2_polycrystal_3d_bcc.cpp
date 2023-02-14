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
    double radius = 0.02;  // 1.15; // particle radius
    UnitCell cell(LatticeType::BCC3D_1Slip, radius);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{0, 1, 0, 1, 0, 1};

    // std::vector<std::array<double, NDIM>> fcc_xyz = createCuboidFCC3D(box, cell, R_matrix);
    // Assembly<n_layer> pt_ass{fcc_xyz, box, cell, BondType::Elastic}; // elastic bond with brittle damage law
    Assembly<n_layer> pt_ass{"../geometry/BCC_Hailong/BCC_10grains.dump", cell, BondType::Elastic}; // read coordinate from dump file

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    double C11{230e3}, C12{135e3}, C44{117e3}; // Elastic constants
    double critical_bstrain = 1.0e-2;          // critical bond strain value at which bond will break
    int nbreak = 20;                           // limit the broken number of bonds in a single iteration, should be an even number

    std::vector<Particle<n_layer> *> top, bottom, left, right, front, back, internal;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        // assign boundary and internal particles
        if (p1->xyz[2] > 1 - 2 * radius)
            top.push_back(p1);
        if (p1->xyz[2] < 2 * radius)
            bottom.push_back(p1);
        if (p1->xyz[1] > 1 - 2 * radius)
            back.push_back(p1);
        if (p1->xyz[1] < 2 * radius)
            front.push_back(p1);
        if (p1->xyz[0] > 1 - 2 * radius)
            right.push_back(p1);
        if (p1->xyz[0] < 2 * radius)
            left.push_back(p1);
        if (p1->nb == cell.nneighbors)
            internal.push_back(p1); // particles with full neighbor list

        // assign material properties
        for (int i = 0; i < n_layer; ++i)
        {
            for (auto bd : p1->bond_layers[i])
            {
                // cast to elastic bond (or other type of bonds)
                ElasticBond<n_layer> *elbd = dynamic_cast<ElasticBond<n_layer> *>(bd);
                elbd->setBondProperty(C11, C12, C44, critical_bstrain, nbreak);
            }
        }
    }

    // simulation settings
    int n_steps = 1;                                    // number of loading steps
    double U_stepx{1e-3}, U_stepy{1e-3}, U_stepz{1e-3}; // step size
    std::vector<LoadStep<n_layer>> load;                // load settings for multiple steps
    for (int i = 0; i < n_steps; i++)
    {
        LoadStep<n_layer> step{1}; // 1 means tension loading, -1 means compression loading

        // boundary conditions
        // left: U_stepx; right: -U_stepx
        // front: U_stepy; back: -U_stepy
        // top: U_stepz; bottom: -U_stepz

        step.dispBCs.push_back(DispBC<n_layer>(left, 'x', -U_stepx));
        step.dispBCs.push_back(DispBC<n_layer>(right, 'x', U_stepx));
        step.dispBCs.push_back(DispBC<n_layer>(front, 'y', -U_stepy));
        step.dispBCs.push_back(DispBC<n_layer>(back, 'y', U_stepy));
        step.dispBCs.push_back(DispBC<n_layer>(top, 'z', U_stepz));
        step.dispBCs.push_back(DispBC<n_layer>(bottom, 'z', -U_stepz));
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
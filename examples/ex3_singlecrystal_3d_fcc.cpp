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

    const int n_layer = 2;             // number of neighbor layers (currently only support 2 layers of neighbors)
    double radius = 0.244948964611159; // 1.15; // particle radius
    UnitCell cell(LatticeType::FCC3D, radius);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{0, 10, 0, 10, 0, 10};

    // std::vector<std::array<double, NDIM>> fcc_xyz = createCuboidFCC3D(box, cell, R_matrix);
    // Assembly<n_layer> pt_ass{fcc_xyz, box, cell, BondType::Elastic}; // elastic bond with brittle damage law
    Assembly<n_layer> pt_ass{"../geometry/FCC_Hailong/FCC_single.dump", cell, BondType::Elastic}; // read coordinate from dump file

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    double C11{108e3}, C12{61e3}, C44{29e3}; // Elastic constants
    //double E0 = 64e3, mu0 = 0.36;      // Young's modulus and Poisson's ratio
    double critical_bstrain = 1.0e-2; // critical bond strain value at which bond will break
    int nbreak = 20;                  // limit the broken number of bonds in a single iteration, should be an even number

    std::vector<Particle<n_layer> *> top, bottom, left, right, front, back, internal;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        // assign boundary and internal particles
        double delta = 4 * radius / sqrt(2);
        if (p1->xyz[2] > box[5] - 0.6 * delta)
        {
            p1->type = 1;
            top.push_back(p1);
        }
        if (p1->xyz[2] < box[4] + 0.62 * delta)
        {
            p1->type = 2;
            bottom.push_back(p1);
        }
        if (p1->xyz[1] > box[3] - 0.6 * delta)
        {
            p1->type = 3;
            back.push_back(p1);
        }
        // if (p1->xyz[1] < 2 * radius)
        //     front.push_back(p1);
        // if (p1->xyz[0] > 1 - 2 * radius)
        //     right.push_back(p1);
        if (p1->xyz[0] < box[0] + 0.7 * delta)
        {
            p1->type = 4;
            left.push_back(p1);
        }
        if (p1->nb == cell.nneighbors)
        {
            internal.push_back(p1); // particles with full neighbor list
        }

        // assign material properties
        for (int i = 0; i < n_layer; ++i)
        {
            for (auto bd : p1->bond_layers[i])
            {
                // cast to elastic bond (or other type of bonds)
                ElasticBond<n_layer> *elbd = dynamic_cast<ElasticBond<n_layer> *>(bd);
                elbd->setBondProperty(C11, C12, C44, critical_bstrain, nbreak);
                //elbd->setBondProperty(E0, mu0, critical_bstrain, nbreak);
            }
        }
    }

    // simulation settings
    int n_steps = 1;       // number of loading steps
    //double U_stepz{-1e-3}; // step size
    double F_stepz{-100};
    std::vector<LoadStep<n_layer>> load; // load settings for multiple steps
    for (int i = 0; i < n_steps; i++)
    {
        LoadStep<n_layer> step{1}; // 1 means tension loading, -1 means compression loading

        // boundary conditions
        // left: U_stepx; right: -U_stepx
        // front: U_stepy; back: -U_stepy
        // top: U_stepz; bottom: -U_stepz

        step.dispBCs.push_back(DispBC<n_layer>(top, 'z', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(back, 'y', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(left, 'x', 0.0));
        //step.dispBCs.push_back(DispBC<n_layer>(bottom, 'z', U_stepz));
        step.forceBCs.push_back(ForceBC<n_layer>(bottom, 0.0, 0.0, F_stepz));
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
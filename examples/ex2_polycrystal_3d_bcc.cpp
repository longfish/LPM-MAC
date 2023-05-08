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
    UnitCell cell(LatticeType::BCC3D_1Slip, radius);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{0, 10, 0, 10, 0, 10};

    // std::vector<std::array<double, NDIM>> fcc_xyz = createCuboidFCC3D(box, cell, R_matrix);
    // Assembly<n_layer> pt_ass{fcc_xyz, box, cell, BondType::Elastic}; // elastic bond with brittle damage law
    Assembly<n_layer> pt_ass{"../benchmark/BCC_Hailong/BCC_5grains.dump", cell, ParticleType::Elastic}; // read coordinate from dump file

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    double C11{230e3}, C12{135e3}, C44{117e3}; // Elastic constants
    std::vector<Particle<n_layer> *> top, bottom, left, right, front, back, internal;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        // assign boundary and internal particles
        if (p1->xyz[2] > box[5] - 2 * radius)
        {
            // p1->type = 11;
            top.push_back(p1);
        }
        if (p1->xyz[2] < box[4] + 2 * radius)
        {
            // p1->type = 12;
            bottom.push_back(p1);
        }
        if (p1->xyz[1] > box[3] - 2 * radius)
        {
            // p1->type = 13;
            back.push_back(p1);
        }
        if (p1->xyz[1] < box[2] + 2 * radius)
        {
            // p1->type = 14;
            front.push_back(p1);
        }
        if (p1->xyz[0] > box[1] - 2 * radius)
        {
            // p1->type = 15;
            right.push_back(p1);
        }
        if (p1->xyz[0] < box[0] + 2 * radius)
        {
            // p1->type = 16;
            left.push_back(p1);
        }
        if (p1->nb == cell.nneighbors)
        {
            internal.push_back(p1); // particles with full neighbor list
        }

        // assign material properties - need to cast to elastic particle
        ParticleElastic<n_layer> *elpt = dynamic_cast<ParticleElastic<n_layer> *>(p1);
        elpt->setParticleProperty(C11, C12, C44);
    }

    // simulation settings
    int n_steps = 1;                                    // number of loading steps
    double U_stepx{1e-2}, U_stepy{1e-2}, U_stepz{1e-2}; // step size
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

        // another boundary condition
        // step.dispBCs.push_back(DispBC<n_layer>(top, 'z', 0));
        // step.dispBCs.push_back(DispBC<n_layer>(top, 'x', 0));
        // step.dispBCs.push_back(DispBC<n_layer>(top, 'y', 0));
        // step.dispBCs.push_back(DispBC<n_layer>(bottom, 'z', -U_stepz));
        // step.forceBCs.push_back(ForceBC<n_layer>(bottom, 0, 0, -200));

        load.push_back(step);
    }

    pt_ass.updateGeometry();
    pt_ass.updateForceState();

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    SolverStatic<n_layer> solv{pt_ass, StiffnessMode::Analytical, SolverMode::CG, "result_position.dump"}; // stiffness mode and solution mode
    solv.solveProblem(load);

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
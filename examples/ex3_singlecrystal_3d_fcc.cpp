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

    const int n_layer = 2;             // number of neighbor layers (currently only support 2 layers of neighbors)
    double radius = 0.244948964611159; // 1.15; // particle radius
    UnitCell cell(LatticeType::FCC3D, radius);

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{0, 10, 0, 10, 0, 10};

    std::vector<std::array<double, NDIM>> fcc_xyz = createCuboidFCC3D(box, cell, R_matrix);
    Assembly<n_layer> pt_ass{fcc_xyz, box, cell, ParticleType::Elastic}; // elastic bond with brittle damage law
    // Assembly<n_layer> pt_ass{"../geometry/FCC_Hailong/FCC_single.dump", cell, BondType::Elastic}; // read coordinate from dump file

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    double C11{108e3}, C12{61e3}, C44{29e3}; // Elastic constants
    double critical_bstrain = 1;             // critical bond strain

    // simulation settings
    int n_steps = 1;                               // number of loading steps
    double cutoff_ratio = 1.5;                     // nonlocal cutoff ratio
    double step_size = -100;                       // step size for force loading
    double nonlocal_L = 0;                         // nonlocal length scale
    int undamaged_pt_type = -1;                    // particles that dont update damage
    int max_iter = 30, start_index = 0;            // maximum Newton iteration number
    double tol_iter = 1e-5;                        // newton iteration tolerance
    std::string dumpFile{"Crystal_position.dump"}; // output file name

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

        // assign material properties - need to cast to elastic particle
        ParticleElastic<n_layer> *elpt = dynamic_cast<ParticleElastic<n_layer> *>(p1);
        elpt->setParticleProperty(nonlocal_L, C11, C12, C44, critical_bstrain);
    }

    std::vector<LoadStep<n_layer>> load; // load settings for multiple steps
    for (int i = 0; i < n_steps; i++)
    {
        LoadStep<n_layer> step;

        // boundary conditions
        // left: U_stepx; right: -U_stepx
        // front: U_stepy; back: -U_stepy
        // top: U_stepz; bottom: -U_stepz

        step.dispBCs.push_back(DispBC<n_layer>(top, LoadMode::Relative, 'z', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(back, LoadMode::Relative, 'y', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(left, LoadMode::Relative, 'x', 0.0));
        // step.dispBCs.push_back(DispBC<n_layer>(bottom, 'z', U_stepz));
        step.forceBCs.push_back(ForceBC<n_layer>(bottom, LoadMode::Relative, 0.0, 0.0, step_size));
        load.push_back(step);
    }

    pt_ass.searchNonlocalNeighbors(cutoff_ratio);
    pt_ass.updateGeometry();
    pt_ass.updateForceState();

    SolverStatic<n_layer> solv{undamaged_pt_type, pt_ass, StiffnessMode::Analytical, SolverMode::PARDISO, dumpFile, max_iter, tol_iter}; // stiffness mode and solution mode

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    solv.solveProblem(load, start_index);

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
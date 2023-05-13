#include "lpm.h"
#include "particle.h"
#include "assembly.h"
#include "utilities.h"
#include "stiffness.h"
#include "load_step.h"
#include "solver_fatigue.h"

const int n_layer = 2; // number of neighbor layers (currently only support 2 layers of neighbors)

void run()
{
    double start = omp_get_wtime(); // record the CPU time, begin

    double radius = 0.2; // particle radius
    UnitCell cell(LatticeType::Hexagon2D, radius);

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{0.0, 40.0, 0.0, 40.0, 0.0, 8.0};                                                             // thickness is used for force calculation
    Assembly<n_layer> pt_ass{"../geometry/geo1_CT_2DHEX.dump", "../geometry/geo1_CT_2DHEX.bond", cell, ParticleType::FatigueHCF}; // read coordinate from local files

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    double E0 = 205e3, mu0 = 0.29;                           // Young's modulus (MPa) and Poisson's ratio
    double sigma_TS = 420;                                   // tensile strength (MPa)
    double f_A = 9.4e-4, f_b = 0.45, f_c = 6.5, f_eta = 0.2; // fatigue parameters
    double critical_bstrain = 1.0e-3;                        // critical bond strain value at which bond will break

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

        // assign material properties, cast to fatigue damage particle type
        ParticleFatigueHCF<n_layer> *ftpt = dynamic_cast<ParticleFatigueHCF<n_layer> *>(p1);
        ftpt->setParticleProperty(E0, mu0, sigma_TS, f_A, f_b, f_c, f_eta);
    }

    // define one cycle of loading
    double f_max = 600, R = 0.1, f_min = R * f_max;         // cyclic force loading (N)
    std::vector<LoadStep<n_layer>> load_cycle;              // load settings
    std::vector<std::vector<LoadStep<n_layer>>> all_cycles; // all loads

    // cycle step-1
    LoadStep<n_layer> step;
    step.dispBCs.push_back(DispBC<n_layer>(mid_group, 'x', 0.0)); // boundary conditions
    step.dispBCs.push_back(DispBC<n_layer>(mid_group, 'y', 0.0));
    step.forceBCs.push_back(ForceBC<n_layer>(top_group, 0.0, f_min, 0.0));
    step.forceBCs.push_back(ForceBC<n_layer>(bottom_group, 0.0, -f_min, 0.0));
    load_cycle.push_back(step);

    // cycle step-2
    step.forceBCs[0].fy = f_max;
    step.forceBCs[1].fy = -f_max;
    load_cycle.push_back(step);

    // initialize 1000000 cycles
    for (int i = 0; i < 1000000; ++i)
        all_cycles.push_back(load_cycle);

    pt_ass.updateForceState();

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    std::string dumpFile{"CT_2DHex_fatigue.dump"};
    pt_ass.writeDump(dumpFile, 0);

    SolverFatigue<n_layer> solv{pt_ass, StiffnessMode::Analytical, SolverMode::CG, dumpFile, 30, 1e-5}; // stiffness mode and solution mode
    solv.solveProblem(all_cycles);                                                                      // initial loading

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
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
    double E0 = 205e3, mu0 = 0.29;    // Young's modulus (MPa) and Poisson's ratio
    double critical_bstrain = 1.0e-3; // critical bond strain value at which bond will break

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
                ftpt->setParticleProperty(E0, mu0, critical_bstrain);
        }
    }

    // define one cycle of loading
    double f_max = 600, R = 0.1; // step size for cyclic force loading (N)
    double f_mid = 0.5 * (f_max + R * f_max), f_amp = 0.5 * (1 - R) * f_max;
    std::vector<LoadStep<n_layer>> load_cycle, load_initial; // load settings
    LoadStep<n_layer> step{100};                             // initialize with number of cycle jumps after the loading
    for (int i = 0; i < 9; ++i)
    {
        step.dispBCs.push_back(DispBC<n_layer>(mid_group, 'x', 0.0)); // boundary conditions
        step.dispBCs.push_back(DispBC<n_layer>(mid_group, 'y', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_group, 'y', 1e-3));
        step.dispBCs.push_back(DispBC<n_layer>(bottom_group, 'y', -1e-3));
        // step.forceBCs.push_back(ForceBC<n_layer>(top_group, 0.0, f_mid, 0.0));
        // step.forceBCs.push_back(ForceBC<n_layer>(bottom_group, 0.0, -f_mid, 0.0));
        load_cycle.push_back(step);
    }
    load_initial.push_back(step);
    // load_cycle[0].forceBCs[0].fy += f_amp, load_cycle[0].forceBCs[1].fy -= f_amp;
    // load_cycle[1].forceBCs[0].fy += (2 * f_amp), load_cycle[1].forceBCs[1].fy -= (2 * f_amp);
    // load_cycle[2].forceBCs[0].fy += f_amp, load_cycle[2].forceBCs[1].fy -= f_amp;
    // // load_cycle[3]
    // load_cycle[4].forceBCs[0].fy -= f_amp, load_cycle[4].forceBCs[1].fy += f_amp;
    // load_cycle[5].forceBCs[0].fy -= (2 * f_amp), load_cycle[5].forceBCs[1].fy += (2 * f_amp);
    // load_cycle[6].forceBCs[0].fy -= f_amp, load_cycle[6].forceBCs[1].fy += f_amp;
    // // load_cycle[7]

    pt_ass.updateForceState();

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    std::string dumpFile{"CT_2DHex_fatigue.dump"};
    pt_ass.writeDump(dumpFile, 0);

    Solver<n_layer> solv{pt_ass, StiffnessMode::Analytical, SolverMode::CG, dumpFile, cell.nneighbors}; // stiffness mode and solution mode
    solv.solveProblem(pt_ass, load_initial);                                                            // initial loading

    // int n_cycles = 5;
    // for (int i = 0; i < n_cycles; ++i)
    //     one_cycle(solv, pt_ass, load_cycle);

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
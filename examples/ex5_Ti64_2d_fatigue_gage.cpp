#include "lpm.h"
#include "particle.h"
#include "assembly.h"
#include "utilities.h"
#include "stiffness.h"
#include "load_step.h"
#include "solver_static.h"
#include "solver_fatigue.h"
#include "load_step_fatigue.h"

void run()
{
    double start = omp_get_wtime(); // record the CPU time, begin

    const int n_layer = 2; // number of neighbor layers (currently only support 2 layers of neighbors)
    double radius = 0.045; // particle radius
    UnitCell cell(LatticeType::Hexagon2D, radius);

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    // import a particle system
    std::array<double, 2 * NDIM> box{-0.0, 1.7, 0.0, 7 + 1, 0.0, 1};
    std::vector<std::array<double, NDIM>> hex_xyz = createPlateHEX2D(box, cell, R_matrix);
    Assembly<n_layer> pt_ass{hex_xyz, box, cell, ParticleType::FatigueHCF}; // fatigue bond

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    bool is_plane_stress = true;
    double E0 = 110e3, mu0 = 0.34; // AM Ti64 alloy, Young's modulus (MPa) and Poisson's ratio
    // double f_A = 8.933e-5, f_B = 2.019, f_d = 0.5; // fatigue parameters
    double f_A = 5.206e-5, f_B = 2.019, f_d = 0.5; // fatigue parameters
    // double f_A = 5.206e-5, f_B = 2.269, f_d = 0.25; // fatigue parameters
    double f_damage_threshold = 1, f_fatigue_limit_ratio = 1.108;
    // f_A = f_A / f_k * (1 - exp(-f_k));

    // fatigue loading parameters
    double f_max = 200 * 1.7,
           R = 0.1,
           f_min = R * f_max,
           f_range = f_max - f_min;                                 // loading displacement definition
    double cutoff_ratio = 1.5;                                      // nonlocal cutoff ratio
    double nonlocal_L = 4;                                          // nonlocal length scale
    double tau = 0.0004;                                            // fatigue time mapping parameter
    int max_iter = 30;                                              // maximum iteration number of Newton-Raphson algorithm
    int start_index = 0;                                            // start index of solution procedure
    double tol_iter = 1e-5;                                         // tolerance of the NR iterations
    int undamaged_pt_type = 4;                                      // undamaged particle type
    std::string dumpFile{"DogBone_2d_fatigue.dump"};                // output file name
    std::string loadFile{"../loading/force_ca_R0.1_Ti64_S90.txt"}; // loading file name

    std::vector<Particle<n_layer> *> top_group, left_group, bottom_group;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        // assign boundary and internal particles
        if (p1->xyz[0] < radius)
        {
            left_group.push_back(p1); // left
            p1->type = 1;
        }
        if (p1->xyz[1] > 8. / 2 + 7. / 2)
        {
            top_group.push_back(p1); // top
            p1->type = 2;
        }
        if (p1->xyz[1] < 8. / 2 - 7. / 2)
        {
            bottom_group.push_back(p1); // bottom
            p1->type = 3;
        }
        if (p1->type == 2 || p1->type == 3)
            p1->type = undamaged_pt_type; // undamaged part

        // assign material properties, cast to fatigue damage particle type
        ParticleFatigueHCF<n_layer> *ftpt = dynamic_cast<ParticleFatigueHCF<n_layer> *>(p1);
        ftpt->setParticleProperty(nonlocal_L, is_plane_stress, E0, mu0, f_A, f_B, f_d, f_damage_threshold, f_fatigue_limit_ratio);
    }

    pt_ass.searchNonlocalNeighbors(cutoff_ratio);
    pt_ass.updateGeometry();
    pt_ass.updateForceState();

    SolverFatigue<n_layer> solv_fatigue{undamaged_pt_type,
                                        pt_ass,
                                        StiffnessMode::Analytical,
                                        SolverMode::CG,
                                        TimeMapMode::Linear,
                                        tau, dumpFile, max_iter, tol_iter}; // stiffness mode and solution mode

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    /*********************************************************************************************/
    // perform static loading, f_min
    std::vector<LoadStep<n_layer>> load_static; // static load settings

    // boundary conditions
    int n_incre_static = 1;
    LoadStep<n_layer> step0;
    step0.dispBCs.push_back(DispBC<n_layer>(left_group, LoadMode::Relative, 'x', 0.0));
    step0.forceBCs.push_back(ForceBC<n_layer>(top_group, LoadMode::Relative, 0.0, f_min / n_incre_static, 0.0));
    // step0.dispBCs.push_back(DispBC<n_layer>(top_group, LoadMode::Relative, 'y', f_min / n_incre_static));
    step0.dispBCs.push_back(DispBC<n_layer>(bottom_group, LoadMode::Relative, 'y', 0.0));
    for (int i = 0; i < n_incre_static; ++i)
        load_static.push_back(step0);

    solv_fatigue.solveProblemStatic(load_static, start_index);
    // std::cout << solv_fatigue.ass.pt_sys[23520]->state_var[2] << std::endl;

    /*********************************************************************************************/
    // perform cyclic loading (incrementally)
    int n_incre_fatigue = 1;
    solv_fatigue.readLoad(loadFile);
    solv_fatigue.solveProblemCyclic(FatigueLoadType::LoadUniaxialForce, n_incre_fatigue, {left_group, top_group, bottom_group});

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
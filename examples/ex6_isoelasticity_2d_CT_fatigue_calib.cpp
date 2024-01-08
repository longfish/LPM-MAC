#include "lpm.h"
#include "particle.h"
#include "assembly.h"
#include "utilities.h"
#include "stiffness.h"
#include "load_step.h"
#include "solver_fatigue.h"

void run()
{
    double start = omp_get_wtime(); // record the CPU time, begin

    const int n_layer = 2; // number of neighbor layers (currently only support 2 layers of neighbors)
    double radius = 0.2;   // particle radius
    UnitCell cell(LatticeType::Hexagon2D, radius);

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    Eigen::MatrixXd R_matrix = createRMatrix(eulerflag, angles);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{0.0, 63.5, 0.0, 60.96, 0.0, 5.0};                                                                      // thickness is used for force calculation
    Assembly<n_layer> pt_ass{"../geometry/geo1_CT_2DHEX_Zhao.dump", "../geometry/geo1_CT_2DHEX_Zhao.bond", cell, ParticleType::FatigueHCF}; // read coordinate from local files

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    bool is_plane_stress = true;
    double E0 = 71.7e3, mu0 = 0.306;                // Al 7075T6, Young's modulus (MPa) and Poisson's ratio
    double f_A = 1.94e-4, f_B = 0.581, f_d = 1.079; // fatigue parameters
    // double f_A = 2.003e-4, f_B = 0.581, f_d = 1.079; // fatigue parameters
    double f_damage_threshold = 1, f_fatigue_limit_ratio = 1.108;

    // fatigue loading parameters
    double force_max = 2500 / (box[5] - box[4]),
           R = 0.1,
           force_min = R * force_max,
           force_range = force_max - force_min;             // loading force definition
    double cutoff_ratio = 1.5;                              // nonlocal cutoff ratio
    double tau = 0.001;                                     // fatigue time mapping parameter
    int max_iter = 30;                                      // maximum iteration number of Newton-Raphson algorithm                                                                                                 /* maximum Newton iteration number */
    int start_index = 0;                                    // start index of solution procedure
    double tol_iter = 1e-5;                                 // tolerance of the NR iterations
    double nonlocal_L = 0.8;                                // nonlocal length scale
    int undamaged_pt_type = 4;                              // undamaged particle type
    std::string dumpFile{"CTZhao_fatigue_calib_L0.8.dump"}; // output file name

    std::vector<Particle<n_layer> *> top_group, bottom_group, mid_group;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        if (p1->id == 12897 || p1->id == 12785)
        {
            mid_group.push_back(p1); // mid, to fix in y direction
            p1->type = 1;
        }
        if (p1->id >= 21217 && p1->id <= 21218)
        {
            top_group.push_back(p1); // top
            p1->type = 2;
        }
        if (p1->id >= 4564 && p1->id <= 4565)
        {
            bottom_group.push_back(p1); // bottom
            p1->type = 3;
        }
        // assign boundary and internal particles
        if (p1->xyz[1] >= 60.96 - 16.46 - 12.7 / 2)
            p1->type = undamaged_pt_type;
        if (p1->xyz[1] <= 16.46 + 12.7 / 2)
            p1->type = undamaged_pt_type; // particles that not update damage
        // if (p1->xyz[0] >= 50.8 - 12.7 / 2 && p1->xyz[0] <= 50.8 + 12.7 / 2 && p1->xyz[1] >= 60.96 - 16.46)
        //     p1->type = undamaged_pt_type;
        // if (p1->xyz[0] >= 50.8 - 12.7 / 2 && p1->xyz[0] <= 50.8 + 12.7 / 2 && p1->xyz[1] <= 16.46)
        //     p1->type = undamaged_pt_type; // particles that not update damage

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

    int n_incre_static = 1;
    LoadStep<n_layer> step0;
    step0.dispBCs.push_back(DispBC<n_layer>(mid_group, LoadMode::Relative, 'y', 0));
    step0.dispBCs.push_back(DispBC<n_layer>(top_group, LoadMode::Relative, 'x', 0.0));
    step0.dispBCs.push_back(DispBC<n_layer>(bottom_group, LoadMode::Relative, 'x', 0.0));
    step0.forceBCs.push_back(ForceBC<n_layer>(top_group, LoadMode::Relative, 0.0, force_min, 0.0));
    step0.forceBCs.push_back(ForceBC<n_layer>(bottom_group, LoadMode::Relative, 0.0, -force_min, 0.0));
    for (int i = 0; i < n_incre_static; ++i)
        load_static.push_back(step0);

    solv_fatigue.solveProblemStatic(load_static, start_index);

    /*********************************************************************************************/
    // perform cyclic loading
    int n_incre_fatigue = 1;
    solv_fatigue.readLoad("../loading/force_ca_R0.1_calib.txt");
    solv_fatigue.solveProblemCyclic(FatigueLoadType::LoadCTForce, n_incre_fatigue, {mid_group, top_group, bottom_group});

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
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
    double radius = 0.2;   // particle radius
    UnitCell cell(LatticeType::SimpleCubic3D, radius);

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    // import a particle system
    Assembly<n_layer> pt_ass{"../geometry/geo2_DogBone_3DSC.dump", "../geometry/geo2_DogBone_3DSC.bond", cell, ParticleType::FatigueHCF}; // read coordinate from local files

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    double E0 = 71.7e3, mu0 = 0.306;              // Al 7075-T651, Young's modulus (MPa) and Poisson's ratio
    double f_A = 2.222e-5, f_B = 2.401, f_k = 10; // fatigue parameters
    double f_damage_threshold = 1, f_fatigue_limit_ratio = 1.108;
    f_A = f_A / f_k * (1 - exp(-f_k));

    // fatigue loading parameters
    double d_max = 0.10757,
           R = 0.0,
           d_min = R * d_max,
           d_range = d_max - d_min;                    // loading displacement definition
    double tau = 0.0001;                               // fatigue time mapping parameter
    int max_iter = 30;                                 // maximum iteration number of Newton-Raphson algorithm                                                                                                 /* maximum Newton iteration number */
    int start_index = 0;                               // start index of solution procedure
    double tol_iter = 1e-5;                            // tolerance of the NR iterations
    int undamaged_pt_type = 3;                         // undamaged particle type
    std::string dumpFile{"DogBone_3DSC_fatigue.dump"}; // output file name

    std::vector<Particle<n_layer> *> top_group, bottom_group;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        // assign boundary and internal particles
        if (p1->xyz[1] < 19.05 - 1.5 * radius || p1->xyz[1] > 19.05 + 1.5 * radius)
            p1->type = undamaged_pt_type; // undamaged part
        if (p1->xyz[1] > 37.1)
        {
            top_group.push_back(p1); // top
            p1->type = 1;
        }
        if (p1->xyz[1] < 1)
        {
            bottom_group.push_back(p1); // bottom
            p1->type = 2;
        }

        // assign material properties, cast to fatigue damage particle type
        ParticleFatigueHCF<n_layer> *ftpt = dynamic_cast<ParticleFatigueHCF<n_layer> *>(p1);
        ftpt->setParticleProperty(E0, mu0, f_A, f_B, f_k, f_damage_threshold, f_fatigue_limit_ratio);
    }

    pt_ass.updateGeometry();
    pt_ass.updateForceState();

    SolverFatigue<n_layer> solv_fatigue{undamaged_pt_type, pt_ass,
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
    step0.dispBCs.push_back(DispBC<n_layer>(top_group, LoadMode::Relative, 'x', 0.0));
    step0.dispBCs.push_back(DispBC<n_layer>(top_group, LoadMode::Relative, 'y', d_min / n_incre_static));
    step0.dispBCs.push_back(DispBC<n_layer>(top_group, LoadMode::Relative, 'z', 0.0));
    step0.dispBCs.push_back(DispBC<n_layer>(bottom_group, LoadMode::Relative, 'x', 0.0));
    step0.dispBCs.push_back(DispBC<n_layer>(bottom_group, LoadMode::Relative, 'y', 0.0));
    step0.dispBCs.push_back(DispBC<n_layer>(bottom_group, LoadMode::Relative, 'z', 0.0));
    for (int i = 0; i < n_incre_static; ++i)
        load_static.push_back(step0);

    solv_fatigue.solveProblemStatic(load_static, start_index);
    // std::cout << solv_fatigue.ass.pt_sys[23520]->state_var[2] << std::endl;

    /*********************************************************************************************/
    // perform cyclic loading (incrementally)
    int n_incre_fatigue = 1;
    solv_fatigue.readLoad("../loading/disp_ca_R0.txt");
    solv_fatigue.solveProblemCyclic(FatigueLoadType::LoadUniaxialDisp, n_incre_fatigue, top_group, bottom_group);

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
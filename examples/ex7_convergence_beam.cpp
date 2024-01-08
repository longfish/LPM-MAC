#include "lpm.h"
#include "particle.h"
#include "assembly.h"
#include "utilities.h"
#include "stiffness.h"
#include "load_step.h"
#include "solver_static.h"
#include "solution.h"

template <typename T>
void writeMatrix(const char *dataName, T *data, int l)
{
    FILE *fpt;
    fpt = fopen(dataName, "w+");
    for (int i = 0; i < l; i++)
    {
        if constexpr (std::is_integral_v<T>)
        { // constexpr only necessary on first statement
            fprintf(fpt, " %lld\n", data[i]);
        }
        else if (std::is_floating_point_v<T>)
        { // automatically constexpr
            fprintf(fpt, " %.5e\n", data[i]);
        }
    }

    fclose(fpt);
}

std::vector<double> calculateBoundary(std::vector<std::array<double, NDIM>> &xyz)
{
    double left = INT_MAX, right = INT_MIN;
    for (auto &p : xyz)
    {
        if (p[0] < left)
            left = p[0];
        if (p[0] > right)
            right = p[0];
    }
    return {left, right};
}

void run()
{
    double start = omp_get_wtime(); // record the CPU time, begin

    const int n_layer = 2; // number of neighbor layers (currently only support 2 layers of neighbors)
    double radius = 0.12;  // particle radius
    UnitCell cell(LatticeType::Hexagon2D, radius);

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    Eigen::MatrixXd R_matrix = createRMatrix(eulerflag, angles);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{0.0, 30, -5, 5, 0.0, 1.0}; // thickness is used for force calculation
    std::vector<std::array<double, NDIM>> hex_xyz = createPlateHEX2D(box, cell, R_matrix);
    Assembly<n_layer> pt_ass{hex_xyz, box, cell, ParticleType::Elastic}; // elastic bond

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // determine the particle boundaries
    std::vector<double> bounds = calculateBoundary(hex_xyz);
    printf("left boundary: %f; right boundary: %f\n\n", bounds[0], bounds[1]);

    // material elastic parameters setting, MPa
    bool is_plane_stress = false;
    double E0 = 69e3, mu0 = 0.33;  // Young's modulus (MPa) and Poisson's ratio
    double f_critical_bstrain = 1; // critical bond strain

    // fatigue loading parameters
    double force_value = 2;                        // loading force definition, Newton
    int n_steps = 1;                               // number of loading steps
    double cutoff_ratio = 1.5;                     // nonlocal cutoff ratio
    double tau = 0.001;                            // fatigue time mapping parameter
    int max_iter = 30;                             // maximum iteration number of Newton-Raphson algorithm                                                                                                 /* maximum Newton iteration number */
    int start_index = 0;                           // start index of solution procedure
    double tol_iter = 1e-5;                        // tolerance of the NR iterations
    double nonlocal_L = EPS;                       // nonlocal length scale
    int undamaged_pt_type = 4;                     // undamaged particle type
    std::string dumpFile{"beam_convergence.dump"}; // output file name

    std::vector<Particle<n_layer> *> left_group, right_group;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        if (p1->xyz[0] < bounds[0] + EPS)
        {
            left_group.push_back(p1); // left
            p1->type = 2;
        }
        if (p1->xyz[0] > bounds[1] - EPS)
        {
            right_group.push_back(p1); // right
            p1->type = 1;
        }
        // if (p1->type == 2 || p1->type == 3)
        //     p1->type = undamaged_pt_type; // undamaged part

        // assign material properties, cast to elastic particle type
        ParticleElastic<n_layer> *ftpt = dynamic_cast<ParticleElastic<n_layer> *>(p1);
        ftpt->setParticleProperty(nonlocal_L, is_plane_stress, E0, mu0, f_critical_bstrain);
    }

    std::vector<LoadStep<n_layer>> load; // load settings for multiple steps
    for (int i = 0; i < n_steps; i++)
    {
        LoadStep<n_layer> step;

        // boundary conditions
        step.dispBCs.push_back(DispBC<n_layer>(right_group, LoadMode::Relative, 'x', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(right_group, LoadMode::Relative, 'y', 0.0));
        step.forceBCs.push_back(ForceBC<n_layer>(left_group, LoadMode::Relative, 0.0, force_value, 0.0));

        load.push_back(step);
    }

    pt_ass.searchNonlocalNeighbors(cutoff_ratio);
    pt_ass.updateGeometry();
    pt_ass.updateForceState();
    pt_ass.updateStateVar();

    SolverStatic<n_layer> solv{undamaged_pt_type, pt_ass, StiffnessMode::Analytical, SolverMode::CG, dumpFile, max_iter, tol_iter}; // stiffness mode and solution mode

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    solv.solveProblem(load, start_index);

    // test the convergence behavior
    Solution<n_layer> sol_conv{pt_ass, ExactSolution::Beam001};
    sol_conv.computeSolBeam001(box[1] - box[0], box[3] - box[2], force_value, E0, mu0);
    printf("L2 error is: %f\n", sol_conv.computeL2Error());
    writeMatrix("exact_sol.txt", sol_conv.exact_u.data(), sol_conv.exact_u.size());
    writeMatrix("approx_sol.txt", sol_conv.approx_u.data(), sol_conv.approx_u.size());
    writeMatrix("error_sol.txt", sol_conv.error.data(), sol_conv.error.size());

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
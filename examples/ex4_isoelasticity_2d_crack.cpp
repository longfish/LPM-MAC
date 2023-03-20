#include "lpm.h"
#include "particle.h"
#include "assembly.h"
#include "utilities.h"
#include "stiffness.h"
#include "load_step.h"
#include "solver.h"

void writeK_global(const char *dataName,double *K_global, int l)
{
    FILE *fpt;
    fpt = fopen(dataName, "w+");
    for (int i = 0; i < l; i++)
    {
        fprintf(fpt, " %8.3f\n", K_global[i]);
    }

    fclose(fpt);
}

void run()
{
    double start = omp_get_wtime(); // record the CPU time, begin

    const int n_layer = 2; // number of neighbor layers (currently only support 2 layers of neighbors)
    double radius = 0.2;   // particle radius
    UnitCell cell(LatticeType::Square2D, radius);

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{0.0, 40.0, 0.0, 40.0, 0.0, 8.0};
    Assembly<n_layer> pt_ass{"../geometry/geo1_CT_2DSquare.dump", "../geometry/geo1_CT_2DSquare.bond", cell, BondType::Elastic}; // read coordinate from local files
    // std::vector<std::array<double, NDIM>> sq_xyz = createPlateSQ2D(box, cell, R_matrix);
    // Assembly<n_layer> pt_ass{sq_xyz, box, cell, BondType::Elastic}; // elastic bond with brittle damage law

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    double E0 = 69e3, mu0 = 0.3;      // Young's modulus and Poisson's ratio
    double critical_bstrain = 1.0e-2; // critical bond strain value at which bond will break
    int nbreak = 20;                  // limit the broken number of bonds in a single iteration, should be an even number

    std::vector<Particle<n_layer> *> top_group, bottom_group, internal_group;
    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        // assign boundary and internal particles
        if (p1->xyz[1] > box[3] - 2 * radius)
        {
            top_group.push_back(p1); // top
            p1->type = 1;
        }
        if (p1->xyz[1] < box[2] + 2 * radius)
        {
            bottom_group.push_back(p1); // bottom
            p1->type = 2;
        }
        if (p1->nb == cell.nneighbors)
            internal_group.push_back(p1); // particles with full neighbor list

        // assign material properties
        for (int i = 0; i < n_layer; ++i)
        {
            for (auto bd : p1->bond_layers[i])
            {
                // cast to elastic bond (or other type of bonds)
                BondElastic<n_layer> *elbd = dynamic_cast<BondElastic<n_layer> *>(bd);
                elbd->setBondProperty(E0, mu0, critical_bstrain, nbreak);
            }
        }
    }

    // simulation settings
    int n_steps = 1;          // number of loading steps
    double step_size = -1e-3; // step size for force or displacement loading

    std::vector<LoadStep<n_layer>> load; // load settings for multiple steps
    for (int i = 0; i < n_steps; i++)
    {
        LoadStep<n_layer> step{1}; // 1 means tension loading, -1 means compression loading

        // boundary conditions
        step.dispBCs.push_back(DispBC<n_layer>(top_group, 'x', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(top_group, 'y', -step_size));
        //step.dispBCs.push_back(DispBC<n_layer>(top_group, 'z', 0.0));
        step.dispBCs.push_back(DispBC<n_layer>(bottom_group, 'y', step_size));
        load.push_back(step);
    }

    pt_ass.updateForceState();

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    Solver<n_layer> solv{pt_ass, StiffnessMode::Analytical, SolverMode::CG, "CT_2DSquare_position.dump"}; // stiffness mode and solution mode
    solv.solveProblem(pt_ass, load);

    // char KFile[] = "result_K.txt";
    // writeK_global(KFile, solv.stiffness.K_global, 12600);

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);
}
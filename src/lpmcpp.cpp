/******************************************************************************
 * Implemented the solving process of nonlocal lattice particle method(LPM) in
 * general cases(2 layers of neighbors)
 * This code can deal with general elastic and plasticity problems under careful
 * modifications
 *
 * Supported lattice types:
 *     2D square, 2D hexagon, 3D simple cubic, 3D BCC & FCC crystal(crystal
 *     plasticity simulation)
 *
 * Contact information:
 *     Changyu Meng, PhD candidate, Arizona State University, cmeng12@asu.edu
 ******************************************************************************/

#include "lpm.h"
#include "particle.h"
#include "assembly.h"
#include "utilities.h"
#include "stiffness.h"
#include "load_step.h"
#include "solver.h"

/************************************************************************/
/****************************** Main procedure **************************/
/************************************************************************/
int main(int argc, char *argv[])
{
    printf("\n==================================================\n");
    printf("            Nonlocal LPM C++ Program              \n");
    printf("==================================================\n");

    const int nt = omp_get_max_threads(); /* maximum number of threads provided by the computer */
    printf("OpenMP with %d threads\n", nt);

    double start = omp_get_wtime(); // record the CPU time, begin

    const int n_layer = 2; // number of neighbor layers (currently only support 2 layers of neighbors)
    double radius = 0.4;   // 1.15; // particle radius
    UnitCell cell(LatticeType::FCC3D, radius);

    // Euler angles setting for system rotation
    // flag is 0 ~ 2 for different conventions, (0: direct rotation; 1: Kocks convention; 2: Bunge convention)
    // angle1, angle2 and an angle3 are Euler angles in degree, double
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    // create a simulation box
    // xmin; xmax; ymin; ymax; zmin; zmax
    std::array<double, 2 * NDIM> box{-0.2, 10.2, -0.2, 10.2, -0.2, 10.2};

    //std::vector<std::array<double, NDIM>> fcc_xyz = createCuboidFCC3D(box, cell, R_matrix);
    //Assembly<n_layer> pt_ass{fcc_xyz, box, cell, BondType::Elastic}; // elastic bond with brittle damage law
    Assembly<n_layer> pt_ass{"../examples/position.dump", cell, BondType::Elastic}; // read coordinate from dump file

    printf("\nParticle number is %d\n", pt_ass.nparticle);

    // material elastic parameters setting, MPa
    double E0 = 146e3, mu0 = 0.3;     // Young's modulus and Poisson's ratio
    double critical_bstrain = 1.0e-2; // critical bond strain value at which bond will break
    int nbreak = 20;                  // limit the broken number of bonds in a single iteration, should be an even number

    for (Particle<n_layer> *p1 : pt_ass.pt_sys)
    {
        // assign boundary and internal particles
        if (p1->xyz[2] > 10.0 - 1.2 * radius)
            p1->type = 1; // top
        if (p1->xyz[2] < 1.2 * radius)
            p1->type = 2; // bottom
        if (p1->nb == cell.nneighbors)
            p1->type = 3; // particles with full neighbor list

        // assign material properties
        for (int i = 0; i < n_layer; ++i)
        {
            for (auto bd : p1->bond_layers[i])
            {
                // cast to elastic bond (or other type of bonds)
                ElasticBond<n_layer> *elbd = dynamic_cast<ElasticBond<n_layer> *>(bd);
                elbd->setBondProperty(E0, mu0, critical_bstrain, nbreak);
            }
        }
    }

    // simulation settings
    int n_steps = 1;                                 // number of loading steps
    double step_size = -2000.0;                      // step size for force or displacement loading

    std::vector<LoadStep> load; // load settings for multiple steps
    for (int i = 0; i < n_steps; i++)
    {
        LoadStep step{1}; // 1 means tension loading, -1 means compression loading

        // boundary conditions
        step.dispBCs.push_back(DispBC(1, 'x', 0.0));
        step.dispBCs.push_back(DispBC(1, 'y', 0.0));
        step.dispBCs.push_back(DispBC(1, 'z', 0.0));
        step.forceBCs.push_back(ForceBC(2, 0.0, 0.0, step_size));
        load.push_back(step);
    }

    pt_ass.updateForceState();

    double initrun = omp_get_wtime();
    printf("Initialization finished in %f seconds\n\n", initrun - start);

    Solver<n_layer> solv{pt_ass, StiffnessMode::Analytical, SolverMode::CG, "result_position.dump"}; // stiffness mode and solution mode
    solv.solveProblem(pt_ass, load);

    double finish = omp_get_wtime();
    printf("Computation time for total steps: %f seconds\n\n", finish - start);

    // frees unused memory allocated by the Intel MKL Memory Allocator
    mkl_free_buffers();

    return 0;
}

/************************************************************************/
/*************************** End main procedures ************************/
/************************************************************************/

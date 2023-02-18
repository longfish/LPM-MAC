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

//#include "ex1_elasticity_3d.cpp"
#include "ex2_polycrystal_3d_bcc.cpp"
//#include "ex3_singlecrystal_3d_fcc.cpp"

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

    run();

    mkl_free_buffers(); // frees unused memory allocated by the Intel MKL Memory Allocator

    return 0;
}

/************************************************************************/
/*************************** End main procedures ************************/
/************************************************************************/

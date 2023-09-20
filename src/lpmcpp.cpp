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

// #include "geo1_CT_configuration.cpp"
// #include "geo2_dogbone_configuration.cpp"
// #include "geo2_dogbone_config2d.cpp"

// #include "ex1_isoelasticity_2d_hex.cpp"
// #include "ex1_isoelasticity_3d_sc.cpp"
// #include "ex2_polycrystal_3d_bcc.cpp"
// #include "ex3_singlecrystal_3d_fcc.cpp"
// #include "ex4_isoelasticity_2d_crack.cpp"
// #include "ex4_isoelasticity_2d_crack_shear.cpp"
// #include "ex4_isoelasticity_2d_crack_shear_sq.cpp"
// #include "ex4_isoelasticity_2d_crack_shear_bond.cpp"
// #include "ex4_isoelasticity_2d_crack_shear_sq_bond.cpp"
// #include "ex5_isoelasticity_2d_fatigue_gage.cpp"
// #include "ex5_isoelasticity_2d_fatigue.cpp"
// #include "ex5_isoelasticity_3d_fatigue.cpp"
// #include "ex5_Ti64_2d_fatigue_gage.cpp"
// #include "ex6_isoelasticity_2d_CT_fatigue_calib.cpp"
// #include "ex6_isoelasticity_2d_CT_fatigue.cpp"
// #include "ex6_Ti64_2d_fatigue_crack_G6a.cpp"
// #include "ex6_Ti64_2d_fatigue_crack_calib.cpp"
// #include "ex7_convergence_beam.cpp"
#include "plasticity/J2_3DSC.cpp"

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

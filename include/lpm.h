#pragma once
#ifndef LPM3D_H
#define LPM3D_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <omp.h>
#include <mkl.h>
#include <mkl_pardiso.h>
#include <mkl_types.h>
#include <mkl_rci.h>
#include <mkl_blas.h>
#include <mkl_spblas.h>
#include <mkl_service.h>

#include "boundary.h"
#include "data_handler.h"
#include "initialization.h"
#include "neighbor.h"
#include "solver.h"
#include "stiffness.h"
#include "utilities.h"
#include "constitutive.h"

#define PI 3.14159265358979323846
#define NDIM 3        /* max number of dimensions, no matter what dimension of the problem */
#define TOL 1e-3      /* tolerance to determine whether two quantities are equal */
#define MAXLINE 500   /* maximum line number */
#define MAXSMALL 20   /* maximum number used for small iteration etc. */
#define MAXITER 100   /* maximum global iteration number */
#define MAXSLIPSYS 50 /* maximum size of slip systems */
#define TOLITER 1e-4  /* newton iteration tolerance */
#define EPS 1e-6      /* small perturbation coefficient used to compute the stiffness matrix */

/* simple functions */
#define LEN(arr) (sizeof(arr) / sizeof(arr[0])) /* Length of an array */
#define MAX(x, y) ((x) < (y) ? (y) : (x))       /* Maximum of two variables */
#define MIN(x, y) ((x) > (y) ? (y) : (x))       /* Minimum of two variables */
#define SIGN(x) ((x) < (0) ? (-1.0) : (1))      /* Sign of variable, zero is considered as positive */

#define SY(x) (620.0 + 3300.0 * (1.0 - exp(-0.4 * x)))                                       /* isotropic hardening function for J2 plasticity, MPa */
#define DAM_PHI(x) (1.0 / damage_L / sqrt(2 * PI) * exp(-0.5 * x * x / damage_L / damage_L)) /* nonlocal helper function for damage accumulation */

struct DispBCs
{
    int type;
    char flag;
    double step;
};

struct ForceBCs
{
    int type;
    char flag1;
    double step1;
    char flag2;
    double step2;
    char flag3;
    double step3;
};

struct UnitCell{
    int lattice, dim;
    int nneighbors, nneighbors1, nneighbors2;
    int nneighbors_AFEM, nneighbors_AFEM1, nneighbors_AFEM2;
    double radius, neighbor1_cutoff, neighbor2_cutoff, particle_volume; /* volume of unit cell */
};

/* Declaration of global variables */
/* int */
extern int nparticle;
extern int nslip_face;
extern int max_nslip_vector, nslipSys, nbreak;

extern int *IK, *JK, *type, *dispBC_index, *fix_index, *nslip_vector, *lacknblist, *pl_flag;
extern int *nb, *nb_initial, *nb_conn, *state_v;
extern int **neighbors, **neighbors1, **neighbors2, **neighbors_AFEM;
extern int **K_pointer, **conn, **nsign, **cp_Jact;

/* double precision float */
extern double cp_tau0[3], cp_taus[3], cp_eta, cp_p, cp_h0, cp_maxloop;
extern double cp_q, cp_dtime, cp_theta, J2_H, J2_xi, J2_C, damage_L;
extern double damage_threshold, damageb_A, damagec_A, critical_bstrain;

extern double *K_global, *plastic_K_global, *residual, *Pin, *Pex, *Pex_temp, *disp, *sigmay, *cp_dA;
extern double *reaction_force, *damage_visual;
extern double *J2_dlambda, *J2_stresseq, *J2_stressm, *J2_triaxiality;

extern double **xyz, **xyz_initial, **xyz_temp, **distance, **distance_initial, **F, **csx, **csy, **csz;
extern double **dL_total, **TdL_total, **csx_initial, **csy_initial, **csz_initial, **Ce;
extern double **slip_normal, **schmid_tensor, **schmid_tensor_local, **cp_RSS, **stress_tensor;
extern double **cp_Cab, **strain_tensor, **dL_ave, **ddL_total, **TddL_total, **F_temp, **ddLp;
extern double **dL, **ddL, **bond_stress, **damage_broken, **damage_w, **bond_stretch, **bond_vector, **bond_force;

extern double **Kn, **Tv, **J2_alpha, **damage_local, **damage_nonlocal, **cp_A, **cp_dgy, **cp_dA_single, **J2_beta_eq;
extern double ***slip_vector, ***dLp, ***J2_beta, ***damage_D, ***cp_gy, ***cp_A_single;

/* function prototype */
void computeBondForceGeneral(int plmode, int temp, struct UnitCell cell);
void initMatrices(struct UnitCell cell);
int updateDamageGeneral(const char *dataName, int tstep, int plmode, struct UnitCell cell);

#endif

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
#include "unit_cell.h"

// Define the format to printf MKL_INT values
#ifndef MKL_ILP64
#define IFORMAT "%i"
#else
#define IFORMAT "%lli"
#endif

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

/* Declaration of global variables */
/* int */
extern int nparticle;
extern int nbreak;

extern MKL_INT *IK, *JK;
extern int *type, *dispBC_index, *fix_index, *pl_flag;
extern int *nb, *nb_initial, *nb_conn;
extern int **neighbors, **neighbors1, **neighbors2;
extern int **K_pointer, **conn, **nsign;

/* double precision float */
extern double critical_bstrain;

extern double *K_global, *residual, *Pin, *Pex, *Pex_temp, *disp;
extern double *reaction_force, *damage_visual;

extern double **xyz, **xyz_initial, **xyz_temp, **distance, **distance_initial, **F, **csx, **csy, **csz;
extern double **dL_total, **TdL_total, **csx_initial, **csy_initial, **csz_initial, **Ce;
extern double **stress_tensor;
extern double **strain_tensor, **ddL_total, **TddL_total, **F_temp;
extern double **dL, **ddL, **bond_stress, **damage_broken, **damage_w;

extern double **Kn, **Tv;
extern double ***damage_D;

#endif

#pragma once
#ifndef LPM3D_H
#define LPM3D_H

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <iterator>

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
#include <boost/range/adaptor/indexed.hpp>

#include <omp.h>
#include <mkl.h>
#include <mkl_pardiso.h>
#include <mkl_types.h>
#include <mkl_rci.h>
#include <mkl_blas.h>
#include <mkl_spblas.h>
#include <mkl_service.h>

using namespace boost::adaptors;

/* macros */
#define PI 3.14159265358979323846
#define EPS  1e-6
#define NDIM 3                             /* max number of dimensions, no matter what dimension of the problem */
#define MAX(x, y) ((x) < (y) ? (y) : (x))  /* Maximum of two variables */
#define MIN(x, y) ((x) > (y) ? (y) : (x))  /* Minimum of two variables */
#define SIGN(x) ((x) < (0) ? (-1.0) : (1)) /* Sign of variable, zero is considered as positive */

// Define the format to printf MKL_INT values
#ifndef MKL_ILP64
#define IFORMAT "%i"
#else
#define IFORMAT "%lli"
#endif

enum class SolverMode : char
{
    CG,
    PARDISO
};

enum class StiffnessMode : char
{
    Analytical,
    FiniteDifference
};

enum class ParticleType : char
{
    Elastic,
    ElasticDamage,
    FatigueHCF
};

enum class LatticeType : char
{
    Square2D,
    Hexagon2D,
    SimpleCubic3D,
    FCC3D,
    BCC3D_1Slip,
    BCC3D_2Slip,
    BCC3D_3Slip
};

#endif

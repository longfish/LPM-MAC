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
#include <queue>
#include <deque>

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
#define EPS 1e-6
#define NDIM 3 /* max number of dimensions, no matter what dimension of the problem */

/* use for nonlocal damage modeling */
// double func_phi(const double x, const double L)
// {
//     return 1.0 / L / sqrt(2 * PI) * exp(-0.5 * x * x / L / L);
// }
double func_phi(const double x, const double L)
{
    double p = 8, q = 2; // nonlocal constants
    return pow(1 / (1 + pow(x / L, p)), q);
}

// Define the format to printf MKL_INT values
#ifndef MKL_ILP64
#define IFORMAT "%i"
#else
#define IFORMAT "%lli"
#endif

enum class ExactSolution : char
{
    Beam001
};

enum class LoadMode : char
{
    Relative,
    Absolute
};

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

enum class TimeMapMode : char
{
    Linear,
    Exponential
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

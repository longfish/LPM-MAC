#pragma once
#ifndef UTILITIES_H
#define UTILITIES_H

#include "unit_cell.h"
#include "lpm.h"

// basic lpm computations
double * createRMatrix(int eulerflag, double angles[]);
std::vector<std::array<double, NDIM>> createCuboidSC3D(double box[], UnitCell cell, double R_matrix[]);

#endif

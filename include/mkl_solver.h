#pragma once
#ifndef MKL_SOLVER_H
#define MKL_SOLVER_H

#include "lpm.h"
#include "unit_cell.h"

void solverPARDISO(UnitCell cell);
void solverCG(UnitCell cell);

#endif

#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include "lpm.h"
#include "unit_cell.h"

void solverPARDISO(UnitCell cell);
void solverCG(UnitCell cell);

#endif

#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include "lpm.h"

struct UnitCell;

void solverPARDISO(struct UnitCell cell);
void solverCG(struct UnitCell cell);

#endif

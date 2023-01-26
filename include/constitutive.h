#pragma once
#ifndef PLASTICITY_H
#define PLASTICITY_H

#include "lpm.h"
#include "unit_cell.h"

void computeBondForceElastic(int i, UnitCell cell);
void computeBondForceIncrementalUpdating(int ii, UnitCell cell);
void updateCrack(UnitCell cell);

int updateBrittleDamage(const char *dataName, int tstep, int nbreak);

#endif

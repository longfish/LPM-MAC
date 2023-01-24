#pragma once
#ifndef PLASTICITY_H
#define PLASTICITY_H

#include "lpm.h"
#include "unit_cell.h"

void switchStateV(int conv_flag, UnitCell cell);

void computeCab(UnitCell cell);
void computeBondForceElastic(int i, UnitCell cell);
void computeBondForceJ2mixedLinear3D(int ii, UnitCell cell);
void computeBondForceJ2nonlinearIso(int ii, UnitCell cell);
void computeBondForceCPMiehe(int ii, UnitCell cell);
void computeBondForceIncrementalUpdating(int ii, UnitCell cell);
void computeBondForceJ2energyReturnMap(int ii, int load_indicator, UnitCell cell);

int updateBrittleDamage(const char *dataName, int tstep, int nbreak);
int updateDuctileDamageBwiseLocal(const char *dataName, int tstep, UnitCell cell);
int updateDuctileDamagePwiseLocal(const char *dataName, int tstep, UnitCell cell);
int updateDuctileDamageBwiseNonlocal(const char *dataName, int tstep, UnitCell cell);
int updateDuctileDamagePwiseNonlocal(const char *dataName, int tstep, UnitCell cell);

void updateCrack(UnitCell cell);

#endif

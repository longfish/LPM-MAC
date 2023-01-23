#pragma once
#ifndef PLASTICITY_H
#define PLASTICITY_H

#include "lpm.h"
struct UnitCell;

void switchStateV(int conv_flag, struct UnitCell cell);

void computeCab(struct UnitCell cell);
void computeBondForceElastic(int i, struct UnitCell cell);
void computeBondForceJ2mixedLinear3D(int ii, struct UnitCell cell);
void computeBondForceJ2nonlinearIso(int ii, struct UnitCell cell);
void computeBondForceCPMiehe(int ii, struct UnitCell cell);
void computeBondForceIncrementalUpdating(int ii, struct UnitCell cell);
void computeBondForceJ2energyReturnMap(int ii, int load_indicator, struct UnitCell cell);

int updateBrittleDamage(const char *dataName, int tstep, int nbreak);
int updateDuctileDamageBwiseLocal(const char *dataName, int tstep, struct UnitCell cell);
int updateDuctileDamagePwiseLocal(const char *dataName, int tstep, struct UnitCell cell);
int updateDuctileDamageBwiseNonlocal(const char *dataName, int tstep, struct UnitCell cell);
int updateDuctileDamagePwiseNonlocal(const char *dataName, int tstep, struct UnitCell cell);

void updateCrack(struct UnitCell cell);

#endif

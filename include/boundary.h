#pragma once
#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "lpm.h"

struct DispBCs;
struct ForceBCs;
struct UnitCell;

void setDispBC_stiffnessUpdate2D(struct UnitCell cell);
void setDispBC_stiffnessUpdate3D(struct UnitCell cell);

void setDispBC(int nboundDisp, struct DispBCs* dBP, struct UnitCell cell);
void setForceBC(int nboundForce, struct ForceBCs* fBP, struct UnitCell cell);

#endif

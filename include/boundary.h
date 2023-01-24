#pragma once
#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "lpm.h"
#include "unit_cell.h"

struct DispBCs;
struct ForceBCs;

void setDispBC_stiffnessUpdate2D(UnitCell cell);
void setDispBC_stiffnessUpdate3D(UnitCell cell);

void setDispBC(int nboundDisp, struct DispBCs* dBP, UnitCell cell);
void setForceBC(int nboundForce, struct ForceBCs* fBP, UnitCell cell);

#endif

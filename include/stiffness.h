#pragma once
#ifndef STIFFNESS_H
#define STIFFNESS_H

#include "lpm.h"
#include "unit_cell.h"

void calcKnTv(int ntype, UnitCell cell);
void updateRR(UnitCell cell);
void calcStiffness2DFiniteDifference(int plmode, UnitCell cell);
void calcStiffness3DFiniteDifference(int plmode, UnitCell cell);

#endif

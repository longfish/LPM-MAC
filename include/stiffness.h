#pragma once
#ifndef STIFFNESS_H
#define STIFFNESS_H

#include "lpm.h"

struct UnitCell;

void calcKnTv(int ntype, struct UnitCell cell);
void updateRR(struct UnitCell cell);
void calcStiffness2DFiniteDifference(int plmode, struct UnitCell cell);
void calcStiffness3DFiniteDifference(int plmode, struct UnitCell cell);

#endif

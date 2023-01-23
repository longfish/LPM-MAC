#pragma once
#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include "lpm.h"
struct UnitCell;

void searchNormalNeighbor(struct UnitCell cell);
int searchAFEMNeighbor(struct UnitCell cell);

#endif

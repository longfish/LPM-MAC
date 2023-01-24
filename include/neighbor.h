#pragma once
#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include "lpm.h"
#include "unit_cell.h"

void searchNormalNeighbor(UnitCell cell);
int searchAFEMNeighbor(UnitCell cell);

#endif

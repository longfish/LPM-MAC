#pragma once
#ifndef UNIT_CELL_H
#define UNIT_CELL_H

#include "lpm.h"

class UnitCell
{
public:
    int lattice, dim;
    int nneighbors, nneighbors1, nneighbors2;
    int nneighbors_AFEM, nneighbors_AFEM1, nneighbors_AFEM2;
    double radius, neighbor1_cutoff, neighbor2_cutoff;
    double particle_volume; /* volume of unit cell */
    double el_mapping[NDIM]; /* mapping from elasticity to Kn, Tv */

    UnitCell(){}
    UnitCell(int p_lattice, double p_radius);
};

#endif
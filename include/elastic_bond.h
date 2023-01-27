#pragma once
#ifndef ELASTIC_BOND_H
#define ELASTIC_BOND_H

#include <vector>
#include <array>

#include "lpm.h"
#include "unit_cell.h"
#include "particle.h"

template <int nlayer>
class ElasticBond
{
public:
    int nb, nb_conn;                                     // number of bonds and connections
    double x, y, z;                                      // coordinates
    double damage_visual;                                // damage value for visualization
    UnitCell cell;                                       // unit cell
    std::array<double, NDIM> disp;                       // displacement vector
    std::array<double, NDIM> P_in, P_ex;                 // internal and external particle force
    std::array<std::vector<double>, nlayer> bond_layers; // an array that store n layers of bonds
    std::array<double, NDIM> stress, strain;             // stress and strain tensor

    ElasticBond();
};

#endif
#pragma once
#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <array>

#include "lpm.h"
#include "unit_cell.h"
#include "bond.h"

template <int nlayer>
class Bond;

template <int nlayer>
class Particle
{
public:
    int nb, nb_conn;                                     // number of bonds and connections
    double x, y, z;                                      // coordinates
    double damage_visual;                                // damage value for visualization
    UnitCell cell;                                       // unit cell
    std::array<double, NDIM> disp;                       // displacement vector
    std::array<double, NDIM> P_in, P_ex;                 // internal and external particle force
    std::array<std::vector<Bond<nlayer>>, nlayer> bond_layers; // an array that store n layers of bonds
    std::array<double, NDIM> stress, strain;             // stress and strain tensor

    //Particle();
};

#endif
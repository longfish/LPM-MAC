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
protected:
    static int _ID;

public:
    int id{0};                                                 // identifier of the particle
    int lattice{0};                                            // lattice type of the current particle
    double radius{0.0};                                        // radius of current particle
    int nb{0}, nconn{0};                                       // number of bonds and connections
    double damage_visual{0.};                                  // damage value for visualization
    int K_pt{0};                                               // number of conn larger than (or equal to) its own index
    UnitCell cell;                                             // unit cell
    std::array<double, NDIM> xyz, Pin, Pex;                    // xyz, internal and external particle force
    std::array<std::vector<Bond<nlayer>>, nlayer> bond_layers; // an array that store n layers of bonds
    std::vector<Bond<nlayer>> bonds;                           // vector that stores all bonds
    std::vector<Particle<nlayer>> conn;                        // store all connections of the particle
    std::array<double, NDIM> stress, strain;                   // stress and strain tensor

    Particle() { id = _ID++; /* id starts from 0 */ }
    Particle(double p_x, double p_y, double p_z, int p_lattice, double p_radius);
    Particle(double p_x, double p_y, double p_z, UnitCell p_cell);
    Particle(const Particle<nlayer> &A);
    Particle<nlayer> &operator=(const Particle<nlayer> &A);

    void moveTo(double new_x, double new_y, double new_z);
    double distanceTo(const Particle<nlayer> &A);
};

template <int nlayer>
int Particle<nlayer>::_ID = 0;

#endif
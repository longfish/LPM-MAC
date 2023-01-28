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
    int id;                                                    // identifier of the particle
    int lattice;                                               // lattice type of the current particle
    double radius;                                             // radius of current particle
    int K_pt{0};                                               // number of conn larger than (or equal to) its own index
    int nb{0}, nconn{0};                                       // number of bonds and connections
    double x, y, z;                                            // coordinates
    double m_damage_visual{0.};                                // damage value for visualization
    UnitCell cell;                                             // unit cell
    std::array<double, NDIM> disp;                             // displacement vector
    std::array<double, NDIM> Pin, Pex;                         // internal and external particle force
    std::array<std::vector<Bond<nlayer>>, nlayer> bond_layers; // an array that store n layers of bonds
    std::vector<Particle<nlayer>> conn;                        // store all connections of the particle
    std::array<double, NDIM> stress, strain;                   // stress and strain tensor

    Particle(double p_x, double p_y, double p_z, int p_lattice, double p_radius);
    Particle(double p_x, double p_y, double p_z, UnitCell p_cell);
    Particle(const Particle<nlayer> &A);
    Particle<nlayer> &operator=(const Particle<nlayer> &A);
};

template <int nlayer>
int Particle<nlayer>::_ID = 0;

template <int nlayer>
Particle<nlayer>::Particle(double p_x, double p_y, double p_z, int p_lattice, double p_radius)
{
    x = p_x;
    y = p_y;
    z = p_z;
    radius = p_radius;
    lattice = p_lattice;
    cell = UnitCell(p_lattice, p_radius);
    id = _ID++; // id starts from 0
}

template <int nlayer>
Particle<nlayer>::Particle(double p_x, double p_y, double p_z, UnitCell p_cell)
{
    x = p_x;
    y = p_y;
    z = p_z;
    cell = p_cell;
    id = _ID++; // id starts from 0
}

template <int nlayer>
Particle<nlayer>::Particle(const Particle<nlayer> &A)
{
    x = A.x;
    y = A.y;
    z = A.z;
    radius = A.radius;
    lattice = A.lattice;
    cell = A.cell;
    id = A.id;
}

template <int nlayer>
Particle<nlayer> &Particle<nlayer>::operator=(const Particle<nlayer> &A)
{
    x = A.x;
    y = A.y;
    z = A.z;
    radius = A.radius;
    lattice = A.lattice;
    cell = A.cell;
    id = A.id;
    return (*this);
}

#endif
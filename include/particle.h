#pragma once
#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <array>

#include "lpm.h"
#include "unit_cell.h"
#include "bond.h"

// template <int nlayer>
// class Bond;

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
    std::array<std::vector<Bond<nlayer>*>, nlayer> bond_layers; // an array that store n layers of bonds
    std::vector<Particle<nlayer>*> neighbors;                   // vector that stores all particles that form bonds
    std::vector<Particle<nlayer>*> conns;                       // store all connections of the particle
    std::array<double, NDIM> stress, strain;                   // stress and strain tensor

    Particle() { id = _ID++; /* id starts from 0 */ }
    Particle(const double &p_x, const double &p_y, const double &p_z, const int &p_lattice, const double &p_radius);
    Particle(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell);
    Particle(const Particle<nlayer> &A);
    Particle<nlayer> &operator=(const Particle<nlayer> &A);

    void moveTo(const double &new_x, const double &new_y, const double &new_z);
    double distanceTo(const Particle<nlayer> &A);
    bool operator==(const Particle<nlayer> &other);
};

template <int nlayer>
int Particle<nlayer>::_ID = 0;

template <int nlayer>
bool Particle<nlayer>::operator==(const Particle<nlayer> &other)
{
    return id == other.id;
}

template <int nlayer>
Particle<nlayer>::Particle(const double &p_x, const double &p_y, const double &p_z, const int &p_lattice, const double &p_radius)
{
    xyz[0] = p_x;
    xyz[1] = p_y;
    xyz[2] = p_z;
    radius = p_radius;
    lattice = p_lattice;
    cell = UnitCell(p_lattice, p_radius);
    id = _ID++;
}

template <int nlayer>
Particle<nlayer>::Particle(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell)
{
    xyz[0] = p_x;
    xyz[1] = p_y;
    xyz[2] = p_z;
    cell = p_cell;
    id = _ID++;
}

template <int nlayer>
Particle<nlayer>::Particle(const Particle<nlayer> &A)
{
    xyz[0] = A.xyz[0];
    xyz[1] = A.xyz[1];
    xyz[2] = A.xyz[2];
    radius = A.radius;
    lattice = A.lattice;
    cell = A.cell;
    id = A.id;
}

template <int nlayer>
Particle<nlayer> &Particle<nlayer>::operator=(const Particle<nlayer> &A)
{
    xyz[0] = A.xyz[0];
    xyz[1] = A.xyz[1];
    xyz[2] = A.xyz[2];
    radius = A.radius;
    lattice = A.lattice;
    cell = A.cell;
    id = A.id;
    return (*this);
}

template <int nlayer>
void Particle<nlayer>::moveTo(const double &new_x, const double &new_y, const double &new_z)
{
    xyz[0] = new_x;
    xyz[1] = new_y;
    xyz[2] = new_z;
}

template <int nlayer>
double Particle<nlayer>::distanceTo(const Particle<nlayer> &A)
{
    return sqrt(pow((xyz[0] - A.xyz[0]), 2) + pow((xyz[1] - A.xyz[1]), 2) + pow((xyz[2] - A.xyz[2]), 2));
}

#endif
#pragma once
#ifndef PARTICLE_H
#define PARTICLE_H

#include <algorithm>
#include <vector>
#include <array>

#include "lpm.h"
#include "unit_cell.h"
#include "bond.h"

template <int nlayer>
class Particle
{
protected:
    static int _ID;

public:
    int id{0};                                                               // identifier of the particle
    int type{0};                                                             // particle type which is needed for boundary condition settings
    int K_pt{0};                                                             // matrix pointer, number of conn larger than (or equal to) its own index
    int nb{0}, nconn{0};                                                     // number of bonds and connections
    double damage_visual{0.};                                                // damage value for visualization
    UnitCell cell{0, 0.0};                                                   // unit cell
    std::array<double, nlayer> dLe_total, TdLe_total, ddL_total, TddL_total; // volumetric bond measure
    std::array<double, NDIM> xyz, xyz_initial, xyz_last;                     // particle coordinates
    std::array<double, NDIM> Pin, Pex;                                       // internal and external particle force
    std::array<std::vector<Bond<nlayer> *>, nlayer> bond_layers;             // an array that store n layers of bonds
    std::vector<Particle<nlayer> *> neighbors;                               // vector that stores all particles that form bonds
    std::vector<Particle<nlayer> *> conns;                                   // store all connections of the particle
    std::array<double, NDIM * NDIM> stress, strain;                          // stress and strain tensor

    Particle() { id = _ID++; /* id starts from 0 */ }
    Particle(const double &p_x, const double &p_y, const double &p_z, const int &p_lattice, const double &p_radius);
    Particle(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell);

    void moveTo(const double &new_x, const double &new_y, const double &new_z);
    void updateParticleForce();
    void updateNeighborsGeometry();
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
    xyz = {p_x, p_y, p_z};
    xyz_initial = xyz;
    cell = UnitCell(p_lattice, p_radius);
    id = _ID++;
}

template <int nlayer>
Particle<nlayer>::Particle(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell)
{
    xyz = {p_x, p_y, p_z};
    xyz_initial = xyz;
    cell = p_cell;
    id = _ID++;
}

template <int nlayer>
void Particle<nlayer>::moveTo(const double &new_x, const double &new_y, const double &new_z)
{
    xyz_last = xyz;
    xyz = {new_x, new_y, new_z};
}

template <int nlayer>
void Particle<nlayer>::updateNeighborsGeometry()
{
    // update all the neighbors information
    for (int i = 0; i < nlayer; ++i)
    {
        dLe_total[i] = 0, TdLe_total[i] = 0;
        ddL_total[i] = 0, TddL_total[i] = 0;
        for (Bond<nlayer> *bd : bond_layers[i])
        {
            bd->updatebGeometry();
            dLe_total += bd->dLe;
            TdLe_total += bd->Tv * bd->dLe;
            ddL_total += bd->ddL;
            TddL_total += bd->Tv * bd->ddL;
        }
    }
}

template <int nlayer>
void Particle<nlayer>::updateParticleForce()
{
    // sum up all the particle forces
    Pin = {0., 0., 0.};
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : bond_layers[i])
        {
            // find the opposite bond
            auto op_bd = std::find_if(bd->p2->bond_layers[i].begin(), bd->p2->bond_layers[i].end(),
                                      [&](Bond<nlayer> *b)
                                      { return b->p2->id == id; });
            Pin[0] += bd->csx * 0.5 * (bd->bforce + (*op_bd)->bforce);
            Pin[1] += bd->csy * 0.5 * (bd->bforce + (*op_bd)->bforce);
            Pin[2] += bd->csz * 0.5 * (bd->bforce + (*op_bd)->bforce);
        }
    }
}

template <int nlayer>
double Particle<nlayer>::distanceTo(const Particle<nlayer> &A)
{
    return sqrt(pow((xyz[0] - A.xyz[0]), 2) + pow((xyz[1] - A.xyz[1]), 2) + pow((xyz[2] - A.xyz[2]), 2));
}

#endif
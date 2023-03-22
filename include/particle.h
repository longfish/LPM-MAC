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
    int id{0};                // identifier of the particle, id starts from 0
    int type{0};              // particle type which is needed to identify phases
    int nconn_largeq{0};      // matrix pointer, number of conn larger than (or equal to) its own index
    int nb{0}, nconn{0};      // number of bonds and connections
    double damage_visual{0.}; // damage value for visualization
    UnitCell cell;            // unit cell

    std::array<int, NDIM> disp_constraint{0};                    // disp BC indicator, 1 means disp BC applied, 0 otherwise
    std::array<double, nlayer> dLe_total, ddL_total;             // volumetric bond measure
    std::array<double, nlayer> cs_sumx, cs_sumy, cs_sumz;        // volumetric bond measure
    std::array<double, NDIM> xyz, xyz_initial, xyz_last;         // particle coordinates
    std::array<double, NDIM> Pin{0}, Pex{0};                     // internal and external particle force
    std::array<std::vector<Bond<nlayer> *>, nlayer> bond_layers; // an array that store n layers of bonds
    std::vector<Particle<nlayer> *> neighbors;                   // vector that stores all particles that form bonds
    std::vector<Particle<nlayer> *> conns;                       // all connections of the particle (include self)
    std::array<double, NDIM * NDIM> stress{0}, strain{0};        // stress and strain tensor

    Particle(const int &p_id) : cell{LatticeType::SimpleCubic3D, 0} { id = p_id; };
    Particle(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell, const int &p_type);
    Particle(const double &p_x, const double &p_y, const double &p_z, const LatticeType &p_lattice, const double &p_radius);
    Particle(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell);

    void moveTo(const double &new_x, const double &new_y, const double &new_z);
    void moveTo(const std::array<double, NDIM> &new_xyz);
    void moveBy(const std::array<double, NDIM> &dxyz);

    bool hasAFEMneighbor(Particle<nlayer> *pj, int layer);

    void updateParticleForce();
    void updateBondsGeometry();
    void updateBondsForce();
    void resumeParticle();
    double distanceTo(Particle<nlayer> *A);
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
bool Particle<nlayer>::hasAFEMneighbor(Particle<nlayer> *pj, int layer)
{
    for (Bond<nlayer> *bd1 : bond_layers[layer])
    {
        if ((bd1->p2->id == pj->id) && !(bd1->broken))
            return true;
        for (Bond<nlayer> *bd2 : bd1->p2->bond_layers[layer])
        {
            if ((bd2->p2->id == pj->id) && !(bd2->broken))
                return true;
        }
    }
    return false;
}

template <int nlayer>
Particle<nlayer>::Particle(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell, const int &p_type)
    : cell{p_cell}
{
    xyz = {p_x, p_y, p_z};
    xyz_initial = xyz;
    type = p_type;
    id = _ID++;
}

template <int nlayer>
Particle<nlayer>::Particle(const double &p_x, const double &p_y, const double &p_z, const LatticeType &p_lattice, const double &p_radius)
    : cell{p_lattice, p_radius}
{
    xyz = {p_x, p_y, p_z};
    xyz_initial = xyz;
    id = _ID++;
}

template <int nlayer>
Particle<nlayer>::Particle(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell)
    : cell{p_cell}
{
    xyz = {p_x, p_y, p_z};
    xyz_initial = xyz;
    id = _ID++;
}

template <int nlayer>
void Particle<nlayer>::moveTo(const double &new_x, const double &new_y, const double &new_z)
{
    xyz_last = xyz;
    xyz = {new_x, new_y, new_z};
}

template <int nlayer>
void Particle<nlayer>::moveTo(const std::array<double, NDIM> &new_xyz)
{
    xyz_last = xyz;
    xyz = new_xyz;
}

template <int nlayer>
void Particle<nlayer>::moveBy(const std::array<double, NDIM> &dxyz)
{
    xyz_last = xyz;
    xyz = {xyz[0] + dxyz[0], xyz[1] + dxyz[1], xyz[2] + dxyz[2]};
}

template <int nlayer>
void Particle<nlayer>::updateBondsGeometry()
{
    // update all the neighbors information
    for (int i = 0; i < nlayer; ++i)
    {
        dLe_total[i] = 0, ddL_total[i] = 0;
        cs_sumx[i] = 0, cs_sumy[i] = 0, cs_sumz[i] = 0;
        for (Bond<nlayer> *bd : bond_layers[i])
        {
            bd->updatebGeometry();
            dLe_total[i] += bd->dLe;
            ddL_total[i] += bd->ddL;
            cs_sumx[i] += bd->csx;
            cs_sumy[i] += bd->csy;
            cs_sumz[i] += bd->csz;
        }
    }
}

template <int nlayer>
void Particle<nlayer>::updateBondsForce()
{
    // update all the bond forces
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : bond_layers[i])
            bd->updatebForce();
    }
}

template <int nlayer>
void Particle<nlayer>::resumeParticle()
{
    // used for calculating stiffness matrix
    updateBondsGeometry();
    updateBondsForce();
}

template <int nlayer>
void Particle<nlayer>::updateParticleForce()
{
    // sum up all bond forces
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
double Particle<nlayer>::distanceTo(Particle<nlayer> *A)
{
    return sqrt(pow((xyz[0] - A->xyz[0]), 2) + pow((xyz[1] - A->xyz[1]), 2) + pow((xyz[2] - A->xyz[2]), 2));
}

#endif
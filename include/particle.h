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
    int id{0};                                                       // identifier of the particle, id starts from 0
    int type{0};                                                     // particle type which is needed to identify phases
    int nconn_largeq{0};                                             // matrix pointer, number of conn larger than (or equal to) its own index
    int nb{0}, nconn{0};                                             // number of bonds and connections
    int ncycle_jump{0};                                              // cycle jumping numbers for fatigue simulation
    double damage_visual{0.};                                        // damage value for visualization
    UnitCell cell;                                                   // unit cell
    std::vector<double> state_var, state_var_last;                   // vectors that store state variables
    std::array<int, NDIM> disp_constraint{0};                        // disp BC indicator, 1 means disp BC applied, 0 otherwise
    std::array<double, nlayer> dLe_total, cs_sumx, cs_sumy, cs_sumz; // volumetric bond measure
    std::array<double, NDIM> xyz, xyz_initial, xyz_last;             // particle coordinates
    std::array<double, NDIM> Pin{0}, Pex{0};                         // internal and external particle force
    std::array<std::vector<Bond<nlayer> *>, nlayer> bond_layers;     // an array that store n layers of bonds
    std::vector<Particle<nlayer> *> neighbors;                       // vector that stores all particles that form bonds
    std::vector<Particle<nlayer> *> conns;                           // all connections of the particle (include self)
    std::vector<double> stress, strain;                              // stress and strain tensor

    Particle(const int &p_id) : cell{LatticeType::SimpleCubic3D, 0} { id = p_id; };
    Particle(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell, const int &p_type);
    Particle(const double &p_x, const double &p_y, const double &p_z, const LatticeType &p_lattice, const double &p_radius);
    Particle(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell);

    void moveTo(const double &new_x, const double &new_y, const double &new_z);
    void moveTo(const std::array<double, NDIM> &new_xyz);
    void moveBy(const std::array<double, NDIM> &dxyz);

    bool hasAFEMneighbor(Particle<nlayer> *pj, int layer);
    bool operator==(const Particle<nlayer> &other);
    double distanceTo(Particle<nlayer> *A);

    void updateParticleForce();
    void updateParticleStress();
    void updateBondsGeometry();
    void updateParticleDamageVisual();
    void storeParticleStateVariables();
    void resetParticleStateVariables();
    void clearParticleStateVariables();
    void resumeParticleState();

    // virtual functions that can be inherited
    virtual void updateBondsForce() {}
    virtual int calcNCycleJump() { return 0; }
    virtual bool updateParticleStateVariables() { return false; }
    virtual bool updateParticleBrokenBonds() { return false; };
};

template <int nlayer>
int Particle<nlayer>::_ID = 0;

template <int nlayer>
bool Particle<nlayer>::operator==(const Particle<nlayer> &other)
{
    return id == other.id;
}

template <int nlayer>
void Particle<nlayer>::storeParticleStateVariables()
{
    // set the last converged state variables to be the current one
    for (int i = 0; i < state_var.size(); ++i)
        state_var_last[i] = state_var[i];

    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : bond_layers[i])
        {
            bd->dis_last = bd->dis;
            bd->bdamage_last = bd->bdamage;
            bd->bforce_last = bd->bforce;
            bd->dLp_last = bd->dLp;
        }
    }
}

template <int nlayer>
void Particle<nlayer>::clearParticleStateVariables()
{
    // clear the current particle state variables (useful in fatigue modeling)
    for (int i = 0; i < state_var.size(); ++i)
    {
        state_var_last[i] = 0;
        state_var[i] = 0;
    }
}

template <int nlayer>
void Particle<nlayer>::resetParticleStateVariables()
{
    // reset back the current state variables to the last converged one
    for (int i = 0; i < state_var.size(); ++i)
        state_var[i] = state_var_last[i];

    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : bond_layers[i])
        {
            bd->dis = bd->dis_last;
            bd->bdamage = bd->bdamage_last;
            bd->bforce = bd->bforce_last;
            bd->dLp = bd->dLp_last;
            // if (id == 5145)
            //     printf("bstrain %f \n", bd->bstrain);
        }
    }
}

template <int nlayer>
void Particle<nlayer>::updateParticleDamageVisual()
{
    // if (id == 5145)
    //     printf("D: %f \n", state_var[0]);

    damage_visual = 0;
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : bond_layers[i])
            damage_visual += bd->bdamage;
    }

    damage_visual = damage_visual / (double)neighbors.size();
}

template <int nlayer>
void Particle<nlayer>::updateBondsGeometry()
{
    // update all the neighbors information
    for (int i = 0; i < nlayer; ++i)
    {
        dLe_total[i] = 0;
        cs_sumx[i] = 0, cs_sumy[i] = 0, cs_sumz[i] = 0;
        for (Bond<nlayer> *bd : bond_layers[i])
        {
            bd->updatebGeometry();
            dLe_total[i] += bd->dLe;
            cs_sumx[i] += bd->csx;
            cs_sumy[i] += bd->csy;
            cs_sumz[i] += bd->csz;
        }
    }
}

template <int nlayer>
void Particle<nlayer>::resumeParticleState()
{
    // this function is mainly for calculating stiffness matrix by finite difference
    for (Particle<nlayer> *pjj : conns)
    {
        pjj->updateBondsGeometry(); // update all bond information, e.g., dL, dL_total
        pjj->resetParticleStateVariables();
    }

    for (Particle<nlayer> *pjj : conns)
        pjj->updateBondsForce(); // update all bond forces
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
void Particle<nlayer>::updateParticleStress()
{
    // initialize the tensor
    stress = std::vector<double>(2 * NDIM, 0.0);

    // compute the tensor using average bond force
    double V_m = cell.particle_volume * nb / cell.nneighbors;
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : bond_layers[i])
        {
            auto op_bd = std::find_if(bd->p2->bond_layers[i].begin(), bd->p2->bond_layers[i].end(),
                                      [&](Bond<nlayer> *b)
                                      { return b->p2->id == id; });
            stress[0] += 0.5 / V_m * (bd->dis) * 0.5 * (bd->bforce + (*op_bd)->bforce) * (bd->csx) * (bd->csx);
            stress[1] += 0.5 / V_m * (bd->dis) * 0.5 * (bd->bforce + (*op_bd)->bforce) * (bd->csy) * (bd->csy);
            stress[2] += 0.5 / V_m * (bd->dis) * 0.5 * (bd->bforce + (*op_bd)->bforce) * (bd->csz) * (bd->csz);
            stress[3] += 0.5 / V_m * (bd->dis) * 0.5 * (bd->bforce + (*op_bd)->bforce) * (bd->csy) * (bd->csz);
            stress[4] += 0.5 / V_m * (bd->dis) * 0.5 * (bd->bforce + (*op_bd)->bforce) * (bd->csx) * (bd->csz);
            stress[5] += 0.5 / V_m * (bd->dis) * 0.5 * (bd->bforce + (*op_bd)->bforce) * (bd->csx) * (bd->csy);
        }
    }
}

template <int nlayer>
bool Particle<nlayer>::hasAFEMneighbor(Particle<nlayer> *pj, int layer)
{
    for (Bond<nlayer> *bd1 : bond_layers[layer])
    {
        if (bd1->p2->id == pj->id)
            return true;
        for (Bond<nlayer> *bd2 : bd1->p2->bond_layers[layer])
        {
            if (bd2->p2->id == pj->id)
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
    xyz_last = xyz;
    type = p_type;
    id = _ID++;
}

template <int nlayer>
Particle<nlayer>::Particle(const double &p_x, const double &p_y, const double &p_z, const LatticeType &p_lattice, const double &p_radius)
    : cell{p_lattice, p_radius}
{
    xyz = {p_x, p_y, p_z};
    xyz_initial = xyz;
    xyz_last = xyz;
    id = _ID++;
}

template <int nlayer>
Particle<nlayer>::Particle(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell)
    : cell{p_cell}
{
    xyz = {p_x, p_y, p_z};
    xyz_initial = xyz;
    xyz_last = xyz;
    id = _ID++;
}

template <int nlayer>
void Particle<nlayer>::moveTo(const double &new_x, const double &new_y, const double &new_z)
{
    xyz = {new_x, new_y, new_z};
}

template <int nlayer>
void Particle<nlayer>::moveTo(const std::array<double, NDIM> &new_xyz)
{
    xyz = new_xyz;
}

template <int nlayer>
void Particle<nlayer>::moveBy(const std::array<double, NDIM> &dxyz)
{
    xyz = {xyz[0] + dxyz[0], xyz[1] + dxyz[1], xyz[2] + dxyz[2]};
}

template <int nlayer>
double Particle<nlayer>::distanceTo(Particle<nlayer> *A)
{
    return sqrt(pow((xyz[0] - A->xyz[0]), 2) + pow((xyz[1] - A->xyz[1]), 2) + pow((xyz[2] - A->xyz[2]), 2));
}

#endif
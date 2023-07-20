#pragma once
#ifndef PARTICLE_ELASTIC_H
#define PARTICLE_ELASTIC_H

#include <vector>
#include <array>

#include "lpm.h"
#include "unit_cell.h"
#include "particle.h"
#include "bond.h"

// implement bond-based damage model

template <int nlayer>
class ParticleElastic : public Particle<nlayer>
{
    double critical_bstrain{0};

public:
    ParticleElastic(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell) : Particle<nlayer>{p_x, p_y, p_z, p_cell}
    {
        this->state_var = {0};
        this->state_var_last = {0};
    }
    ParticleElastic(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell, const int &p_type) : Particle<nlayer>{p_x, p_y, p_z, p_cell, p_type}
    {
        this->state_var = {0};
        this->state_var_last = {0};
    }

    void updateBondsForce();
    bool updateParticleStaticDamage();
    void setParticleProperty(double p_nonlocalL, bool is_plane_stress, double p_E, double p_mu, double p_critical_bstrain);
    void setParticleProperty(double p_nonlocalL, double p_C11, double p_C12, double p_C44, double p_critical_bstrain);
};

template <int nlayer>
void ParticleElastic<nlayer>::updateBondsForce()
{
    // for elastic bonds, a trial elastic calculation is enough
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
        {
            bd->bforce = 2. * bd->Kn * bd->dLe + 2. * bd->Tv * this->dLe_total[bd->layer];
            bd->bforce *= (1.0 - bd->bdamage);
        }
    }
}

template <int nlayer>
bool ParticleElastic<nlayer>::updateParticleStaticDamage()
{
    bool any_broken{false};
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
        {
            if (bd->bstrain >= critical_bstrain && abs(bd->bdamage - 1.0) > EPS)
            {
                any_broken = any_broken || true;
                bd->bdamage = 1;
            }
        }
    }
    return any_broken;
}

template <int nlayer>
void ParticleElastic<nlayer>::setParticleProperty(double p_nonlocalL, bool is_plane_stress, double p_E, double p_mu, double p_critical_bstrain)
{
    critical_bstrain = p_critical_bstrain;
    this->nonlocal_L = p_nonlocalL;

    double KnTv[NDIM]{0};
    std::vector<double> Ce(NDIM);
    if (!is_plane_stress)
        Ce = {p_E * (1.0 - p_mu) / (1.0 + p_mu) / (1.0 - 2.0 * p_mu),
              p_E * p_mu / (1.0 + p_mu) / (1.0 - 2.0 * p_mu),
              p_E / 2.0 / (1.0 + p_mu)}; // C11, C12, C44
    else
        Ce = {p_E / (1.0 - p_mu) / (1.0 + p_mu),
              p_E * p_mu / (1.0 + p_mu) / (1.0 - p_mu),
              p_E / 2.0 / (1.0 + p_mu)}; // C11, C12, C44

    if (this->cell.dim == 2)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NDIM, 1, NDIM, 1.0, this->cell.el_mapping.data(), 3, Ce.data(), 1, 0.0, KnTv, 1);
    else
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NDIM, 1, NDIM, this->cell.radius, this->cell.el_mapping.data(), 3, Ce.data(), 1, 0.0, KnTv, 1);

    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
        {
            if (this->cell.lattice == LatticeType::Hexagon2D)
            {
                bd->Kn = KnTv[0];
                bd->Tv = KnTv[1];
            }
            else
            {
                bd->Kn = KnTv[bd->layer]; // layer is 0 or 1
                bd->Tv = KnTv[2];
            }
        }
    }
}

template <int nlayer>
void ParticleElastic<nlayer>::setParticleProperty(double p_nonlocalL, double p_C11, double p_C12, double p_C44, double p_critical_bstrain)
{
    critical_bstrain = p_critical_bstrain;
    this->nonlocal_L = p_nonlocalL;

    double Ce[NDIM]{p_C11, p_C12, p_C44}; // C11, C12, C44
    double KnTv[NDIM]{0};

    if (this->cell.dim == 2)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NDIM, 1, NDIM, 1.0, this->cell.el_mapping.data(), 3, Ce, 1, 0.0, KnTv, 1);
    else
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NDIM, 1, NDIM, this->cell.radius, this->cell.el_mapping.data(), 3, Ce, 1, 0.0, KnTv, 1);

    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
        {
            if (this->cell.lattice == LatticeType::Hexagon2D)
            {
                bd->Kn = KnTv[0];
                bd->Tv = KnTv[1];
            }
            else
            {
                bd->Kn = KnTv[bd->layer]; // layer is 0 or 1
                bd->Tv = KnTv[2];
            }
        }
    }
}

#endif
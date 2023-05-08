#pragma once
#ifndef PARTICLE_ELASTIC_DAMAGE_H
#define PARTICLE_ELASTIC_DAMAGE_H

#include <vector>
#include <array>

#include "lpm.h"
#include "unit_cell.h"
#include "particle.h"
#include "bond.h"

// Elastic plane strain or 3D material
// Update bdamage and bforce after geometry calculation 
// Three particle-wise state variables: [0]damage


template <int nlayer>
class ParticleFatigueHCF : public Particle<nlayer>
{
public:
    double kappa0{0}, alpha{0}, beta{0}; // damage parameters
    double comp_tensile_ratio{1};
    double damage_threshold{0.35}, critical_bstrain{1.2e-3};

    ParticleFatigueHCF(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell) : Particle<nlayer>{p_x, p_y, p_z, p_cell}
    {
        this->state_var = std::vector<double>(3, 0.);
        this->state_var_last = std::vector<double>(3, 0.);
    }

    ParticleFatigueHCF(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell, const int &p_type) : Particle<nlayer>{p_x, p_y, p_z, p_cell, p_type}
    {
        this->state_var = std::vector<double>(3, 0.);
        this->state_var_last = std::vector<double>(3, 0.);
    }

    double calcEqStrain();

    bool updateParticleStateVariables();
    bool updateParticleBrokenBonds();

    void updateBondsForce();
    void setParticleProperty(double p_E, double p_mu, double p_kappa0, double p_alpha, double p_beta, double p_comp_tensile_ratio);
    void setParticleProperty(double p_C11, double p_C12, double p_C44, double p_kappa0, double p_alpha, double p_beta, double p_comp_tensile_ratio);
};

template <int nlayer>
bool ParticleFatigueHCF<nlayer>::updateParticleBrokenBonds()
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
                //--(this->nb);
            }
        }
    }
    return any_broken;
}

template <int nlayer>
void ParticleFatigueHCF<nlayer>::updateBondsForce()
{
    // calculate the current bond force using current damage value
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
        {
            bd->bforce = 2. * bd->Kn * bd->dLe + 2. * bd->Tv * this->dLe_total[bd->layer];           // trial elastic bforce
            bd->bdamage = std::max(bd->bdamage, std::max(this->state_var[0], bd->p2->state_var[0])); // update the bond-wise damage
            bd->bforce *= (1.0 - bd->bdamage);
        }
    }
}

template <int nlayer>
double ParticleFatigueHCF<nlayer>::calcEqStrain()
{

}

template <int nlayer>
bool ParticleFatigueHCF<nlayer>::updateParticleStateVariables()
{

}

template <int nlayer>
void ParticleFatigueHCF<nlayer>::setParticleProperty(double p_E, double p_mu, double p_kappa0, double p_alpha, double p_beta, double p_comp_tensile_ratio)
{
    kappa0 = p_kappa0;
    alpha = p_alpha;
    beta = p_beta;
    comp_tensile_ratio = p_comp_tensile_ratio;
    this->state_var[2] = kappa0;
    this->state_var_last[2] = kappa0;

    double Ce[NDIM]{p_E * (1.0 - p_mu) / (1.0 + p_mu) / (1.0 - 2.0 * p_mu),
                    p_E * p_mu / (1.0 + p_mu) / (1.0 - 2.0 * p_mu),
                    p_E / 2.0 / (1.0 + p_mu)}; // C11, C12, C44
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

template <int nlayer>
void ParticleFatigueHCF<nlayer>::setParticleProperty(double p_C11, double p_C12, double p_C44, double p_kappa0, double p_alpha, double p_beta, double p_comp_tensile_ratio)
{
    kappa0 = p_kappa0;
    alpha = p_alpha;
    beta = p_beta;
    comp_tensile_ratio = p_comp_tensile_ratio;
    this->state_var[2] = kappa0;
    this->state_var_last[2] = kappa0;

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
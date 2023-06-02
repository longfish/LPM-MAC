#pragma once
#ifndef PARTICLE_FATIGUE_HCF_H
#define PARTICLE_FATIGUE_HCF_H

#include <vector>
#include <array>

#include "lpm.h"
#include "unit_cell.h"
#include "particle.h"
#include "bond.h"

// Elastic plane strain or 3D material
// Update bdamage and bforce after geometry calculation
// Three particle-wise state variables: [0]dDdN(G), [1]min_energy, [2]max_energy, [3]his_max_energy

template <int nlayer>
class ParticleFatigueHCF : public Particle<nlayer>
{
public:
    double A{0}, B{0}, k{0};
    double E{0}, mu{0}, fatigue_limit_ratio{1.108};
    double damage_threshold{1};

    ParticleFatigueHCF(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell) : Particle<nlayer>{p_x, p_y, p_z, p_cell}
    {
        this->state_var = {0, INT_MAX, 0, 0};
        this->state_var_last = {0, INT_MAX, 0, 0};
    }

    ParticleFatigueHCF(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell, const int &p_type) : Particle<nlayer>{p_x, p_y, p_z, p_cell, p_type}
    {
        this->state_var = {0, INT_MAX, 0, 0};
        this->state_var_last = {0, INT_MAX, 0, 0};
    }

    double calcEqEnergy();

    void updateParticleStateVariables();
    bool updateParticleBrokenBonds();
    bool updateParticleFatigueDamage(double &dt);

    void updateBondsForce();
    void setParticleProperty(double p_E, double p_mu, double p_A, double p_B, double p_k, double p_d_thres, double p_f_lmt_ratio);
};

template <int nlayer>
bool ParticleFatigueHCF<nlayer>::updateParticleBrokenBonds()
{
    bool any_broken{false};
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
        {
            if (abs(bd->bdamage - 1.0) > EPS)
            {
                any_broken = any_broken || true;
                bd->bdamage = 1;
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
            bd->bforce = 2. * bd->Kn * bd->dLe + 2. * bd->Tv * this->dLe_total[bd->layer]; // trial elastic bforce
            // bd->bdamage = std::max(bd->bdamage, std::max(this->damage, bd->p2->damage));   // update the bond-wise damage
            if (abs(this->damage - 1.0) < EPS || abs(bd->p2->damage - 1.0) < EPS)
                bd->bdamage = 1; // update the bond-wise damage
            bd->bforce *= (1.0 - bd->bdamage);
        }
    }
}

template <int nlayer>
double ParticleFatigueHCF<nlayer>::calcEqEnergy()
{
    // compute dilatational stretch for each layer
    std::vector<double> dLe_dil(nlayer, 0);
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
            dLe_dil[i] += bd->dLe;
        dLe_dil[i] /= this->bond_layers[i].size();
    }

    // compute energy terms
    double V_m = this->cell.particle_volume * this->nb / this->cell.nneighbors; // modified particle volume
    double energy_total{0}, energy_dil{0}, energy_dis{0}, energy_uni{0};
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
        {
            double dLe_dis = bd->dLe - dLe_dil[i];
            energy_dis += 0.5 * (bd->Kn * dLe_dis * dLe_dis) / V_m;
            energy_total += 0.5 * (bd->Kn * bd->dLe * bd->dLe) / V_m; // first add spring energy
        }
        energy_total += 0.5 * this->TdLe_total[i] * this->dLe_total[i] / V_m; // then add volumetric energy
    }
    energy_dil = energy_total - energy_dis;

    double k = 1 / 3 / (1 - 2 * mu) * pow(503 / 570, 2); // 276/434
    if (energy_dis > (2 + 2 * mu) * energy_dil / (1 - 2 * mu))
    {
        energy_uni = 3 * energy_dil / (1 - 2 * mu);
        energy_dis -= (2 + 2 * mu) * energy_dil / (1 - 2 * mu);
        energy_dil = 0;
    }
    else
    {
        energy_uni = 3 * energy_dis / (2 + 2 * mu);
        energy_dil -= (1 - 2 * mu) * energy_dis / (2 + 2 * mu);
        energy_dis = 0;
    }

    // return energy_total;
    return fatigue_limit_ratio * energy_dis + energy_uni + k * energy_dil;
}

template <int nlayer>
void ParticleFatigueHCF<nlayer>::updateParticleStateVariables()
{
    // compute equivalent energy and store the maximum and minimum ones
    double curr_energy = calcEqEnergy();

    this->state_var[1] = std::min(curr_energy, this->state_var[1]);
    this->state_var[2] = std::max(curr_energy, this->state_var[2]);
    this->state_var[3] = std::max(this->state_var[2], this->state_var[3]);
}

template <int nlayer>
bool ParticleFatigueHCF<nlayer>::updateParticleFatigueDamage(double &dt)
{
    // calculate damage rate
    if (this->state_var[2] - this->state_var[1] > 0)
        this->state_var[0] = A * this->state_var[3] * exp(k * this->damage) * (this->state_var[2] - this->state_var[1]);
    else
        this->state_var[0] = 0;

    // update damage value
    this->damage += this->state_var[0] * dt;
    if (this->damage >= damage_threshold)
        this->damage = 1;

    // clear the unnecessary state variables
    this->state_var[0] = 0;
    this->state_var[1] = INT_MAX;
    this->state_var[2] = 0;

    return abs(this->damage - this->damage_last) > 1e-3; // return true if damage is updated noticeably
}

template <int nlayer>
void ParticleFatigueHCF<nlayer>::setParticleProperty(double p_E, double p_mu, double p_A, double p_B, double p_k, double p_d_thres, double p_f_lmt_ratio)
{
    E = p_E;
    mu = p_mu;
    A = p_A;
    B = p_B;
    k = p_k;
    damage_threshold = p_d_thres;
    fatigue_limit_ratio = p_f_lmt_ratio;

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

#endif
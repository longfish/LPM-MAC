#pragma once
#ifndef PARTICLE_ELASTIC_DAMAGE_H
#define PARTICLE_ELASTIC_DAMAGE_H

#include <vector>
#include <array>

#include "lpm.h"
#include "unit_cell.h"
#include "particle.h"
#include "bond.h"

// Update bdamage and bforce after calculation
// Two particle-wise state variables: [0]strain energy

// update bond force
// 1. update the particle geometry (loop for all particles)
// 2. update the particle damage (loop for all particles)
//     i. calculate the strain energy and store it as a state variable - state_var
//     ii. calculate the damage dot - Ddot
//     iii. using explicit Euler, D = D + Ddot*delta_energy
// 3. update bforce (after calculating the bdamage of all bonds)

template <int nlayer>
class ParticleElasticDamage : public Particle<nlayer>
{
public:
    double k0{0}, k1{0}; // damage parameters
    double damage_threshold{0}, comp_tensile_ratio{0};

    ParticleElasticDamage(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell) : Particle<nlayer>{p_x, p_y, p_z, p_cell}
    {
        this->state_var = std::vector<double>(1, 0.);
        this->state_var_last = std::vector<double>(1, 0.);
    }

    ParticleElasticDamage(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell, const int &p_type) : Particle<nlayer>{p_x, p_y, p_z, p_cell, p_type}
    {
        this->state_var = std::vector<double>(1, 0.);
        this->state_var_last = std::vector<double>(1, 0.);
    }

    double calcStrainEnergyDensity();
    double calcEqStrain();

    void updateParticleStateVariables();
    bool updateParticleStaticDamage();
    bool updateParticleBrokenBonds();

    void updateBondsForce();
    void setParticleProperty(double p_nonlocalL, bool is_plane_stress, double p_E, double p_mu, double p_k0, double p_k1, double p_ct_ratio, double p_d_thres);
    void setParticleProperty(double p_nonlocalL, double p_C11, double p_C12, double p_C44, double p_k0, double p_k1, double p_ct_ratio, double p_d_thres);
};

template <int nlayer>
bool ParticleElasticDamage<nlayer>::updateParticleBrokenBonds()
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
                //--(this->nb);
            }
        }
    }
    return any_broken;
}

template <int nlayer>
void ParticleElasticDamage<nlayer>::updateBondsForce()
{
    // calculate the current bond force using current damage value
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
        {
            bd->bforce = 2. * bd->Kn * bd->dLe + 2. * bd->Tv * this->dLe_total[bd->layer]; // trial elastic bforce
            bd->bdamage = std::max(bd->bdamage, std::max(this->damage, bd->p2->damage));   // update the bond-wise damage
            bd->bforce *= (1.0 - bd->bdamage);
        }
    }
}

template <int nlayer>
double ParticleElasticDamage<nlayer>::calcStrainEnergyDensity()
{
    // compute energy terms
    double V_m = this->cell.particle_volume * this->nb / this->cell.nneighbors; // modified particle volume
    double energy_total{0};
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
            energy_total += 0.5 * (bd->Kn * bd->dLe * bd->dLe) / V_m; // first add spring energy

        energy_total += 0.5 * this->TdLe_total[i] * this->dLe_total[i] / V_m; // then add volumetric energy
    }

    return energy_total;
}

double func_g(const double k1)
{
    return 1 / k1;
}

// double func_g(const double k, const double k0, const double undamaged_pt_type)
// {
//     if (k >= undamaged_pt_type)
//         return 0;
//     else
//         return k0 * undamaged_pt_type / (k0 + undamaged_pt_type) * pow(k, -2); // the k here is kappa
// }

// double func_g(const double k, const double k0, const double a, const double b)
// {
//     return k0 / k / k * (1. + a * ((1. + b * k) * std::exp(b * (k0 - k)) - 1.)); // the k here is kappa
// }

// double func_D(const double k, const double k0, const double a, const double b)
// {
//     return 1.0 - k0 / k * (1.0 - a + a * std::exp(-b * (k - k0)));
// }

double func_eq_strain(const double &comp_tensile_ratio, const double &I1, const double &J2)
{
    double k = comp_tensile_ratio; // the k has different meaning compared with k in func_g
    return (1. - 1. / k) * I1 + 1. / k * std::sqrt((k - 1.) * (k - 1.) * I1 * I1 - k * J2);
}

template <int nlayer>
double ParticleElasticDamage<nlayer>::calcEqStrain()
{
    double I1{0}, J2{0};
    for (int i = 0; i < nlayer; ++i)
    {
        double bl = this->cell.neighbor_cutoff[i]; // bond length of the current layer
        I1 += this->dLe_total[i] / bl;
        for (Bond<nlayer> *bd : this->bond_layers[i])
            J2 += bd->bstrain * bd->bstrain;
    }
    J2 = 0.5 * I1 * I1 / this->nb - 0.5 * J2;
    return func_eq_strain(comp_tensile_ratio, I1, J2);
}

template <int nlayer>
void ParticleElasticDamage<nlayer>::updateParticleStateVariables()
{
    // update local state variables
    this->state_var[0] = calcStrainEnergyDensity();              // local strain measures
    double delta = this->state_var[0] - this->state_var_last[0]; // state_var change
    double f = this->state_var[0] - k1 * this->damage - k0;      // damage surface
    if (f > 0 && delta > 0)
        this->Ddot_local = func_g(k1) * delta; // local Ddot
    else
        this->Ddot_local = 0;
}

template <int nlayer>
bool ParticleElasticDamage<nlayer>::updateParticleStaticDamage()
{
    // update nonlocal damage
    this->damage += this->Ddot_nonlocal;

    if (this->damage >= damage_threshold)
        this->damage = 1;

    return abs(this->damage - this->damage_last) > 1e-3; // return true if damage is updated noticeably
}

template <int nlayer>
void ParticleElasticDamage<nlayer>::setParticleProperty(double p_nonlocalL, bool is_plane_stress, double p_E, double p_mu, double p_k0, double p_k1, double p_ct_ratio, double p_d_thres)
{
    k0 = p_k0;
    k1 = p_k1;
    damage_threshold = p_d_thres;
    comp_tensile_ratio = p_ct_ratio;
    this->nonlocal_L = p_nonlocalL;

    Eigen::VectorXd KnTv(NDIM), Ce(NDIM);
    if (!is_plane_stress)
    {
        Ce(0) = p_E * (1.0 - p_mu) / (1.0 + p_mu) / (1.0 - 2.0 * p_mu);
        Ce(1) = p_E * p_mu / (1.0 + p_mu) / (1.0 - 2.0 * p_mu);
        Ce(2) = p_E / 2.0 / (1.0 + p_mu); // C11, C12, C44
    }
    else
    {
        Ce(0) = p_E / (1.0 - p_mu) / (1.0 + p_mu);
        Ce(1) = p_E * p_mu / (1.0 + p_mu) / (1.0 - p_mu);
        Ce(2) = p_E / 2.0 / (1.0 + p_mu);
    }; // C11, C12, C44

    if (this->cell.dim == 2)
        KnTv = this->cell.el_mapping * Ce;
    else
        KnTv = this->cell.radius * this->cell.el_mapping * Ce;

    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
        {
            if (this->cell.lattice == LatticeType::Hexagon2D)
            {
                bd->Kn = KnTv(0);
                bd->Tv = KnTv(1);
            }
            else
            {
                bd->Kn = KnTv(bd->layer); // layer is 0 or 1
                bd->Tv = KnTv(2);
            }
        }
    }
}

template <int nlayer>
void ParticleElasticDamage<nlayer>::setParticleProperty(double p_nonlocalL, double p_C11, double p_C12, double p_C44, double p_k0, double p_k1, double p_ct_ratio, double p_d_thres)
{
    k0 = p_k0;
    k1 = p_k1;
    damage_threshold = p_d_thres;
    comp_tensile_ratio = p_ct_ratio;
    this->nonlocal_L = p_nonlocalL;

    Eigen::VectorXd KnTv(NDIM), Ce(NDIM);
    Ce << p_C11, p_C12, p_C44; // C11, C12, C44

    if (this->cell.dim == 2)
        KnTv = this->cell.el_mapping * Ce;
    else
        KnTv = this->cell.radius * this->cell.el_mapping * Ce;

    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
        {
            if (this->cell.lattice == LatticeType::Hexagon2D)
            {
                bd->Kn = KnTv(0);
                bd->Tv = KnTv(1);
            }
            else
            {
                bd->Kn = KnTv(bd->layer); // layer is 0 or 1
                bd->Tv = KnTv(2);
            }
        }
    }
}

#endif
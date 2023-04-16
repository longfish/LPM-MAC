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
// Update bdamage and bforce after calculation
// Three particle-wise state variables: [0]damage, [1] eq_strain, [2]kappa (max eq_strain)

// update bond force
// 1. update the particle geometry (loop for all particles)
// 2. update the particle damage (loop for all particles)
//     i. calculate the eq strain and store it as a state variable - state_var
//     ii. d_eq = state_var - state_var_last
//     iii. using explicit Euler, D = D + g*d_eq
// 3. update bforce (after calculating the bdamage of all bonds)

template <int nlayer>
class ParticleElasticDamage : public Particle<nlayer>
{
public:
    double kappa0{0}, alpha{0}, beta{0}; // damage parameters
    double comp_tensile_ratio{1};
    double damage_threshold{0.99}; // should be less than 1 - e-6

    ParticleElasticDamage(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell) : Particle<nlayer>{p_x, p_y, p_z, p_cell}
    {
        this->state_var = std::vector<double>(3, 0.);
        this->state_var_last = std::vector<double>(3, 0.);
    }

    ParticleElasticDamage(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell, const int &p_type) : Particle<nlayer>{p_x, p_y, p_z, p_cell, p_type}
    {
        this->state_var = std::vector<double>(3, 0.);
        this->state_var_last = std::vector<double>(3, 0.);
    }

    double calcEqStrain();

    void updateParticleStateVar();
    void updateBondsForce();
    void setParticleProperty(double p_E, double p_mu, double p_kappa0, double p_alpha, double p_beta, double p_comp_tensile_ratio);
    void setParticleProperty(double p_C11, double p_C12, double p_C44, double p_kappa0, double p_alpha, double p_beta, double p_comp_tensile_ratio);
};

template <int nlayer>
void ParticleElasticDamage<nlayer>::updateBondsForce()
{
    // calculate the current bond force using current damage value
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
        {
            bd->bforce = 2. * bd->Kn * bd->dLe + 2. * bd->Tv * this->dLe_total[bd->layer]; // trial elastic bforce
            bd->bdamage = std::max(this->state_var[0], bd->p2->state_var[0]);              // update the bond-wise damage
            bd->bforce *= (1.0 - bd->bdamage);
        }
    }
}

double func_eq_strain(const double &comp_tensile_ratio, const double &I1, const double &J2)
{
    double k = comp_tensile_ratio; // the k has different meaning compared with k in func_g
    return (1. - 1. / k) * I1 + 1. / k * std::sqrt((k - 1.) * (k - 1.) * I1 * I1 - k * J2);
}

double func_g(const double k, const double k0, const double a, const double b)
{
    return k0 / k / k * (1. + a * ((1. + b * k) * std::exp(b * (k0 - k)) - 1.)); // the k here is kappa
}

double func_D(const double k, const double k0, const double a, const double b)
{
    return 1.0 - k0 / k * (1.0 - a + a * std::exp(-b * (k - k0)));
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
            J2 += bd->dLe * bd->dLe / bl / bl;
    }
    J2 = 0.5 * I1 * I1 / this->nb - 0.5 * J2;
    return func_eq_strain(comp_tensile_ratio, I1, J2);
}

template <int nlayer>
void ParticleElasticDamage<nlayer>::updateParticleStateVar()
{
    this->state_var[1] = calcEqStrain();                        // local equivalent strain
    double d_eq = this->state_var[1] - this->state_var_last[1]; // incremental equivalent strain
    if (this->state_var[1] - this->state_var[2] > 0 && d_eq > 0 && this->state_var[0] < 1 - EPS)
    {
        // double dp = this->state_var_last[0] + d_eq * func_g(this->state_var_last[2], kappa0, alpha, beta); // update damage using explicit Euler
        this->state_var[2] = std::max(this->state_var[1], this->state_var[2]); // update kappa
        this->state_var[0] = func_D(this->state_var[2], kappa0, alpha, beta);
    }
    else
    {
        this->state_var[0] = this->state_var_last[0];
        this->state_var[2] = this->state_var_last[2];
    }

    if (this->state_var[0] > damage_threshold)
        this->state_var[0] = 1.0;
}

template <int nlayer>
void ParticleElasticDamage<nlayer>::setParticleProperty(double p_E, double p_mu, double p_kappa0, double p_alpha, double p_beta, double p_comp_tensile_ratio)
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
void ParticleElasticDamage<nlayer>::setParticleProperty(double p_C11, double p_C12, double p_C44, double p_kappa0, double p_alpha, double p_beta, double p_comp_tensile_ratio)
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
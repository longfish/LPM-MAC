#pragma once
#ifndef PARTICLE_J2PLASTICITY_H
#define PARTICLE_J2PLASTICITY_H

#include <vector>
#include <array>

#include "lpm.h"
#include "unit_cell.h"
#include "particle.h"
#include "bond.h"

// Elastoplastic isotropic material with a mixed linear hardening law
// State variables (7):
////// alpha, beta[2*NDIM]

template <int nlayer>
class ParticleJ2Plasticity : public Particle<nlayer>
{
public:
    double E{0}, mu{0}, sigmay{0}, xi{0}, H{0}, A{0};
    double damage_threshold{0}, critical_bstrain{0};

    ParticleJ2Plasticity(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell) : Particle<nlayer> { p_x, p_y, p_z, p_cell }
    {
        this->state_var = std::vector<double>(7, 0.);
        this->state_var_last = std::vector<double>(7, 0.);
    }

    ParticleJ2Plasticity(const double &p_x, const double &p_y, const double &p_z, const UnitCell &p_cell, const int &p_type) : Particle<nlayer> { p_x, p_y, p_z, p_cell, p_type }
    {
        this->state_var = std::vector<double>(7, 0.);
        this->state_var_last = std::vector<double>(7, 0.);
    }

    void updateParticleStateVariables();
    bool updateParticleStaticDamage();
    bool updateParticleBrokenBonds();

    void updateBondsForce();
    void setParticleProperty(double p_nonlocalL, bool is_plane_stress, double p_E, double p_mu, double p_sigmay, double p_xi, double p_H, double p_A, double p_critical_bstrain);
};

template <int nlayer>
bool ParticleJ2Plasticity<nlayer>::updateParticleBrokenBonds()
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
void ParticleJ2Plasticity<nlayer>::updateBondsForce()
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
void ParticleJ2Plasticity<nlayer>::updateParticleStateVariables()
{
    double sigma_m = 0.0, sigma_eq = 0.0, triaxiality = 0.0, dlambda = 0.0; // plastic multiplier
    std::vector<double> stress_trial = this->stress;                        // trial stress tensor
    std::vector<double> dplstrain(2 * NDIM, 0.0);                           // delta plastic strain

    /* update stress tensor to be trial devitoric stress tensor */
    sigma_m = 1.0 / 3.0 * (stress_trial[0] + stress_trial[1] + stress_trial[2]);
    for (int j = 0; j < NDIM; j++)
        stress_trial[j] -= sigma_m;

    /* substract trial stress tensor with the back stress tensor, beta */
    for (int j = 0; j < 2 * NDIM; j++)
        stress_trial[j] -= this->state_var[j + 1];

    /* von Mises equivalent stress */
    for (int j = 0; j < 2 * NDIM; j++)
    {
        if (j < NDIM)
            sigma_eq += stress_trial[j] * stress_trial[j]; // s11, s22, s33
        else
            sigma_eq += 2.0 * stress_trial[j] * stress_trial[j]; // s23, s13, s12
    }
    sigma_eq = sqrt(3.0 / 2.0 * sigma_eq);

    /* test the trial yield function */
    double yield_func = sigma_eq - (sigmay + (1.0 - xi) * H * this->state_var[0]);
    double G = E / 2 / (1 + mu); // shear modulus
    if (yield_func > 0.0)
        dlambda = yield_func / (3 * G + H);

    this->state_var[0] += dlambda;

    /* incremental plastic strain tensor */
    for (int j = 0; j < 2 * NDIM; j++)
    {
        if (fabs(sigma_eq) > EPS)
        {
            dplstrain[j] = dlambda * 1.5 * stress_trial[j] / sigma_eq;
            this->state_var[j + 1] += 2. / 3. * xi * H * dplstrain[j];
        }
    }

    /* incremental plastic bond stretch */
    for (int i = 0; i < nlayer; ++i)
    {
        for (Bond<nlayer> *bd : this->bond_layers[i])
        {
            double ddLp = bd->dis * (dplstrain[0] * bd->csx * bd->csx +
                                     dplstrain[1] * bd->csy * bd->csy +
                                     dplstrain[2] * bd->csz * bd->csz +
                                     2 * dplstrain[3] * bd->csy * bd->csz +
                                     2 * dplstrain[4] * bd->csx * bd->csz +
                                     2 * dplstrain[5] * bd->csx * bd->csy);
            bd->dLp += ddLp;
        }
    }

    if (sigma_eq > EPS)
        triaxiality = sigma_m / sigma_eq;
    this->Ddot_local = dlambda * (1.0 + A * triaxiality); // update local damage rate
}

template <int nlayer>
bool ParticleJ2Plasticity<nlayer>::updateParticleStaticDamage()
{
    // update nonlocal damage
    this->damage += this->Ddot_nonlocal;

    if (this->damage >= damage_threshold)
        this->damage = 1;

    return abs(this->damage - this->damage_last) > 1e-3; // return true if damage is updated noticeably
}

template <int nlayer>
void ParticleJ2Plasticity<nlayer>::setParticleProperty(double p_nonlocalL, bool is_plane_stress, double p_E, double p_mu, double p_sigmay, double p_xi, double p_H, double p_A, double p_critical_bstrain)
{
    E = p_E;
    mu = p_mu;
    sigmay = p_sigmay;
    xi = p_xi;
    H = p_H;
    A = p_A;
    critical_bstrain = p_critical_bstrain;
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

#endif
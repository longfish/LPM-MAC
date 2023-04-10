#pragma once
#ifndef BOND_BRITTLE_FATIGUE_WIM_H
#define BOND_BRITTLE_FATIGUE_WIM_H

#include <vector>
#include <array>

#include "lpm.h"
#include "unit_cell.h"
#include "particle.h"
#include "bond.h"

// elastic plane strain or 3D material

template <int nlayer>
class BondBrittleFatigueWim : public Bond<nlayer>
{
public:
    double A{0}, b{0}, c{0}; // fatigue damage constants

    BondBrittleFatigueWim(Particle<nlayer> *p_p1, Particle<nlayer> *p_p2) : Bond<nlayer>{p_p1, p_p2} {}
    BondBrittleFatigueWim(Particle<nlayer> *p_p1, Particle<nlayer> *p_p2, int p_layer, double p_dis) : Bond<nlayer>{p_p1, p_p2, p_layer, p_dis} {}

    bool calcbDamageIndicator();
    void updatebBroken();
    void updatebForce();
    void setBondProperty(double p_E, double p_mu, double p_A, double p_b, double p_c);
};

template <int nlayer>
void BondBrittleFatigueWim<nlayer>::updatebForce()
{
    // for elastic bonds, a trial elastic calculation is enough
    this->bforce_last = this->bforce;
    this->bforce += 2. * this->Kn * this->ddL + 2. * this->Tv * this->p1->ddL_total[this->layer];
    this->bforce *= (1.0 - this->bdamage);
}

template <int nlayer>
bool BondBrittleFatigueWim<nlayer>::calcbDamageIndicator()
{
    this->d_indicator = this->dLe / this->dis_initial; // note the dLe of broken bonds should be already zero
    if (this->d_indicator >= 1 - EPS)
        return true; // potential broken bond
    return false;
}

template <int nlayer>
void BondBrittleFatigueWim<nlayer>::updatebBroken()
{
    if (this->d_indicator >= 1 - EPS)
    {
        this->bdamage = 1.0;
        this->Kn *= (1.0 - this->bdamage);
        this->Tv *= (1.0 - this->bdamage);
        this->damaged = true;
        this->broken = true; // damage is equivalent to broken for brittle materials
    }
}

template <int nlayer>
void BondBrittleFatigueWim<nlayer>::setBondProperty(double p_E, double p_mu, double p_A, double p_b, double p_c)
{
    A = p_A;
    b = p_b;
    c = p_c;

    double Ce[NDIM]{p_E * (1.0 - p_mu) / (1.0 + p_mu) / (1.0 - 2.0 * p_mu),
                    p_E * p_mu / (1.0 + p_mu) / (1.0 - 2.0 * p_mu),
                    p_E / 2.0 / (1.0 + p_mu)}; // C11, C12, C44
    double KnTv[NDIM]{0};

    if (this->p1->cell.dim == 2)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NDIM, 1, NDIM, 1.0, this->p1->cell.el_mapping.data(), 3, Ce, 1, 0.0, KnTv, 1);
    else
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NDIM, 1, NDIM, this->p1->cell.radius, this->p1->cell.el_mapping.data(), 3, Ce, 1, 0.0, KnTv, 1);

    if (this->p1->cell.lattice == LatticeType::Hexagon2D)
    {
        this->Kn = KnTv[0];
        this->Tv = KnTv[1];
    }
    else
    {
        this->Kn = KnTv[this->layer]; // layer is 0 or 1
        this->Tv = KnTv[2];
    }
}

#endif
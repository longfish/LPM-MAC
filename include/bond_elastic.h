#pragma once
#ifndef ELASTIC_BOND_H
#define ELASTIC_BOND_H

#include <vector>
#include <array>

#include "lpm.h"
#include "unit_cell.h"
#include "particle.h"
#include "bond.h"

// elastic plane strain or 3D material

template <int nlayer>
class BondElastic : public Bond<nlayer>
{
public:
    double cr_bstrain{0}; // critical bond strain value at which bond will break

    BondElastic(Particle<nlayer> *p_p1, Particle<nlayer> *p_p2) : Bond<nlayer>{p_p1, p_p2} {}
    BondElastic(Particle<nlayer> *p_p1, Particle<nlayer> *p_p2, int p_layer, double p_dis) : Bond<nlayer>{p_p1, p_p2, p_layer, p_dis} {}

    void updatebDamage();
    void updatebForce();
    void setBondProperty(double p_E, double p_mu, double p_cr_bstrain);
    void setBondProperty(double p_C11, double p_C12, double p_C44, double p_cr_bstrain);
};

template <int nlayer>
void BondElastic<nlayer>::updatebForce()
{
    // for elastic bonds, a trial elastic calculation is enough
    this->bforce_last = this->bforce;
    this->bforce += 2. * this->Kn * this->ddL + 2. * this->Tv * this->p1->ddL_total[this->layer];
    this->bforce *= (1.0 - this->bdamage);
}

template <int nlayer>
void BondElastic<nlayer>::updatebDamage()
{
    this->bstrain = this->dLe / this->dis_initial;
    if (this->bstrain >= cr_bstrain && !(this->damaged))
    {
        this->bdamage = 1.0;
        this->Kn *= (1.0 - this->bdamage);
        this->Tv *= (1.0 - this->bdamage);
        this->damaged = true;
        this->broken = true; // damage is equivalent to broken for brittle materials
    }
}

template <int nlayer>
void BondElastic<nlayer>::setBondProperty(double p_E, double p_mu, double p_cr_bstrain)
{
    cr_bstrain = p_cr_bstrain;

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

template <int nlayer>
void BondElastic<nlayer>::setBondProperty(double p_C11, double p_C12, double p_C44, double p_cr_bstrain)
{
    cr_bstrain = p_cr_bstrain;

    double Ce[NDIM]{p_C11, p_C12, p_C44}; // C11, C12, C44
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
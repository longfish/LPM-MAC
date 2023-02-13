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
class ElasticBond : public Bond<nlayer>
{
public:
    double E{0}, mu{0};   // Youngs modulus and Poisson ratio
    int nbreak{0};        // limit the broken number of bonds in a single iteration, should be an even number
    double cr_bstrain{0}; // critical bond strain value at which bond will break

    ElasticBond(Particle<nlayer> *p_p1, Particle<nlayer> *p_p2) : Bond<nlayer>{p_p1, p_p2} {}
    ElasticBond(Particle<nlayer> *p_p1, Particle<nlayer> *p_p2, int p_layer, double p_dis) : Bond<nlayer>{p_p1, p_p2, p_layer, p_dis} {}

    void updatebForce();
    void setBondProperty(double p_E, double p_mu, double p_cr_bstrain, int p_nbreak);
};

template <int nlayer>
void ElasticBond<nlayer>::updatebForce()
{
    // for elastic bonds, a trial elastic calculation is enough
    this->bforce_last = this->bforce;
    this->bforce += 2. * this->Kn * this->ddL + this->p1->TddL_total[this->layer] + this->Tv * this->p1->ddL_total[this->layer];
}

template <int nlayer>
void ElasticBond<nlayer>::setBondProperty(double p_E, double p_mu, double p_cr_bstrain, int p_nbreak)
{
    E = p_E;
    mu = p_mu;
    nbreak = p_nbreak;
    cr_bstrain = p_cr_bstrain;

    double Ce[NDIM]{E * (1.0 - mu) / (1.0 + mu) / (1.0 - 2.0 * mu),
                    E * mu / (1.0 + mu) / (1.0 - 2.0 * mu),
                    E / 2.0 / (1.0 + mu)}; // C11, C12, C44
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
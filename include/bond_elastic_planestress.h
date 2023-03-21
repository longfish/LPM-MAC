#pragma once
#ifndef ELASTIC_PLANESTRESS_BOND_H
#define ELASTIC_PLANESTRESS_BOND_H

#include <vector>
#include <array>

#include "lpm.h"
#include "unit_cell.h"
#include "particle.h"
#include "bond.h"

// elastic plane stress 2D material

template <int nlayer>
class BondElasticPlaneStress : public Bond<nlayer>
{
public:
    double E{0}, mu{0};   // Youngs modulus and Poisson ratio
    double cr_bstrain{0}; // critical bond strain value at which bond will break

    BondElasticPlaneStress(Particle<nlayer> *p_p1, Particle<nlayer> *p_p2) : Bond<nlayer>{p_p1, p_p2} {}
    BondElasticPlaneStress(Particle<nlayer> *p_p1, Particle<nlayer> *p_p2, int p_layer, double p_dis) : Bond<nlayer>{p_p1, p_p2, p_layer, p_dis} {}

    void setBondProperty(double p_E, double p_mu, double p_cr_bstrain)
    {
        E = p_E;
        mu = p_mu;
        cr_bstrain = p_cr_bstrain;

        double Ce[NDIM]{E / (1.0 - mu) / (1.0 + mu),
                        E * mu / (1.0 - mu) / (1.0 + mu),
                        E / 2.0 / (1.0 + mu)}; // C11, C12, C44
        double KnTv[NDIM]{0};

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NDIM, 1, NDIM, 1.0, this->p1->cell.el_mapping.data(), 3, Ce, 1, 0.0, KnTv, 1);
        if (this->p1->cell.dim != 2)
            printf("Warning: BondElasticPlaneStress bond only work in 2D!\n");

        if (this->p1->cell.lattice == LatticeType::Hexagon2D)
        {
            this->Kn = KnTv[0];
            this->Tv = KnTv[1];
        }
        else
        {
            this->Kn = KnTv[this->layer];
            this->Tv = KnTv[2];
        }
    }
};

#endif
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

    void setBondProperty(double p_E, double p_mu, double p_cr_bstrain, int p_nbreak)
    {
        E = p_E;
        mu = p_mu;
        nbreak = p_nbreak;
        cr_bstrain = p_cr_bstrain;

        double Ce[NDIM]{E * (1.0 - mu) / (1.0 + mu) / (1.0 - 2.0 * mu),
                        E * mu / (1.0 + mu) / (1.0 - 2.0 * mu),
                        E / 2.0 / (1.0 + mu)}; // C11, C12, C44
        double KnTv[NDIM]{0};
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, NDIM, 1, NDIM, 1.0, this->p1->cell.el_mapping, 3, Ce, 1, 0.0, KnTv, 1);

        if (this->p1->cell.lattice == 1)
        {
            Kn = KnTv[0];
            Tv = KnTv[1];
        }
        else
        {
            Kn = KnTv[this->layer];
            Tv = KnTv[2];
        }
    }
};

#endif
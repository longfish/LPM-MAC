#pragma once
#ifndef BOND_H
#define BOND_H

#include "lpm.h"
#include "unit_cell.h"
#include "particle.h"

template <int nlayer>
class Particle;

template <int nlayer>
class Bond
{
public:
    int m_layer;                           // index of current bond layer
    double m_Kn, m_Tv;                     // LPM coefficient
    double m_csx, m_csy, m_csz;            // damage value for visualization
    double m_bforce, m_bstress, m_bstrain; // bond-wise quantities
    Particle<nlayer> m_p1, m_p2;           // particles

    Bond(const Particle<nlayer> &p1, const Particle<nlayer> &p2, int layer)
    {
        m_p1 = p1;
        m_p2 = p2;
        m_layer = layer;
    }
};

#endif
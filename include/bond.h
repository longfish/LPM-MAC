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
    int layer{-1};                   // index of current bond layer
    double distance;                 // bond length
    double Kn, Tv;                   // LPM coefficient
    double csx, csy, csz;            // damage value for visualization
    double bforce, bstress, bstrain; // bond-wise quantities
    Particle<nlayer> p1, p2;         // particles

    Bond(const Particle<nlayer> &p_p1, const Particle<nlayer> &p_p2)
    {
        p1 = p_p1;
        p2 = p_p2;
        distance = p1.distanceTo(p2);
        if ((distance < 1.01 * p1.cell.neighbor1_cutoff) && (p1.id != p2.id))
            layer = 0;
        else if ((distance > 1.01 * p1.cell.neighbor1_cutoff) && (distance < 1.01 * p1.cell.neighbor2_cutoff))
            layer = 1;

        csx = (p1.xyz[0] - p2.xyz[0]) / distance;
        csy = (p1.xyz[1] - p2.xyz[1]) / distance;
        csz = (p1.xyz[2] - p2.xyz[2]) / distance;
    }
};

#endif
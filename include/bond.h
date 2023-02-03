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
    int layer{-1};                              // index of current bond layer
    double dis_initial{0}, dis_last{0}, dis{0}; // initial, last, and current bond length
    double Kn{0}, Tv{0};                        // LPM coefficient
    double csx{0}, csy{0}, csz{0};              // damage value for visualization
    double bforce{0}, bstress{0}, bstrain{0};   // bond-wise quantities
    Particle<nlayer> *p1, *p2;                  // particles are not owned by the bond (only store the location)

    virtual void updatebForce() { bforce = 0; }
    virtual void updatebStress() { bstress = 0; }
    virtual void updatebStrain() { bstrain = 0; }

    Bond(Particle<nlayer> *p_p1, Particle<nlayer> *p_p2)
    {
        p1 = p_p1;
        p2 = p_p2;
        dis = p1->distanceTo(*p2);
        dis_initial = dis;
        if ((dis < 1.01 * p1->cell.neighbor1_cutoff) && (p1->id != p2->id))
            layer = 0;
        else if ((dis > 1.01 * p1->cell.neighbor1_cutoff) && (dis < 1.01 * p1->cell.neighbor2_cutoff))
            layer = 1;

        csx = (p1->xyz[0] - p2->xyz[0]) / dis;
        csy = (p1->xyz[1] - p2->xyz[1]) / dis;
        csz = (p1->xyz[2] - p2->xyz[2]) / dis;
    }

    Bond(Particle<nlayer> *p_p1, Particle<nlayer> *p_p2, int p_layer, double p_dis)
    {
        p1 = p_p1;
        p2 = p_p2;
        dis = p_dis;
        dis_initial = p_dis;
        layer = p_layer;

        csx = (p1->xyz[0] - p2->xyz[0]) / dis;
        csy = (p1->xyz[1] - p2->xyz[1]) / dis;
        csz = (p1->xyz[2] - p2->xyz[2]) / dis;
    }
};

#endif
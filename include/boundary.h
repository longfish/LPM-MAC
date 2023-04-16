#pragma once
#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>

#include "particle.h"

template <int nlayer>
class DispBC
{
public:
    std::vector<Particle<nlayer> *> group;
    char flag;
    double step;

    DispBC(std::vector<Particle<nlayer> *> p_group, char p_flag, double p_step)
        : group{p_group}, flag{p_flag}, step{p_step} {}

    void cutHalf()
    {
        step /= 2.0;
    }
};

template <int nlayer>
class ForceBC
{
public:
    std::vector<Particle<nlayer> *> group;
    double fx;
    double fy;
    double fz;

    ForceBC(std::vector<Particle<nlayer> *> p_group, double p_fx, double p_fy, double p_fz)
        : group{p_group}, fx{p_fx}, fy{p_fy}, fz{p_fz} {}

    void cutHalf()
    {
        fx /= 2.0, fy /= 2.0, fz /= 2.0;
    }
};

#endif
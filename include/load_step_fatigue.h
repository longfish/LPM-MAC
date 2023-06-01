#pragma once
#ifndef LOAD_STEP_FATIGUE_H
#define LOAD_STEP_FATIGUE_H

#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include "load_step.h"
#include "boundary.h"

enum class FatigueLoadType : char
{
    LoadUniaxialDisp
};

// create the template for fatigue loading
template <int nlayer>
class LoadUniaxialDisp : public LoadStep<nlayer>
{
public:
    LoadUniaxialDisp(double x, double y, double z, std::vector<Particle<nlayer> *> &group1, std::vector<Particle<nlayer> *> &group2)
    {
        this->dispBCs.push_back(DispBC<nlayer>(group1, LoadMode::Relative, 'x', x));
        this->dispBCs.push_back(DispBC<nlayer>(group1, LoadMode::Relative, 'y', y));
        this->dispBCs.push_back(DispBC<nlayer>(group1, LoadMode::Relative, 'z', z));
        this->dispBCs.push_back(DispBC<nlayer>(group2, LoadMode::Relative, 'x', 0.0));
        this->dispBCs.push_back(DispBC<nlayer>(group2, LoadMode::Relative, 'y', 0.0));
        this->dispBCs.push_back(DispBC<nlayer>(group2, LoadMode::Relative, 'z', 0.0));
    }
};

#endif
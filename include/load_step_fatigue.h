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
    LoadUniaxialDisp,
    LoadDogBoneDisp,
    LoadCTForce
};

// create the template for fatigue loading
template <int nlayer>
class LoadUniaxialDisp : public LoadStep<nlayer>
{
public:
    LoadUniaxialDisp(double dy,
                    std::vector<Particle<nlayer> *> &left_group,
                    std::vector<Particle<nlayer> *> &top_group,
                    std::vector<Particle<nlayer> *> &bottom_group)
    {
        this->dispBCs.push_back(DispBC<nlayer>(left_group, LoadMode::Relative, 'x', 0.0));
        this->dispBCs.push_back(DispBC<nlayer>(top_group, LoadMode::Relative, 'y', dy));
        this->dispBCs.push_back(DispBC<nlayer>(bottom_group, LoadMode::Relative, 'y', 0.0));
    }
};

template <int nlayer>
class LoadDogBoneDisp : public LoadStep<nlayer>
{
public:
    LoadDogBoneDisp(double x, double y, double z,
                    std::vector<Particle<nlayer> *> &group1,
                    std::vector<Particle<nlayer> *> &group2)
    {
        this->dispBCs.push_back(DispBC<nlayer>(group1, LoadMode::Relative, 'x', x));
        this->dispBCs.push_back(DispBC<nlayer>(group1, LoadMode::Relative, 'y', y));
        this->dispBCs.push_back(DispBC<nlayer>(group1, LoadMode::Relative, 'z', z));
        this->dispBCs.push_back(DispBC<nlayer>(group2, LoadMode::Relative, 'x', 0.0));
        this->dispBCs.push_back(DispBC<nlayer>(group2, LoadMode::Relative, 'y', 0.0));
        this->dispBCs.push_back(DispBC<nlayer>(group2, LoadMode::Relative, 'z', 0.0));
    }
};

template <int nlayer>
class LoadCTForce : public LoadStep<nlayer>
{
public:
    LoadCTForce(double fy,
                std::vector<Particle<nlayer> *> &mid_group,
                std::vector<Particle<nlayer> *> &top_group,
                std::vector<Particle<nlayer> *> &bottom_group)
    {
        this->dispBCs.push_back(DispBC<nlayer>(mid_group, LoadMode::Relative, 'y', 0));
        this->dispBCs.push_back(DispBC<nlayer>(top_group, LoadMode::Relative, 'x', 0.0));
        this->dispBCs.push_back(DispBC<nlayer>(bottom_group, LoadMode::Relative, 'x', 0.0));
        this->forceBCs.push_back(ForceBC<nlayer>(top_group, LoadMode::Relative, 0.0, fy, 0.0));
        this->forceBCs.push_back(ForceBC<nlayer>(bottom_group, LoadMode::Relative, 0.0, -fy, 0.0));
    }
};

#endif
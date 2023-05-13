#pragma once
#ifndef LOAD_STEP_H
#define LOAD_STEP_H

#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include "boundary.h"

template <int nlayer>
class LoadStep
{
public:
    std::vector<DispBC<nlayer>> dispBCs;
    std::vector<ForceBC<nlayer>> forceBCs;

    LoadStep() {}

    void loadCutHalf()
    {
        for (DispBC<nlayer> &d : dispBCs)
            d.cutHalf();

        for (ForceBC<nlayer> &f : forceBCs)
            f.cutHalf();
    }
};

#endif
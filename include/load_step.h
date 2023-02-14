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
    int lind; // load indicator
    std::vector<DispBC<nlayer>> dispBCs;
    std::vector<ForceBC<nlayer>> forceBCs;

    LoadStep(int p_lind)
        : lind{p_lind} {}
};

#endif
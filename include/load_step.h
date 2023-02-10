#pragma once
#ifndef LOAD_STEP_H
#define LOAD_STEP_H

#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include "lpm.h"

class LoadStep
{
public:
    int lind; // load indicator
    std::vector<DispBC> dispBCs;
    std::vector<ForceBC> forceBCs;

    LoadStep(int p_lind)
        : lind{p_lind} {}
};

#endif
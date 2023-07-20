#pragma once
#ifndef SOLUTION_H
#define SOLUTION_H

#include "lpm.h"
#include "assembly.h"

// a class that can manipulate the LPM solution

template <int nlayer>
class Solution
{
public:
    ExactSolution exact_sol_mode;
    Assembly<nlayer> ass;
    std::vector<std::array<double, NDIM>> position;
    std::vector<double> exact_u, approx_u, error;

    Solution(Assembly<nlayer> &p_ass, const ExactSolution &p_exact_sol)
        : ass{p_ass}, exact_sol_mode{p_exact_sol} {}

    double computeL2Error();
    void computeSolBeam001(const double &L, const double &D, const double &P, const double &E, const double &mu);
};

template <int nlayer>
void Solution<nlayer>::computeSolBeam001(const double &L, const double &D, const double &P, const double &E, const double &mu)
{
    // A 2D cantilever beam with one end fixed
    double G = E / 2 / (1 + mu); // shear modulus
    double I = D * D * D / 12;
    for (Particle<nlayer> *pt : ass.pt_sys)
    {
        std::array<double, NDIM> pos = {pt->xyz_initial[0], pt->xyz_initial[1], pt->xyz_initial[2]};

        // Approximate solution
        approx_u.push_back(pt->xyz[0] - pt->xyz_initial[0]);
        approx_u.push_back(pt->xyz[1] - pt->xyz_initial[1]);

        // Exact solution
        double u1 = -P * pow(pos[0], 2) * pos[1] / 2 / E / I + P * pow(pos[1], 3) * (1 / 6 / G / I - mu / 6 / E / I) + pos[1] * (P * L * L / 2 / E / I - P * D * D / 8 / G / I);
        double u2 = mu * P * pos[0] * pow(pos[1], 2) / 2 / E / I + P * pow(pos[0], 3) / 6 / E / I - P * L * L * pos[0] / 2 / E / I + P * L * L * L / 3 / E / I;
        exact_u.push_back(u1);
        exact_u.push_back(u2);
    }
}

template <int nlayer>
double Solution<nlayer>::computeL2Error()
{
    // compute norm square
    double error_sq = 0, exact_sq = 0;
    for (int i = 0; i < exact_u.size(); ++i)
    {
        error_sq += pow(exact_u[i] - approx_u[i], 2);
        exact_sq += pow(exact_u[i], 2);
    }

    if (exact_sq < EPS)
        return -1;
    else
        return sqrt(error_sq / exact_sq);
}
#endif
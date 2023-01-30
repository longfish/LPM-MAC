#include "problem.h"

template <int nlayer>
void LPMProblem<nlayer>::createParticles(std::vector<std::array<double, NDIM>> p_xyz, UnitCell p_cell)
{
    for (auto xyz : p_xyz)
    {
        auto pt = Particle(xyz[0], xyz[1], xyz[2], p_cell);
        psystem.push_back(pt);
    }
}

template <int nlayer>
void LPMProblem<nlayer>::createBonds()
{
    for (Particle<nlayer> p1 : psystem)
    {
        for (Particle<nlayer> p2 : psystem)
        {
            Bond<nlayer> bd { p1, p2 };
            if (bd.layer == 0)
            {
                p1.bond_layers[0].push_back(bd);
                p1.bonds.push_back(bd);
            }
            else if (bd.layer == 1)
            {
                p1.bond_layers[1].push_back(bd);
                p1.bonds.push_back(bd);
            }
        }
    }
}

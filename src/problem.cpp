#include <algorithm>
#include "problem.h"

// template <int nlayer>
// LPMProblem<nlayer>::LPMProblem(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell)
// {
//     createParticles(p_xyz, p_cell);
//     createBonds();
//     createConnections();
//     nparticle = ptsystem.size();
// }

template <int nlayer>
void LPMProblem<nlayer>::createParticles(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell)
{
    for (auto xyz : p_xyz)
    {
        Particle<nlayer> pt{xyz[0], xyz[1], xyz[2], p_cell};
        ptsystem.push_back(pt);
    }
}

template <int nlayer>
void LPMProblem<nlayer>::createBonds()
{
    for (Particle<nlayer> p1 : ptsystem)
    {
        for (Particle<nlayer> p2 : ptsystem)
        {
            Bond<nlayer> bd{p1, p2};
            if (bd.layer == 0)
            {
                p1.bond_layers[0].push_back(bd);
                p1.bonds.push_back(p2);
            }
            else if (bd.layer == 1)
            {
                p1.bond_layers[1].push_back(bd);
                p1.bonds.push_back(p2);
            }
        }
        p1.nb = p1.bonds.size();
    }
}

template <int nlayer>
void LPMProblem<nlayer>::createConnections()
{
    for (Particle<nlayer> p1 : ptsystem)
    {
        for (int i = 0; i < nlayer; i++)
        {
            // loop forward bond particles
            for (Bond<nlayer> bd_fw : p1.bond_layers[i])
            {
                p1.conns.push_back(bd_fw.p2);

                // loop backward bond particles
                for (Bond<nlayer> bd_bw : bd_fw.p2.bond_layers[i])
                {
                    if (std::find(p1.bonds.begin(), p1.bonds.end(), bd_bw.p2) == p1.bonds.end())
                        p1.conns.push_back(bd_bw.p2);
                }
            }
        }
        std::sort(p1.conns.begin(), p1.conns.end(), [](Particle<nlayer> &a, Particle<nlayer> &b)
                  { return a.id < b.id; });
    }
}
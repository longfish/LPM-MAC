#pragma once
#ifndef LPM_SYSTEM_H
#define LPM_SYSTEM_H

#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include "lpm.h"
#include "unit_cell.h"
#include "bond.h"
#include "elastic_bond.h"
#include "elastic_planestress_bond.h"

template <int nlayer>
class Bond;

template <int nlayer>
class ElasticBond;

template <int nlayer>
class Particle;

template <int nlayer>
class Assembly
{
public:
    BondType btype;                         // bond type
    int nparticle;                          // number of particles
    std::vector<Particle<nlayer> *> pt_sys; // system of particles

    Assembly(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell, BondType p_btype); // Construct a particle system

    void createParticles(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell);
    void createBonds();
    void createConnections();

    void updateForceState(); // update bond force and particle forces
    void writeDump(const char *dataName, int step, char flag, double box[]);
};

template <int nlayer>
Assembly<nlayer>::Assembly(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell, BondType p_btype)
{
    btype = p_btype;
    createParticles(p_xyz, p_cell);
    createBonds();
    createConnections();
    nparticle = pt_sys.size();
}

template <int nlayer>
void Assembly<nlayer>::updateForceState()
{
    // update particle geometry and bond force
    for (Particle<nlayer> *pt : pt_sys)
    {
        pt->updateBondsGeometry();
        pt->updateBondsForce();
    }

    for (Particle<nlayer> *pt : pt_sys)
        pt->updateParticleForce();
}

template <int nlayer>
void Assembly<nlayer>::createParticles(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell)
{
    for (auto xyz : p_xyz)
    {
        Particle<nlayer> *pt = new Particle<nlayer>(xyz[0], xyz[1], xyz[2], p_cell);
        pt_sys.push_back(pt);
    }
}

template <int nlayer>
void Assembly<nlayer>::createBonds()
{
#pragma omp parallel for
    for (Particle<nlayer> *p1 : pt_sys)
    {
        for (Particle<nlayer> *p2 : pt_sys)
        {
            double distance = p1->distanceTo(*p2);
            if ((distance < 1.01 * p1->cell.neighbor2_cutoff) && (p1->id != p2->id))
            {
                int layer = 1;
                if (distance < 1.01 * p1->cell.neighbor1_cutoff)
                    layer = 0;

                Bond<nlayer> *bd = nullptr; // create bonds
                if (btype == BondType::Elastic)
                    bd = new ElasticBond<nlayer>(p1, p2, layer, distance);

                p1->bond_layers[layer].push_back(bd);
                p1->neighbors.push_back(p2);
            }
        }
        p1->nb = p1->neighbors.size();
    }
}

template <int nlayer>
void Assembly<nlayer>::createConnections()
{
#pragma omp parallel for
    for (Particle<nlayer> *p1 : pt_sys)
    {
        for (int i = 0; i < nlayer; i++)
        {
            // loop forward bond particles
            for (auto *bd_fw : p1->bond_layers[i])
            {
                p1->conns.push_back(bd_fw->p2);

                // loop backward bond particles
                for (auto *bd_bw : bd_fw->p2->bond_layers[i])
                    p1->conns.push_back(bd_bw->p2);
            }
        }
        std::sort(p1->conns.begin(), p1->conns.end(), [](Particle<nlayer> *a, Particle<nlayer> *b)
                  { return a->id < b->id; });
        p1->conns.erase(unique(p1->conns.begin(), p1->conns.end()), p1->conns.end());
        p1->nconn = p1->conns.size();

        // get the number of conns with id larger or equal to self
        auto curr = std::find(p1->conns.begin(), p1->conns.end(), p1);
        p1->nconn_largeq = (int)std::distance(curr, p1->conns.end());
    }
}

template <int nlayer>
void Assembly<nlayer>::writeDump(const char *dataName, int step, char flag, double box[])
{
    FILE *fpt;
    if (flag == 's')
        fpt = fopen(dataName, "w+");
    else
        fpt = fopen(dataName, "a+");

    fprintf(fpt, "ITEM: TIMESTEP\n");
    fprintf(fpt, "%d\n", step);
    fprintf(fpt, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fpt, "%d\n", nparticle);
    fprintf(fpt, "ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(fpt, "%8.8f %8.8f\n", box[0], box[1]);
    fprintf(fpt, "%8.8f %8.8f\n", box[2], box[3]);
    fprintf(fpt, "%8.8f %8.8f\n", box[4], box[5]);

    fprintf(fpt, "ITEM: ATOMS id type x y z dx dy dz s11 s22 s33 s23 s13 s12 damage\n");
    for (auto pt : pt_sys)
    {
        fprintf(fpt, "%d %d %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",
                pt->id, pt->type,
                pt->xyz[0], pt->xyz[1], pt->xyz[2],
                pt->xyz[0] - pt->xyz_initial[0], pt->xyz[1] - pt->xyz_initial[1], pt->xyz[2] - pt->xyz_initial[2],
                pt->stress[0], pt->stress[1], pt->stress[2], pt->stress[3], pt->stress[4], pt->stress[5],
                pt->damage_visual);
    }

    fclose(fpt);
}

#endif
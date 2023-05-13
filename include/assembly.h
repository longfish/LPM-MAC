#pragma once
#ifndef LPM_SYSTEM_H
#define LPM_SYSTEM_H

#include <vector>
#include <array>
#include <algorithm>
#include <string>
#include <map>

#include "lpm.h"
#include "unit_cell.h"
#include "bond.h"
#include "particle_elastic.h"
#include "particle_elastic_damage.h"
#include "particle_fatigue_hcf.h"

template <int nlayer>
class Bond;

template <int nlayer>
class Particle;

/////////////////////////////////
// Declare particle class
template <int nlayer>
class ParticleElastic;

template <int nlayer>
class ParticleElasticDamage;

template <int nlayer>
class ParticleFatigueHCF;
/////////////////////////////////

template <int nlayer>
class Assembly
{
public:
    ParticleType ptype;                     // particle type
    int nparticle;                          // number of particles
    std::vector<Particle<nlayer> *> pt_sys; // system of particles
    std::array<double, 2 * NDIM> box;       // simulation box

    Assembly(std::vector<std::array<double, NDIM>> &p_xyz, std::array<double, 2 * NDIM> &p_box, UnitCell &p_cell, const ParticleType &p_ptype); // Construct a particle system from scratch
    Assembly(const std::string &dumpFile, UnitCell &p_cell, const ParticleType &p_ptype);                                                       // Assemble the particle system from the dump file
    Assembly(const std::string &dumpFile, const std::string &bondFile, UnitCell &p_cell, const ParticleType &p_ptype);                          // Assemble the particle system from the dump file

    void createParticles(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell);
    void createBonds();
    void updateConnections();
    void updateGeometry();
    void resetStateVar(bool reset_xyz);
    void storeStateVar();
    void clearStateVar(bool clear_damage);
    void updateForceState(); // update bond force and particle forces

    bool updateStateVar();
    bool updateBrokenBonds();

    std::map<int, Particle<nlayer> *> toMap();
    void readBond(const std::string &bondFile);
    void writeBond(const std::string &bondFile);
    void readDump(const std::string &dumpFile, UnitCell &cell);
    void writeDump(const std::string &dumpFile, int step);
    void writeConfigurationDump(const std::string &dumpFile);
};

template <int nlayer>
Assembly<nlayer>::Assembly(const std::string &dumpFile, UnitCell &p_cell, const ParticleType &p_ptype)
{
    ptype = p_ptype;
    printf("Reading data ...\n");
    readDump(dumpFile, p_cell); // extract the particle system
    createBonds();
    updateConnections();
}

template <int nlayer>
Assembly<nlayer>::Assembly(const std::string &dumpFile, const std::string &bondFile, UnitCell &p_cell, const ParticleType &p_ptype)
{
    ptype = p_ptype;
    printf("Reading data ...\n");
    readDump(dumpFile, p_cell);
    readBond(bondFile);

    for (Particle<nlayer> *p1 : pt_sys)
        p1->nb = p1->neighbors.size();

    updateConnections();
}

template <int nlayer>
Assembly<nlayer>::Assembly(std::vector<std::array<double, NDIM>> &p_xyz, std::array<double, 2 * NDIM> &p_box, UnitCell &p_cell, const ParticleType &p_ptype)
{
    ptype = p_ptype;
    box = p_box;
    createParticles(p_xyz, p_cell);
    createBonds();
    updateConnections();
    nparticle = pt_sys.size();
}

template <int nlayer>
std::map<int, Particle<nlayer> *> Assembly<nlayer>::toMap()
{
    std::map<int, Particle<nlayer> *> map;
    for (Particle<nlayer> *pt : pt_sys)
        map[pt->id] = pt;
    return map;
}

template <int nlayer>
void Assembly<nlayer>::updateGeometry()
{
    // update particle geometry
    for (Particle<nlayer> *pt : pt_sys)
        pt->updateBondsGeometry();
}

template <int nlayer>
void Assembly<nlayer>::storeStateVar()
{
    for (Particle<nlayer> *pt : pt_sys)
    {
        pt->xyz_last = pt->xyz;
        pt->storeParticleStateVariables();
    }
}

template <int nlayer>
void Assembly<nlayer>::resetStateVar(bool reset_xyz)
{
    for (Particle<nlayer> *pt : pt_sys)
    {
        if (reset_xyz)
            pt->xyz = pt->xyz_last;
        pt->resetParticleStateVariables();
    }
}

template <int nlayer>
void Assembly<nlayer>::clearStateVar(bool clear_damage)
{
    for (Particle<nlayer> *pt : pt_sys)
    {
        double dam = pt->state_var[0];
        pt->clearParticleStateVariables();
        if (!clear_damage)
            pt->state_var[0] = dam;
    }
}

template <int nlayer>
bool Assembly<nlayer>::updateStateVar()
{
    bool any_damaged{false};
    for (Particle<nlayer> *pt : pt_sys)
        any_damaged = pt->updateParticleStateVariables() || any_damaged;

    return any_damaged;
}

template <int nlayer>
bool Assembly<nlayer>::updateBrokenBonds()
{
    bool any_broken{false};
    for (Particle<nlayer> *pt : pt_sys)
        any_broken = pt->updateParticleBrokenBonds() || any_broken;

    return any_broken;
}

template <int nlayer>
void Assembly<nlayer>::updateForceState()
{
    // this function should be after the geometry updating (ie., above function)
    for (Particle<nlayer> *pt : pt_sys)
        pt->updateBondsForce();

    for (Particle<nlayer> *pt : pt_sys)
    {
        pt->updateParticleForce();
        pt->updateParticleStress();
        pt->updateParticleDamageVisual();
    }
}

template <int nlayer>
void Assembly<nlayer>::createParticles(std::vector<std::array<double, NDIM>> &p_xyz, UnitCell &p_cell)
{
    for (auto xyz : p_xyz)
    {
        Particle<nlayer> *pt = nullptr;
        if (ptype == ParticleType::Elastic)
            pt = new ParticleElastic<nlayer>(xyz[0], xyz[1], xyz[2], p_cell);
        if (ptype == ParticleType::ElasticDamage)
            pt = new ParticleElasticDamage<nlayer>(xyz[0], xyz[1], xyz[2], p_cell);
        if (ptype == ParticleType::FatigueHCF)
            pt = new ParticleFatigueHCF<nlayer>(xyz[0], xyz[1], xyz[2], p_cell);

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
            double distance = p1->distanceTo(p2);
            if ((distance < 1.01 * p1->cell.neighbor_cutoff[1]) && (p1->id != p2->id))
            {
                int layer = 1;
                if (distance < 1.01 * p1->cell.neighbor_cutoff[0])
                    layer = 0;

                Bond<nlayer> *bd = new Bond<nlayer>(p1, p2, layer, distance); // create bonds
                p1->bond_layers[layer].push_back(bd);
                p1->neighbors.push_back(p2);
            }
        }
        p1->nb = p1->neighbors.size();
    }
}

template <int nlayer>
void Assembly<nlayer>::updateConnections()
{
#pragma omp parallel for
    for (Particle<nlayer> *p1 : pt_sys)
    {
        p1->conns.clear(); // clear the connection list
        for (int i = 0; i < nlayer; i++)
        {
            // loop forward bond particles
            for (auto *bd_fw : p1->bond_layers[i])
            {
                p1->conns.push_back(bd_fw->p2);

                // loop backward bond particles
                for (auto *bd_bw : bd_fw->p2->bond_layers[i])
                {
                    p1->conns.push_back(bd_bw->p2);
                }
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
void Assembly<nlayer>::writeBond(const std::string &bondFile)
{
    // please note that the bond has directionality, i.e., bond_12 != bond_21

    FILE *fpt = fopen(bondFile.c_str(), "w+");
    for (Particle<nlayer> *pt : pt_sys)
    {
        for (int i = 0; i < nlayer; ++i)
        {
            for (Bond<nlayer> *bd : pt->bond_layers[i])
                fprintf(fpt, "%d %d %d %d\n", bd->id, bd->layer, bd->p1->id, bd->p2->id);
        }
    }
}

template <int nlayer>
void Assembly<nlayer>::writeConfigurationDump(const std::string &dumpFile)
{
    FILE *fpt = fopen(dumpFile.c_str(), "w+");

    fprintf(fpt, "ITEM: TIMESTEP\n");
    fprintf(fpt, "%d\n", 0);
    fprintf(fpt, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fpt, "%d\n", nparticle);
    fprintf(fpt, "ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(fpt, "%8.8f %8.8f\n", box[0], box[1]);
    fprintf(fpt, "%8.8f %8.8f\n", box[2], box[3]);
    fprintf(fpt, "%8.8f %8.8f\n", box[4], box[5]);

    fprintf(fpt, "ITEM: ATOMS id type x y z\n");
    for (auto pt : pt_sys)
    {
        fprintf(fpt, "%d %d %.4e %.4e %.4e \n",
                pt->id, pt->type,
                pt->xyz[0], pt->xyz[1], pt->xyz[2]);
    }

    fclose(fpt);
}

template <int nlayer>
void Assembly<nlayer>::writeDump(const std::string &dumpFile, int step)
{
    FILE *fpt = fopen(dumpFile.c_str(), "a+");

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

template <int nlayer>
void Assembly<nlayer>::readDump(const std::string &dumpFile, UnitCell &cell)
{
    // Only support dump data type: particle id, type, x, y, z
    // Please note that the particle id needs to start from 0

    FILE *fpt;
    fpt = fopen(dumpFile.c_str(), "r+"); /* read-only */

    if (fpt == NULL)
    {
        printf("\'%s\' does not exist!\n", dumpFile.c_str());
        exit(1);
    }

    int SKIP = 3, nonsense;    // number of lines that need to skip
    const int MAXLENGTH = 500; /* maximum length of a line */
    char *tmp;
    char line[MAXLENGTH]; /* stores lines as they are read in */

    for (int n = 0; n < SKIP; n++) /* skip beginning lines and get total particle number */
        tmp = fgets(line, MAXLENGTH, fpt);

    nonsense = fscanf(fpt, "%d\n", &(nparticle));
    tmp = fgets(line, MAXLENGTH, fpt);

    nonsense = fscanf(fpt, "%lf %lf\n", &(box[0]), &(box[1]));
    nonsense = fscanf(fpt, "%lf %lf\n", &(box[2]), &(box[3]));
    nonsense = fscanf(fpt, "%lf %lf\n", &(box[4]), &(box[5]));
    tmp = fgets(line, MAXLENGTH, fpt);

    /* store into position variable xyz */
    int id, type;
    double xyz[NDIM];
    while (fscanf(fpt, "%d %d %lf %lf %lf", &(id), &(type), &(xyz[0]), &(xyz[1]), &(xyz[2])) > 0)
    {
        Particle<nlayer> *pt = nullptr;
        if (ptype == ParticleType::Elastic)
            pt = new ParticleElastic<nlayer>(xyz[0], xyz[1], xyz[2], cell, type);
        if (ptype == ParticleType::ElasticDamage)
            pt = new ParticleElasticDamage<nlayer>(xyz[0], xyz[1], xyz[2], cell, type);
        if (ptype == ParticleType::FatigueHCF)
            pt = new ParticleFatigueHCF<nlayer>(xyz[0], xyz[1], xyz[2], cell, type);

        pt->id = id;
        pt_sys.push_back(pt);
    }

    fclose(fpt);
}

template <int nlayer>
void Assembly<nlayer>::readBond(const std::string &bondFile)
{
    FILE *fpt;
    fpt = fopen(bondFile.c_str(), "r+"); /* read-only */

    if (fpt == NULL)
    {
        printf("\'%s\' does not exist!\n", bondFile.c_str());
        exit(1);
    }

    std::map<int, Particle<nlayer> *> pt_map = toMap();

    int id, layer, p1id, p2id;
    while (fscanf(fpt, "%d %d %d %d", &(id), &(layer), &(p1id), &(p2id)) == 4)
    {
        // find the two particles
        auto p1 = pt_map[p1id];
        auto p2 = pt_map[p2id];

        // auto p1 = std::find_if(pt_set.begin(), pt_set.end(), [&](const Particle<nlayer> *p)
        //                        { return p->id == p1id; });
        // auto p2 = std::find_if(pt_set.begin(), pt_set.end(), [&](const Particle<nlayer> *p)
        //                        { return p->id == p2id; });

        double distance = p1->distanceTo(p2);
        Bond<nlayer> *bd = new Bond<nlayer>(p1, p2, layer, distance); // create bonds

        p1->bond_layers[layer].push_back(bd);
        p1->neighbors.push_back(p2);
    }

    fclose(fpt);
}

#endif
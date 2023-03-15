#include "lpm.h"
#include "particle.h"
#include "assembly.h"
#include "utilities.h"
#include "stiffness.h"
#include "load_step.h"
#include "solver.h"

/* write the particle-wise information into LAMMPS dump file type */
void writeDump(const std::string &dumpFile, std::vector<std::array<double, NDIM>> &xyz, std::array<double, 2 * NDIM> &box)
{
    FILE *fpt = fopen(dumpFile.c_str(), "w+");

    fprintf(fpt, "ITEM: TIMESTEP\n");
    fprintf(fpt, "%d\n", 0);
    fprintf(fpt, "ITEM: NUMBER OF ATOMS\n");
    fprintf(fpt, "%ld\n", xyz.size());
    fprintf(fpt, "ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(fpt, "%8.8f %8.8f\n", box[0], box[1]);
    fprintf(fpt, "%8.8f %8.8f\n", box[2], box[3]);
    fprintf(fpt, "%8.8f %8.8f\n", box[4], box[5]);

    fprintf(fpt, "ITEM: ATOMS id type x y z\n");
    for (size_t i = 0; i < xyz.size(); ++i)
    {
        fprintf(fpt, "%ld %d %.4e %.4e %.4e \n", i, 0,
                xyz[i][0], xyz[i][1], xyz[i][2]);
    }

    fclose(fpt);
}

/* write the bond-wise information */
void writeBond(const std::string &bondFile, std::vector<std::vector<int>> &bonds)
{
    FILE *fpt = fopen(bondFile.c_str(), "w+");

    for (size_t i = 0; i < bonds.size(); ++i)
        fprintf(fpt, "%ld %d %d %d \n", i, bonds[i][0], bonds[i][1], bonds[i][2]);

    fclose(fpt);
}

/* test if the point is inside the circle (ignore z direction) */
bool inCircle(const std::array<double, NDIM> &pt, const std::array<double, NDIM> &circ_ct, const double r)
{
    double dist_sq = pow(pt[0] - circ_ct[0], 2.0) + pow(pt[1] - circ_ct[1], 2.0);
    return dist_sq <= r * r;
}

/* test if the point is inside the notch (ignore z direction) */
bool inNotch(const std::array<double, NDIM> &pt, const std::array<double, NDIM> &tip, const double hf_width)
{
    // first ensure the point is inside the vertical range
    if (abs(pt[1] - tip[1]) <= hf_width)
    {
        // then the slope should be between [-1,1]
        if ((pt[1] - tip[1]) >= (tip[0] - pt[0]) && (pt[1] - tip[1]) <= (-tip[0] + pt[0]))
            return true;
    }
    return false;
}

/* test if the point is valid or not */
bool isValid(const std::array<double, NDIM> &pt)
{
    if (inCircle(pt, {40.0 - 10.5, 9.2, 0.0}, 9.5 / 2.0))
        return false; // top circle
    if (inCircle(pt, {40.0 - 10.5, 40.0 - 9.2, 0.0}, 9.5 / 2.0))
        return false; // bottom circle
    if (inCircle(pt, {23.0 - 8.3, 20.0 + 8.1, 0.0}, 3.5))
        return false; // random circle (CT-1, refer to Zhang's PD validation paper)
    if (inNotch(pt, {23.0, 20.0, 0.0}, 3.0 / 2.0))
        return false; // notch

    return true;
}

/* search the neighbor for the particle system, without a crack */
std::vector<std::vector<int>> searchNeighbor(std::vector<std::array<double, NDIM>> &xyz, struct UnitCell cell)
{
    int np = xyz.size();
    std::vector<std::vector<int>> neighbors(np, std::vector<int>(cell.nneighbors, -1));
    std::vector<std::vector<int>> neighbors_layer(np, std::vector<int>(cell.nneighbors, -1));

#pragma omp parallel for
    for (int i = 0; i < np; i++)
    {
        int index{0};
        for (int j = 0; j < np; j++)
        {
            double dis = sqrt(pow((xyz[j][0] - xyz[i][0]), 2) + pow((xyz[j][1] - xyz[i][1]), 2) + pow((xyz[j][2] - xyz[i][2]), 2));

            if ((dis < 1.01 * cell.neighbor2_cutoff) && (j != i)) /* The first nearest neighbors */
            {
                // test if the bond is valid
                std::array<double, NDIM> mid{(xyz[j][0] - xyz[i][0]) / 2.0, (xyz[j][1] - xyz[i][1]) / 2.0, 0.0};
                if (isValid(mid))
                {
                    int layer = 1;
                    if (dis < 1.01 * cell.neighbor1_cutoff)
                        layer = 0;
                    neighbors[i][index] = j;
                    neighbors_layer[i][index] = layer;
                    index++;
                }
            }
        }
    }

    std::vector<std::vector<int>> bonds;
    for (int i = 0; i < np; i++)
    {
        int nb = cell.nneighbors - std::count(neighbors[i].begin(), neighbors[i].end(), -1);
        for (int j = 0; j < nb; j++)
        {
            bonds.push_back({neighbors_layer[i][j], i, neighbors[i][j]});
        }
    }

    return bonds;
}

void run()
{
    printf("\nCreating a CT model ...\n");

    const int n_layer = 2; // number of neighbor layers (currently only support 2 layers of neighbors)
    double radius = 0.2;   // particle radius
    UnitCell cell(LatticeType::Square2D, radius);

    // Euler angles setting for system rotation
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    std::array<double, 2 * NDIM> box{0.0, 40.0, 0.0, 40.0, 0.0, 8.0}; // thickness is used for force calculation
    std::vector<std::array<double, NDIM>> sq_xyz = createPlateSQ2D(box, cell, R_matrix), sq_CT;
    for (std::array<double, NDIM> &pt : sq_xyz)
        if (isValid(pt))
            sq_CT.push_back(pt);

    std::cout << "\nTotal particle number is: " << sq_CT.size() << std::endl;
    writeDump("CT_2DSquare.dump", sq_CT, box);

    std::vector<std::vector<int>> sq_bonds = searchNeighbor(sq_CT, cell);
    writeBond("CT_2DSquare.bond", sq_bonds);

    printf("\nDone.\n");
}
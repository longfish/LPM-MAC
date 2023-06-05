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

/* test if the point is inside the rectanglular defined with two diagonal points (ignore z direction) */
bool inRect(const std::array<double, NDIM> &pt, const std::array<double, NDIM> &pt1, const std::array<double, NDIM> &pt2)
{
    return pt[0] >= pt1[0] && pt[1] >= pt1[1] && pt[0] <= pt2[0] && pt[1] <= pt2[1];
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

bool isValid(const std::array<double, NDIM> &pt); // declare the isValid function

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

            if ((dis < 1.01 * cell.neighbor_cutoff[1]) && (j != i)) /* The second nearest neighbors */
            {
                // test if the bond is valid
                // divide the bond into 10 pieces, if any pieces are inside the invalid region, then no bond
                bool valid = true;
                int n_interval = 10;
                double dx = (xyz[j][0] - xyz[i][0]) / n_interval;
                double dy = (xyz[j][1] - xyz[i][1]) / n_interval;
                double dz = (xyz[j][2] - xyz[i][2]) / n_interval;
                for (int k = 0; k < n_interval; ++k)
                {
                    std::array<double, NDIM> pt1{xyz[i][0] + k * dx, xyz[i][1] + k * dy, xyz[i][2] + k * dz};
                    std::array<double, NDIM> pt2{xyz[i][0] + (k + 1) * dx, xyz[i][1] + (k + 1) * dy, xyz[i][2] + (k + 1) * dz};
                    std::array<double, NDIM> mid{(pt1[0] + pt2[0]) / 2.0, (pt1[1] + pt2[1]) / 2.0, (pt1[2] + pt2[2]) / 2.0};
                    valid = isValid(mid) && valid;
                }

                if (valid)
                {
                    int layer = 1;
                    if (dis < 1.01 * cell.neighbor_cutoff[0])
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

/* test if the point is valid or not */
bool isValid(const std::array<double, NDIM> &pt)
{
    if (inCircle(pt, {0.0, 2 * 5.08 + 17.78, 0.0}, 5.08))
        return false; // top left circle
    if (inCircle(pt, {5.08 * 3, 2 * 5.08 + 17.78, 0.0}, 5.08))
        return false; // top right circle
    if (inCircle(pt, {0.0, 2 * 5.08, 0.0}, 5.08))
        return false; // bottom left circle
    if (inCircle(pt, {5.08 * 3, 2 * 5.08, 0.0}, 5.08))
        return false; // bottom right circle
    if (inRect(pt, {0, 2 * 5.08, 0.0}, {5.08, 2 * 5.08 + 17.78, 0.0}))
        return false; // left region
    if (inRect(pt, {5.08 * 2, 2 * 5.08, 0.0}, {5.08 * 3, 2 * 5.08 + 17.78, 0.0}))
        return false; // right region
    // if (inNotch(pt, {63.5 - 22.7, 30.5, 0.0}, 1.0))
    //     return false; // notch for no-hole plate

    return true;
}

void run()
{
    printf("\nCreating a 2D dog bone model ...\n");

    const int n_layer = 2; // number of neighbor layers (currently only support 2 layers of neighbors)
    double radius = 0.16;  // particle radius
    UnitCell cell(LatticeType::Hexagon2D, radius);

    // Euler angles setting for system rotation
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    std::array<double, 2 * NDIM> box{0.0, 5.08 * 3, 0.0, 17.78 + 4 * 5.08, 0.0, 4.8};
    std::vector<std::array<double, NDIM>> sq_xyz = createPlateHEX2D(box, cell, R_matrix), sq_DogBone;

    for (std::array<double, NDIM> &pt : sq_xyz)
        if (isValid(pt))
            sq_DogBone.push_back(pt);

    std::cout << "\nTotal particle number is: " << sq_DogBone.size() << std::endl;
    writeDump("../geometry/geo2_DogBone_2dhex.dump", sq_DogBone, box);

    std::vector<std::vector<int>> sq_bonds = searchNeighbor(sq_DogBone, cell);
    writeBond("../geometry/geo2_DogBone_2dhex.bond", sq_bonds);

    printf("\nDone.\n");
}
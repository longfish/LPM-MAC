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
        // then the slope should be between [-1/sqrt(3), 1/sqrt(3)]
        if ((pt[1] - tip[1]) >= 1. / sqrt(3.) * (tip[0] - pt[0]) && (pt[1] - tip[1]) <= 1. / sqrt(3.) * (-tip[0] + pt[0]))
            return true;
    }
    return false;
}

bool inReverseNotch(const std::array<double, NDIM> &pt, const std::array<double, NDIM> &tip, const double hf_width)
{
    // first ensure the point is inside the vertical range
    if (abs(pt[1] - tip[1]) <= hf_width)
    {
        // then the slope should be between [-1/sqrt(3), 1/sqrt(3)]
        if ((pt[1] - tip[1]) <= 1. / sqrt(3.) * (tip[0] - pt[0]) && (pt[1] - tip[1]) >= 1. / sqrt(3.) * (-tip[0] + pt[0]))
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
                // divide the bond into 10 pieces, if any pieces are inside the invalid region, then the bond itself is invalid
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
    // if (inReverseNotch(pt, {0.5, 0.504, 0.0}, 8e-3 / 2.0))
    //     return false; // notch used for shear loading

    // if (inCircle(pt, {50.8, 16.46, 0.0}, 12.7 / 2))
    //     return false; // bottom circle
    // if (inCircle(pt, {50.8, 60.96 - 16.46, 0.0}, 12.7 / 2))
    //     return false; // top circle
    // // if (inCircle(pt, {18, 35., 0.0}, 4))
    // //     return false; // random circle
    // if (inNotch(pt, {50.8 - 6, 60.96 / 2, 0.0}, 2 / 2.0))
    //     return false; // notch used for calibration

    if (inCircle(pt, {40.0, 14., 0.0}, 5))
        return false; // bottom circle
    if (inCircle(pt, {40.0, 50. - 14, 0.0}, 5))
        return false; // top circle
    if (inCircle(pt, {18, 35., 0.0}, 4))
        return false; // random circle
    if (inNotch(pt, {32.0, 25., 0.0}, 2.08 / 2.0))
        return false; // notch
    // if (inNotch(pt, {63.5 - 22.7, 30.5, 0.0}, 1.0))
    //     return false; // notch for no-hole plate

    return true;
}

void run()
{
    printf("\nCreating a CT model ...\n");

    const int n_layer = 2; // number of neighbor layers (currently only support 2 layers of neighbors)
    double radius = 0.16;  // particle radius
    UnitCell cell(LatticeType::Hexagon2D, radius);

    // Euler angles setting for system rotation
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    // std::array<double, 2 * NDIM> box{0.0, 1, 0.0, 1, 0.0, 8.0}; // thickness is used for force calculation
    // std::array<double, 2 * NDIM> box{0.0, 63.5, 0.0, 60.96, 0.0, 5.0}; // thickness is used for force calculation
    // std::array<double, 2 * NDIM> box{0.0, 50.08, 0.0, 50.8, 0.0, 5.04}; // thickness is used for force calculation
    std::array<double, 2 * NDIM> box{0.0, 50, 0.0, 50, 0.0, 5.0}; // Lu CT specimen
    std::vector<std::array<double, NDIM>> hex_xyz = createPlateHEX2D(box, cell, R_matrix), hex_CT;
    // std::vector<std::array<double, NDIM>> sq_xyz = createPlateSQ2D(box, cell, R_matrix), sq_CT;

    for (std::array<double, NDIM> &pt : hex_xyz)
        if (isValid(pt))
            hex_CT.push_back(pt);

    std::cout << "\nTotal particle number is: " << hex_CT.size() << std::endl;
    writeDump("../geometry/geo1_CT_2DHEX_Lu.dump", hex_CT, box);

    std::vector<std::vector<int>> hex_bonds = searchNeighbor(hex_CT, cell);
    writeBond("../geometry/geo1_CT_2DHEX_Lu.bond", hex_bonds);

    printf("\nDone.\n");
}
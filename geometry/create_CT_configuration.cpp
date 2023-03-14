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

                int layer = 1;
                if (dis < 1.01 * cell.neighbor1_cutoff)
                    layer = 0;
                neighbors[i][index] = j;
                neighbors_layer[i][index] = layer;
                index++;
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

// /* remove particles inside the circle through c (x, y or z) direction */
// void inCircle(double *pc, double ra, char c)
// {
//     double **xyz_t = allocDouble2D(nparticle, 3, 0.);

//     int nparticle_t = 0;
//     for (int i = 0; i < nparticle; i++)
//     {
//         double dis2 = 0.0;
//         if (c == 'x')
//             dis2 = pow(xyz[i][1] - pc[1], 2.0) + pow(xyz[i][2] - pc[2], 2.0);
//         else if (c == 'y')
//             dis2 = pow(xyz[i][0] - pc[0], 2.0) + pow(xyz[i][2] - pc[2], 2.0);
//         else if (c == 'z')
//             dis2 = pow(xyz[i][0] - pc[0], 2.0) + pow(xyz[i][1] - pc[1], 2.0);
//         if (dis2 > ra * ra)
//         {
//             for (int j = 0; j < 3; j++)
//                 xyz_t[nparticle_t][j] = xyz[i][j];

//             nparticle_t++;
//         }
//     }

//     freeDouble2D(xyz, nparticle);
//     xyz = allocDouble2D(nparticle_t, 3, 0.);

//     /* transfer to new position array */
//     for (int i = 0; i < nparticle_t; i++)
//     {
//         for (int j = 0; j < 3; j++)
//             xyz[i][j] = xyz_t[i][j];
//     }

//     freeDouble2D(xyz_t, nparticle);

//     nparticle = nparticle_t;
// }

// /* remove particles inside the circle part through z direction, theta1 < theta2 */
// /* theta1 and theta2 should be both positive or both negative */
// void removeCirclePartZ(double *pc, double ra, double theta1, double theta2)
// {
//     double **xyz_t = allocDouble2D(nparticle, 3, 0.);

//     int nparticle_t = 0;
//     for (int i = 0; i < nparticle; i++)
//     {
//         double dx = 0, dy = 0;
//         dx = xyz[i][0] - pc[0];
//         dy = xyz[i][1] - pc[1];
//         if (fabs(dx) < EPS && fabs(dy) < EPS)
//             continue;
//         else
//         {
//             double theta = atan2(dy, dx);
//             double dis2 = pow(xyz[i][0] - pc[0], 2.0) + pow(xyz[i][1] - pc[1], 2.0);
//             if (dis2 > ra * ra || (theta >= theta1 && theta <= theta2))
//             {
//                 for (int j = 0; j < 3; j++)
//                     xyz_t[nparticle_t][j] = xyz[i][j];

//                 nparticle_t++;
//             }
//         }
//     }

//     freeDouble2D(xyz, nparticle);
//     xyz = allocDouble2D(nparticle_t, 3, 0.);

//     /* transfer to new position array */
//     for (int i = 0; i < nparticle_t; i++)
//     {
//         for (int j = 0; j < 3; j++)
//             xyz[i][j] = xyz_t[i][j];
//     }

//     freeDouble2D(xyz_t, nparticle);

//     nparticle = nparticle_t;
// }

// // remove a block with size: xlo=r0, xhi=r1, ylo=r2, yhi=r3, zlo=r4, zhi=r5
// void removeBlock(double r0, double r1, double r2, double r3, double r4, double r5)
// {
//     double **xyz_t = allocDouble2D(nparticle, 3, 0.);

//     int nparticle_t = 0;
//     for (int i = 0; i < nparticle; i++)
//     {
//         if (xyz[i][0] < r0 || xyz[i][0] > r1 || xyz[i][1] < r2 || xyz[i][1] > r3 || xyz[i][2] < r4 || xyz[i][2] > r5)
//         {
//             for (int j = 0; j < 3; j++)
//                 xyz_t[nparticle_t][j] = xyz[i][j];

//             nparticle_t++;
//         }
//     }

//     freeDouble2D(xyz, nparticle);
//     xyz = allocDouble2D(nparticle_t, 3, 0.);

//     /* transfer to new position array */
//     for (int i = 0; i < nparticle_t; i++)
//     {
//         for (int j = 0; j < 3; j++)
//             xyz[i][j] = xyz_t[i][j];
//     }

//     freeDouble2D(xyz_t, nparticle);

//     nparticle = nparticle_t;
// }

// // create a cylinder from the original cuboid configuration
// void createCylinderz(double *pc, double ra)
// {
//     double **xyz_t = allocDouble2D(nparticle, 3, 0.);

//     int nparticle_t = 0;
//     for (int i = 0; i < nparticle; i++)
//     {
//         if (pow(xyz[i][0] - pc[0], 2.0) + pow(xyz[i][1] - pc[1], 2.0) < ra * ra)
//         {
//             for (int j = 0; j < 3; j++)
//                 xyz_t[nparticle_t][j] = xyz[i][j];

//             nparticle_t++;
//         }
//     }

//     freeDouble2D(xyz, nparticle);
//     xyz = allocDouble2D(nparticle_t, 3, 0.);

//     /* transfer to new position array */
//     for (int i = 0; i < nparticle_t; i++)
//     {
//         for (int j = 0; j < 3; j++)
//             xyz[i][j] = xyz_t[i][j];
//     }

//     freeDouble2D(xyz_t, nparticle);

//     nparticle = nparticle_t;
// }

// // remove a ring from the initial configuration, normal direction of the ring is along z direction
// void removeRingz(double *pc, double R, double r)
// {
//     double **xyz_t = allocDouble2D(nparticle, 3, 0.);

//     int nparticle_t = 0;
//     for (int i = 0; i < nparticle; i++)
//     {
//         double h = xyz[i][2] - pc[2];
//         if (fabs(h) < r)
//         {
//             if (pow(xyz[i][0] - pc[0], 2.0) + pow(xyz[i][1] - pc[1], 2.0) < pow(R - sqrt(r * r - h * h), 2.0))
//             {
//                 for (int j = 0; j < 3; j++)
//                     xyz_t[nparticle_t][j] = xyz[i][j];

//                 nparticle_t++;
//             }
//         }
//         else
//         {
//             for (int j = 0; j < 3; j++)
//                 xyz_t[nparticle_t][j] = xyz[i][j];

//             nparticle_t++;
//         }
//     }

//     freeDouble2D(xyz, nparticle);
//     xyz = allocDouble2D(nparticle_t, 3, 0.);

//     /* transfer to new position array */
//     for (int i = 0; i < nparticle_t; i++)
//     {
//         for (int j = 0; j < 3; j++)
//             xyz[i][j] = xyz_t[i][j];
//     }

//     freeDouble2D(xyz_t, nparticle);

//     nparticle = nparticle_t;
// }

// /* Create an initial crack which starts at edge, with length a, width 2w and height h, (perpendicular to xOy plane) */
// void createCrack(double a1, double a2, double w, double h)
// {
//     double **xyz_t = allocDouble2D(nparticle, 3, 0.);

//     int nparticle_t = 0;
//     for (int i = 0; i < nparticle; i++)
//     {
//         if (xyz[i][0] < a1 || xyz[i][0] > a2 || xyz[i][1] > h + w || xyz[i][1] < h - w)
//         {
//             // for (int j = 0; j < 3; j++)
//             //     xyz_t[nparticle_t][j] = xyz[i][j];

//             // nparticle_t++;

//             // with cycle
//             double dis2 = pow(xyz[i][0] - a2, 2.0) + pow(xyz[i][1] - h, 2.0);
//             if (dis2 > w * w)
//             {
//                 for (int j = 0; j < 3; j++)
//                     xyz_t[nparticle_t][j] = xyz[i][j];

//                 nparticle_t++;
//             }
//         }
//     }
//     freeDouble2D(xyz, nparticle);
//     xyz = allocDouble2D(nparticle_t, 3, 0.);

//     /* transfer to new position array */
//     for (int i = 0; i < nparticle_t; i++)
//     {
//         for (int j = 0; j < 3; j++)
//             xyz[i][j] = xyz_t[i][j];
//     }
//     freeDouble2D(xyz_t, nparticle);
//     nparticle = nparticle_t;
// }

// /* define the initial crack, length=|a1-a2|, height=h */
// void defineCrack(double a1, double a2, double h)
// {
//     for (int i = 0; i < nparticle; i++)
//     {
//         damage_visual[i] = 0.0;
//         for (int j = 0; j < nb_initial[i]; j++)
//         {
//             if ((xyz[neighbors[i][j]][1] - h) * (xyz[i][1] - h) < 0 && xyz[i][0] > a1 && xyz[i][0] < a2 &&
//                 xyz[neighbors[i][j]][0] > a1 && xyz[neighbors[i][j]][0] < a2)
//             {
//                 damage_D[i][j][0] = 1.0;
//                 damage_D[i][j][1] = 1.0;
//                 damage_broken[i][j] = 0.0;
//                 damage_w[i][j] = 0.0;
//                 nb[i] -= 1;
//             }
//             damage_visual[i] += damage_broken[i][j];
//         }
//         damage_visual[i] = 1 - damage_visual[i] / nb_initial[i];
//     }
// }

void run()
{
    printf("\nCreating a CT model ...\n");

    const int n_layer = 2; // number of neighbor layers (currently only support 2 layers of neighbors)
    double radius = 0.25;  // particle radius
    UnitCell cell(LatticeType::Square2D, radius);

    // Euler angles setting for system rotation
    int eulerflag = 0; // direct rotation
    double angles[] = {PI / 180.0 * 0.0, PI / 180.0 * 0.0, PI / 180.0 * 0.0};
    double *R_matrix = createRMatrix(eulerflag, angles);

    std::array<double, 2 * NDIM> box{0.0, 40.0, 0.0, 40.0, 0.0, 8.0}; // thickness is used for force calculation
    std::vector<std::array<double, NDIM>> sq_xyz = createPlateSQ2D(box, cell, R_matrix);
    writeDump("CT_2DSquare.dump", sq_xyz, box);

    std::vector<std::vector<int>> sq_bonds = searchNeighbor(sq_xyz, cell);
    writeBond("CT_2DSquare.bond", sq_bonds);

    printf("\nDone.\n");
}
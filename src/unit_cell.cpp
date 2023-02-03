#include <math.h>
#include "unit_cell.h"

UnitCell::UnitCell(int p_lattice, double p_radius)
{
    lattice = p_lattice;
    radius = p_radius;

    if (lattice == 0) /* 2D square lattice (double layer neighbor) */
    {
        dim = 2;
        nneighbors = 8;
        nneighbors1 = 4;
        nneighbors2 = 4;
        neighbor1_cutoff = 2.0 * radius;
        neighbor2_cutoff = 2.0 * sqrt(2.0) * radius;
        nneighbors_AFEM = 16;
        nneighbors_AFEM1 = 12;
        nneighbors_AFEM2 = 12;
        particle_volume = 4 * pow(radius, 2); /* volume (area) of unit cell */
        el_mapping[0] = 1 / 2.0, el_mapping[1] = -1 / 2.0, el_mapping[2] = 0.0;
        el_mapping[3] = 0.0, el_mapping[4] = 0.0, el_mapping[5] = 1. / 2;
        el_mapping[6] = 0.0, el_mapping[7] = 1 / 12.0, el_mapping[8] = -1 / 12.0;
    }

    if (lattice == 1) /* 2D hexagon lattice (double layer neighbor) */
    {
        dim = 2;
        nneighbors = 12;
        nneighbors1 = 6;
        nneighbors2 = 6;
        nneighbors_AFEM = 30;
        nneighbors_AFEM1 = 18;
        nneighbors_AFEM2 = 18;
        neighbor1_cutoff = 2.0 * radius;
        neighbor2_cutoff = 2.0 * sqrt(3.0) * radius;
        particle_volume = 2 * sqrt(3) * pow(radius, 2); /* volume (area) of unit cell */
        el_mapping[0] = sqrt(3.0) / 12.0, el_mapping[1] = -sqrt(3.0) / 12.0, el_mapping[2] = 0.0;
        el_mapping[3] = -sqrt(3.0) / 144.0, el_mapping[4] = sqrt(3.0) / 48.0, el_mapping[5] = 0.0;
        el_mapping[6] = 0.0, el_mapping[7] = 0.0, el_mapping[8] = 1.0;
    }

    if (lattice == 2) /* simple cubic */
    {
        dim = 3;
        nneighbors = 18;
        nneighbors1 = 6;
        nneighbors2 = 12;
        neighbor1_cutoff = 2.0 * radius;
        neighbor2_cutoff = 2.0 * sqrt(2.0) * radius;
        nneighbors_AFEM = 60;
        nneighbors_AFEM1 = 24;
        nneighbors_AFEM2 = 54;
        particle_volume = pow(2 * radius, 3); /* volume of unit cell */
        el_mapping[0] = 1, el_mapping[1] = -1, el_mapping[2] = -1.0;
        el_mapping[3] = 0.0, el_mapping[4] = 0.0, el_mapping[5] = 1.;
        el_mapping[6] = 0.0, el_mapping[7] = 1 / 18.0, el_mapping[8] = -1 / 18.0;
    }

    if (lattice == 3) /* face centered cubic  */
    {
        dim = 3;
        nneighbors = 18;
        nneighbors1 = 12;
        nneighbors2 = 6;
        neighbor1_cutoff = 2.0 * radius;
        neighbor2_cutoff = 2.0 * sqrt(2.0) * radius;
        nneighbors_AFEM = 60;
        nneighbors_AFEM1 = 54;
        nneighbors_AFEM2 = 24;
        particle_volume = 4.0 * sqrt(2.0) * pow(radius, 3); /* volume of unit cell 1, rhombic dodecahedron */
        el_mapping[0] = 0.0, el_mapping[1] = 0.0, el_mapping[2] = sqrt(2.0);
        el_mapping[3] = sqrt(2.0) / 4, el_mapping[4] = -sqrt(2.0) / 4, el_mapping[5] = -sqrt(2.0) / 4;
        el_mapping[6] = 0.0, el_mapping[7] = sqrt(2.0) / 24, el_mapping[8] = -sqrt(2.0) / 24;
    }

    if (lattice == 4) /* body centered cubic  */
    {
        dim = 3;
        nneighbors = 14;
        nneighbors1 = 8;
        nneighbors2 = 6;
        neighbor1_cutoff = 2.0 * radius;
        neighbor2_cutoff = 4.0 / sqrt(3.0) * radius;
        nneighbors_AFEM = 40;
        nneighbors_AFEM1 = 34;
        nneighbors_AFEM2 = 24;
        particle_volume = 32.0 * sqrt(3.0) / 9.0 * pow(radius, 3); /* volume of unit cell 1, truncated octahedron */
        el_mapping[0] = 0.0, el_mapping[1] = 0.0, el_mapping[2] = sqrt(3.0);
        el_mapping[3] = 1/sqrt(3.0) , el_mapping[4] = -1/sqrt(3.0) , el_mapping[5] = 0.;
        el_mapping[6] = 0.0, el_mapping[7] = sqrt(3.0) / 14, el_mapping[8] = sqrt(3.0) / 14;
    }
}

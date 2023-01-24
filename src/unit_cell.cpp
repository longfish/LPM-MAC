#include <math.h>
#include "unit_cell.h"

UnitCell::UnitCell(int lattice, double radius)
{
    UnitCell::lattice = lattice;
    UnitCell::radius = radius;

    if (lattice == 0) /* 2D square lattice (double layer neighbor) */
    {
        UnitCell::dim = 2;
        UnitCell::nneighbors = 8;
        UnitCell::nneighbors1 = 4;
        UnitCell::nneighbors2 = 4;
        UnitCell::neighbor1_cutoff = 2.0 * radius;
        UnitCell::neighbor2_cutoff = 2.0 * sqrt(2.0) * radius;
        UnitCell::nneighbors_AFEM = 16;
        UnitCell::nneighbors_AFEM1 = 12;
        UnitCell::nneighbors_AFEM2 = 12;
        UnitCell::particle_volume = 4 * pow(radius, 2); /* volume (area) of unit cell */
    }

    if (lattice == 1) /* 2D hexagon lattice (double layer neighbor) */
    {
        UnitCell::dim = 2;
        UnitCell::nneighbors = 12;
        UnitCell::nneighbors1 = 6;
        UnitCell::nneighbors2 = 6;
        UnitCell::nneighbors_AFEM = 30;
        UnitCell::nneighbors_AFEM1 = 18;
        UnitCell::nneighbors_AFEM2 = 18;
        UnitCell::neighbor1_cutoff = 2.0 * radius;
        UnitCell::neighbor2_cutoff = 2.0 * sqrt(3.0) * radius;
        UnitCell::particle_volume = 2 * sqrt(3) * pow(radius, 2); /* volume (area) of unit cell */
    }

    if (lattice == 2) /* simple cubic */
    {
        UnitCell::dim = 3;
        UnitCell::nneighbors = 18;
        UnitCell::nneighbors1 = 6;
        UnitCell::nneighbors2 = 12;
        UnitCell::neighbor1_cutoff = 2.0 * radius;
        UnitCell::neighbor2_cutoff = 2.0 * sqrt(2.0) * radius;
        UnitCell::nneighbors_AFEM = 60;
        UnitCell::nneighbors_AFEM1 = 24;
        UnitCell::nneighbors_AFEM2 = 54;
        UnitCell::particle_volume = pow(2 * radius, 3); /* volume of unit cell */
    }

    if (lattice == 3) /* face centered cubic  */
    {
        UnitCell::dim = 3;
        UnitCell::nneighbors = 18;
        UnitCell::nneighbors1 = 12;
        UnitCell::nneighbors2 = 6;
        UnitCell::neighbor1_cutoff = 2.0 * radius;
        UnitCell::neighbor2_cutoff = 2.0 * sqrt(2.0) * radius;
        UnitCell::nneighbors_AFEM = 60;
        UnitCell::nneighbors_AFEM1 = 54;
        UnitCell::nneighbors_AFEM2 = 24;
        UnitCell::particle_volume = 4.0 * sqrt(2.0) * pow(radius, 3); /* volume of unit cell 1, rhombic dodecahedron */
    }

    if (lattice == 4) /* body centered cubic  */
    {
        UnitCell::dim = 3;
        UnitCell::nneighbors = 14;
        UnitCell::nneighbors1 = 8;
        UnitCell::nneighbors2 = 6;
        UnitCell::neighbor1_cutoff = 2.0 * radius;
        UnitCell::neighbor2_cutoff = 4.0 / sqrt(3.0) * radius;
        UnitCell::nneighbors_AFEM = 40;
        UnitCell::nneighbors_AFEM1 = 34;
        UnitCell::nneighbors_AFEM2 = 24;
        UnitCell::particle_volume = 32.0 * sqrt(3.0) / 9.0 * pow(radius, 3); /* volume of unit cell 1, truncated octahedron */
    }
}

#include "particle.h"

template <int nlayer>
Particle<nlayer>::Particle(double p_x, double p_y, double p_z, int p_lattice, double p_radius)
{
    xyz[0] = p_x;
    xyz[1] = p_y;
    xyz[2] = p_z;
    radius = p_radius;
    lattice = p_lattice;
    cell = UnitCell(p_lattice, p_radius);
    id = _ID++;
}

template <int nlayer>
Particle<nlayer>::Particle(double p_x, double p_y, double p_z, UnitCell p_cell)
{
    xyz[0] = p_x;
    xyz[1] = p_y;
    xyz[2] = p_z;
    cell = p_cell;
    id = _ID++;
}

template <int nlayer>
Particle<nlayer>::Particle(const Particle<nlayer> &A)
{
    xyz[0] = A.xyz[0];
    xyz[1] = A.xyz[1];
    xyz[2] = A.xyz[2];
    radius = A.radius;
    lattice = A.lattice;
    cell = A.cell;
    id = A.id;
}

template <int nlayer>
Particle<nlayer> &Particle<nlayer>::operator=(const Particle<nlayer> &A)
{
    xyz[0] = A.xyz[0];
    xyz[1] = A.xyz[1];
    xyz[2] = A.xyz[2];
    radius = A.radius;
    lattice = A.lattice;
    cell = A.cell;
    id = A.id;
    return (*this);
}

template <int nlayer>
void Particle<nlayer>::moveTo(double new_x, double new_y, double new_z)
{
    xyz[0] = new_x;
    xyz[1] = new_y;
    xyz[2] = new_z;
}

template <int nlayer>
double Particle<nlayer>::distanceTo(const Particle<nlayer> &A)
{
    return sqrt(pow((xyz[0] - A.xyz[0]), 2) + pow((xyz[1] - A.xyz[1]), 2) + pow((xyz[2] - A.xyz[2]), 2));
}
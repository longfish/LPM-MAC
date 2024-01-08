#pragma once
#ifndef UNIT_CELL_H
#define UNIT_CELL_H

#include "lpm.h"

class UnitCell
{
public:
    LatticeType lattice;
    int dim;
    int nneighbors, nneighbors1, nneighbors2;
    int nneighbors_AFEM, nneighbors_AFEM1, nneighbors_AFEM2;
    double radius;
    double particle_volume;     /* volume of unit cell */
    Eigen::MatrixXd el_mapping; /* mapping from elasticity to Kn, Tv */
    std::vector<double> neighbor_cutoff;

    UnitCell(LatticeType p_lattice, double p_radius);
    void setSquare();
    void setHexagon();
    void setSimpleCubic();
    void setFCC();
    void setBCC();
};

UnitCell::UnitCell(LatticeType p_lattice, double p_radius)
{
    lattice = p_lattice;
    radius = p_radius;

    if (lattice == LatticeType::Square2D) /* 2D square lattice (double layer neighbor) */
        setSquare();
    else if (lattice == LatticeType::Hexagon2D) /* 2D hexagon lattice (double layer neighbor) */
        setHexagon();
    else if (lattice == LatticeType::SimpleCubic3D) /* simple cubic */
        setSimpleCubic();
    else if (lattice == LatticeType::FCC3D) /* face centered cubic  */
        setFCC();
    else if (lattice == LatticeType::BCC3D_1Slip || lattice == LatticeType::BCC3D_2Slip || lattice == LatticeType::BCC3D_3Slip) /* body centered cubic  */
        setBCC();
}

void UnitCell::setSquare()
{
    double data[] = {1. / 2., -1. / 2., 0.,
                     0., 0., 1. / 2.,
                     0., 1. / 12., -1. / 12.};
    el_mapping = Eigen::Map<const Eigen::MatrixXd>(data, NDIM, NDIM);
    dim = 2;
    nneighbors = 8;
    nneighbors1 = 4;
    nneighbors2 = 4;
    neighbor_cutoff = {2.0 * radius, 2.0 * sqrt(2.0) * radius};
    nneighbors_AFEM = 16;
    nneighbors_AFEM1 = 12;
    nneighbors_AFEM2 = 12;
    particle_volume = 4 * pow(radius, 2); /* volume (area) of unit cell */
}

void UnitCell::setHexagon()
{
    double data[] = {sqrt(3.0) / 12.0, -sqrt(3.0) / 12.0, 0.,
                     -sqrt(3.0) / 144.0, sqrt(3.0) / 48.0, 0.0,
                     0.0, 0.0, 1.0};
    el_mapping = Eigen::Map<const Eigen::MatrixXd>(data, NDIM, NDIM);
    dim = 2;
    nneighbors = 12;
    nneighbors1 = 6;
    nneighbors2 = 6;
    nneighbors_AFEM = 30;
    nneighbors_AFEM1 = 18;
    nneighbors_AFEM2 = 18;
    neighbor_cutoff = {2.0 * radius, 2.0 * sqrt(3.0) * radius};
    particle_volume = 2 * sqrt(3) * pow(radius, 2); /* volume (area) of unit cell */
}

void UnitCell::setSimpleCubic()
{
    double data[] = {1, -1, -1.0,
                     0.0, 0.0, 1.,
                     0.0, 1 / 18.0, -1 / 18.0};
    el_mapping = Eigen::Map<const Eigen::MatrixXd>(data, NDIM, NDIM);
    dim = 3;
    nneighbors = 18;
    nneighbors1 = 6;
    nneighbors2 = 12;
    neighbor_cutoff = {2.0 * radius, 2.0 * sqrt(2.0) * radius};
    nneighbors_AFEM = 60;
    nneighbors_AFEM1 = 24;
    nneighbors_AFEM2 = 54;
    particle_volume = pow(2 * radius, 3); /* volume of unit cell */
}

void UnitCell::setFCC()
{
    double data[] = {0.0, 0.0, sqrt(2.0),
                     sqrt(2.0) / 4, -sqrt(2.0) / 4, -sqrt(2.0) / 4,
                     0.0, sqrt(2.0) / 24, -sqrt(2.0) / 24};
    el_mapping = Eigen::Map<const Eigen::MatrixXd>(data, NDIM, NDIM);
    dim = 3;
    nneighbors = 18;
    nneighbors1 = 12;
    nneighbors2 = 6;
    neighbor_cutoff = {2.0 * radius, 2.0 * sqrt(2.0) * radius};
    nneighbors_AFEM = 60;
    nneighbors_AFEM1 = 54;
    nneighbors_AFEM2 = 24;
    particle_volume = 4.0 * sqrt(2.0) * pow(radius, 3); /* volume of unit cell 1, rhombic dodecahedron */
}

void UnitCell::setBCC()
{
    double data[] = {0.0, 0.0, sqrt(3.0),
                     1 / sqrt(3.0), -1 / sqrt(3.0), 0.,
                     0.0, sqrt(3.0) / 14, sqrt(3.0) / 14};
    el_mapping = Eigen::Map<const Eigen::MatrixXd>(data, NDIM, NDIM);
    dim = 3;
    nneighbors = 14;
    nneighbors1 = 8;
    nneighbors2 = 6;
    neighbor_cutoff = {2.0 * radius, 4.0 / sqrt(3.0) * radius};
    nneighbors_AFEM = 40;
    nneighbors_AFEM1 = 34;
    nneighbors_AFEM2 = 24;
    particle_volume = 32.0 * sqrt(3.0) / 9.0 * pow(radius, 3); /* volume of unit cell 1, truncated octahedron */
}

#endif
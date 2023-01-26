#pragma once
#ifndef INITIALIZATION_H
#define INITIALIZATION_H

#include "lpm.h"
#include "unit_cell.h"

void createCuboid(double box[], UnitCell cell, double R_matrix[]);
double * createRMatrix(int eulerflag, double angles[]);
void createCylinderz(double *pc, double ra);
void moveParticle(double box[], double *movexyz);
void removeCircle(double *pc, double ra, char c);
void removeCirclePartZ(double *pc, double ra, double theta1, double theta2);
void removeBlock(double r0, double r1, double r2, double r3, double r4, double r5);
void removeRingz(double *pc, double R, double r);
void createCrack(double a1, double a2, double w, double h);
void defineCrack(double a1, double a2, double h);

#endif

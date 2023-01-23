#pragma once
#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include "lpm.h"

struct UnitCell;

// basic lpm computations
void setTypeRect(double r0, double r1, double r2, double r3, double r4, double r5, int t);
void setTypeCircle(double x, double y, double r, int t);
void setTypeFullNeighbor(int k, struct UnitCell cell);
void computedL();
void computeStress(struct UnitCell cell);
void computeStrain(struct UnitCell cell);
void computeBondForce(struct UnitCell cell);
void computeBondStretch(struct UnitCell cell);

// basic operations
int findMaxInt(int *arr, int len);
int findMaxDouble(double *arr, int len);
int getLine(char *line, int max, FILE *fp);
int sumInt(int *arr, int len);
int countEqual(int *arr, int len, int b);
int countNEqual(int *arr, int len, int b);
int countLarger(int *arr, int len, int b);
int findNElement(int *in, int *out, int len, int a);
double sumDouble(double *arr, int len);
void sortInt(int *arr, int len);
void printDouble(double *target, int m, int n);
void printInt(int *target, int n);

/* allocate memory and initialize specific value into the array */
int *allocInt1D(int num, int a);
int **allocInt2D(int row, int column, int a);
int ***allocInt3D(int row, int column, int layer, int a);
double *allocDouble1D(int num, double a);
double **allocDouble2D(int row, int column, double a);
double ***allocDouble3D(int row, int column, int layer, double a);

void freeInt2D(int **arr, int row);
void freeDouble2D(double **arr, int row);
void freeDouble3D(double ***arr, int row, int column);

void copyInt2D(int **target, int **source, int row, int col);
void copyDouble1D(double *target, double *source, int row);
void copyDouble2D(double **target, double **source, int row, int col);
void copyDouble3D(double ***target, double ***source, int row, int col, int deep);

#endif

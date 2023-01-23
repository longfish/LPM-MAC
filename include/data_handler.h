#pragma once
#ifndef DATA_HANDLER_H
#define DATA_HANDLER_H

#include "lpm.h"
struct UnitCell;

//void readLammps(const char *dataName, int skip);
void writeDump(const char *dataName, int step, char flag, double box[], double plmode);
void writeRSS(const char *dataName, int step);
void writeK_global(const char *dataName, int l);
void writeDlambda(const char *dataName, int m, int n, int globalStep, int iterStep);
void writeDisp(const char *dataName, char c, int t1, int tstep);
void writeForce(const char *dataname, char c, double p, int tstep, double box[]);
void writeReaction(const char *dataName, char c, int t1, int tstep);
void writeStress(const char *dataName, int t1, int tstep);
void writeStrain(const char *dataName, int t1, int tstep);
void writeBondstretch(const char *dataName, int step);
void writeBondforce(const char *dataName, int step, struct UnitCell cell);
void writeJact(const char *dataName, int step);
void writeInternalForce(const char *dataName, int step);
void writeDamage(const char *dataName, int step);
void writeNeighbor(const char *dataName);
void writeCab(const char *dataName, int ii);

#endif

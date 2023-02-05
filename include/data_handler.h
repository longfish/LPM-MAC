#pragma once
#ifndef DATA_HANDLER_H
#define DATA_HANDLER_H

#include "lpm.h"
#include "unit_cell.h"

//void readLammps(const char *dataName, int skip);
void writeDump(const char *dataName, int step, char flag, double box[]);
void writeK_global(const char *dataName, int l);
void writeDisp(const char *dataName, char c, int t1, int tstep);
void writeForce(const char *dataname, char c, double p, int tstep, double box[]);
void writeReaction(const char *dataName, char c, int t1, int tstep);
void writeStress(const char *dataName, int t1, int tstep);
void writeStrain(const char *dataName, int t1, int tstep);

void writeInternalForce(const char *dataName, int step);
void writeDamage(const char *dataName, int step);
void writeNeighbor(const char *dataName);
void writeConnection(const char *dataName);
void writeKnTv(const char *dataName);

#endif

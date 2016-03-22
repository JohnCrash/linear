#ifndef __MISC_H__
#define __MISC_H__
#include "linear.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#define NN 3

DYNFUNC void disablePrint(int b);
DYNFUNC int printDiffent(const char * s, real * A, real *B);
DYNFUNC int printDiffent1(const char * s, real * b, real *x);
DYNFUNC void printMat(const char * s, real * A);
DYNFUNC void printMat3(real * P, real * A, real *B, int n);
DYNFUNC void printVec(const char * s, real *v);
DYNFUNC void copyMatrix(real * des, real * src);
DYNFUNC real * makeMatrix();
DYNFUNC real * makeRandMatrix();
DYNFUNC real * makeRandVec();

DYNFUNC real randNegative();
/*
 * 允许有正负数
 */
DYNFUNC real * makeRandMatrix2();
DYNFUNC real * makeRandVec2();
DYNFUNC real randomReal();

DYNFUNC real * makeRandSPDMatrix();
DYNFUNC real * makeRandSPDMatrixNUB(int nub);

DYNFUNC void freeMatrix(real * A);

DYNFUNC double getClock();
#endif
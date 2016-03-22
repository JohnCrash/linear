#ifndef __MISC_H__
#define __MISC_H__
#include "linear.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#define NN 3

void disablePrint(int b);
int printDiffent(const char * s,real * A, real *B);
int printDiffent1(const char * s,real * b,real *x);
void printMat(const char * s,real * A);
void printMat3(real * P,real * A,real *B,int n);
void printVec(const char * s,real *v);
void copyMatrix(real * des,real * src);
real * makeMatrix();
real * makeRandMatrix();
real * makeRandVec();

real randNegative();
/*
 * 允许有正负数
 */
real * makeRandMatrix2();
real * makeRandVec2();
real randomReal();

real * makeRandSPDMatrix();
real * makeRandSPDMatrixNUB(int nub);

void freeMatrix(real * A);

double getClock();
#endif
#ifndef __MISC_H__
#define __MISC_H__
#include "linear.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#define N 5

void disablePrint(int b);
int printDiffent(const char * s,real * A, real *B);
void printMat(const char * s,real * A);
void copyMatrix(real * des,real * src);
real * makeMatrix();
real * makeRandMatrix();
void freeMatrix(real * A);

double getClock();
#endif
#ifndef _LCP_H_
#define _LCP_H_
#include "linear.h"

/*
 * 使用Gauss-Seidel迭代法解LCP
 * x'(Ax+b) = 0互补是互补条件，x'是x的转置
 */
int lcp_pgs(real * A,real *b,real *x,int n);
int Solve_GaussSeidel(real * A, real * b, real *x,int n,int kMax);

#endif

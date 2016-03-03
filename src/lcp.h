#ifndef _LCP_H_
#define _LCP_H_
#include "linear.h"
#include <vector>

/*
 * 将线性互补问题全部解穷举出来（因此算法复杂度为O(pow(2,n)*pow(n,3))）
 * y = Ax+b;x,y>=0;x'y=0
 * 成功返回1并将全部的解都放入到xs中,失败返回0
 * 注意返回的解x,y是混合放置的，后面接互补集合，互补集合为1对于x,0对应y
 * 如a b c d 1 0 1 1,那么x=a 0 c d,y=0 c 0 0
 */
int lcp(real *A,real *b,std::vector<real *>& xs,int n);
/*
 * 混合线性互补问题
 * 标准互补问题的基础上前nub个x变量是自由变量(可以为任何数)，而对于的y=0
 * 如果说标准互补需要x=0,y>=0或者x>=0,y=0，而混合互补中有nub个变量是x是自由数,y=0
 */
int mlcp(real * A,real *b,std::vector<real *>& xs,int nub,int n);

int freeLcpSolve(std::vector<real *>& xs);
/*
 * 使用Gauss-Seidel迭代法解LCP
 * x'(Ax+b) = 0互补是互补条件，x'是x的转置
 * nMax最大迭代次数,acc收敛精度
 */
int lcp_pgs(real * A,real *b,real *x,int n,int nMax,real acc);
int mlcp_pgs(real * A,real *b,real *x,int nub,int n,int nMax,real acc);

int Solve_GaussSeidel(real * A, real * b, real *x,int n,int kMax);
 
void sub_line(real * M,int i,int j,int n,int skip);
void add_line(real * M,int i,int j,int n,int skip);
void multiply_line(real * M,real d,int i,int n,int skip);
void negative_line(real * M,int j,int n,int skip);
void elimination(real *M,int eli,int row,int col,int n,int skip);

int lcp_pivot( real *A,real *b,real *x,int n);

int lcp_lemke(real * A,real *b,real *x,int n);

typedef int ( *lcpSolver)(real *A,real *b,real *x,int n);

int mlcpSolver(real * A,real *b,real *x,int nub,int n,lcpSolver lcpfunc);

void printM(real * M,int * N,int n,int skip);
#endif

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
 * 成功返回1,失败返回0
 */
int lcp(real *A,real *b,std::vector<real *>& xs,int n);
/*
 * 混合线性互补问题
 * 标准互补问题的基础上前nub个x变量是自由变量(可以为任何数)，而对应的y=0
 * 如果说标准互补需要x=0,y>=0或者x>=0,y=0，
 * 而混合互补中有nub个变量是x是自由数(可以为任何数)而必须y=0
 * 下面函数使用穷举法
 * 成功返回1,失败返回0
 */
int mlcp(real * A,real *b,std::vector<real *>& xs,int nub,int n);

/* 释放数组xs */
int freeLcpSolve(std::vector<real *>& xs);

/*
 * 使用Gauss-Seidel迭代法解LCP
 * y=Ax+b,y'x=0,x>=0,y>=0
 * A是一个nxn矩阵，b是一个nx1向量
 * nMax最大迭代次数,acc收敛精度
 * x是一个nx1向量，成功返回x
 * 成功返回1,失败返回0
 */
int lcp_pgs(real * A,real *b,real *x,int n,int nMax,real acc);
int mlcp_pgs(real * A,real *b,real *x,int nub,int n,int nMax,real acc);

int Solve_GaussSeidel(real * A, real * b, real *x,int n,int kMax);
 
 /*
  * sub_line,add_line 表示M是一个nxn矩阵,行距是skip(每行的元素数)
  * 将第i行减去或者加上第j行
  */
void sub_line(real * M,int i,int j,int n,int skip);
void add_line(real * M,int i,int j,int n,int skip);
/*
 * M是一个nxn矩阵，skip是行距(每行的元素数)
 * 将第i行的全部元素都乘以d (e*=d)
 */
void multiply_line(real * M,real d,int i,int n,int skip);
/*
 * M是一个nxn矩阵，skip是行距(每行的元素数)
 * 将第j行的元素都取反(e=-e)
 */
void negative_line(real * M,int j,int n,int skip);
/*
 * M是一个nxn矩阵，skip是行距(每行的元素数)
 * 
 */
void elimination(real *M,int eli,int row,int col,int n,int skip);

/*
 * 使用Principad pivot method求解LCP
 * y=Ax+b,y'x=0,x>=0,y>=0
 * A是一个nxn矩阵，b是一个nx1向量
 * x是一个2*n个元素的real数组，前n个返回x,后n个返回y
 * 函数成功返回1，失败返回0
 */
int lcp_pivot( real *A,real *b,real *x,int n);
real *mallocPivotMatrix(real * A,real *b,int n,int m,int *pskip);
/*
 * 使用lemke method求解LCP
 * y=Ax+b,y'x=0,x>=0,y>=0
 * A是一个nxn矩阵，b是一个nx1向量
 * x是一个2*n个元素的real数组，前n个返回x,后n个返回y
 * 函数成功返回1，失败返回0 
 */
int lcp_lemke(real * A,real *b,real *x,int n);
real *mallocLemkeMatrix(real * A,real *b,int n,int m,int *pskip);

/*
 * 将问题放入到一个统一的矩阵M中，进行运算。
 * 下面是M矩阵的结构,其中a是辅助变量，lemke算法需要，pivot不需要这一列。
 * |	y	|	x	| a | b |
 * |	I	|	-A	|	| b | I是nxn单位矩阵,A，b就是y=Ax+b
 * |		|		| 0	|	| m行参与pivot的行
 * 在运算停止时将解放入最后一行。
 * x是一个2xn数组，前面的n个元素是解的x部分，剩下的是y部分。
 * 函数成功返回1，失败返回0 
 */
int lcp_pivotBlock(real *M,real *x,int n,int m,int skip);
int lcp_lemkeBlock(real *M,real *x,int n,int m,int skip);

enum LCPSolver{
	PGS,
	LEMKE,
	PIVOT,
};

/*
 * 求解混合互补问题mlcp
 * A是一个nxn矩阵，b是一个nx1向量
 * x是一个2*n个元素的real数组，前n个返回x,后n个返回y
 * nub表示A的前nub个列，x是任意数，y=0
 */
int mlcpSolver(real * A,real *b,real *x,int nub,int n,LCPSolver solver);

void printM(real * M,int * N,int n,int skip);
#endif

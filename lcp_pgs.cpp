#include "lcp.h"

/*
 *使用Gauss-Seidel迭代法解LCP
 * x'(Ax+b)=0 x'表示x转置
 * x>=0 and Ax+b>=0
 * x是迭代初始值
 * 解法不停的迭代下面公式直到x(k+1)收敛
 * max(0,(M'(Nx(k)-b))) = x(k+1) ,M'是M的逆，x(k)第k次迭代
 */
int lcp_pgs(real * A,real *b,real *x,int n)
{
}

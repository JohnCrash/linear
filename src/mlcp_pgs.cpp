#include "lcp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "misc.h"

/*
 *使用Gauss-Seidel迭代法解LCP
 * Gauss-Seidel迭代成立的条件是:
 * A是一个对称正定矩阵或者A是一个严格主元占优的
 * x'(Ax+b)=0 x'表示x转置
 * x>=0 and Ax+b>=0
 * x是迭代初始值
 * 解法不停的迭代下面公式直到x(k+1)收敛
 * max(0,(M'(Nx(k)-b))) = x(k+1) ,M'是M的逆，x(k)第k次迭代
 */
int mlcp_pgs(real * A,real *b,real *x,int nub,int n,int nMax,real acc)
{
	int i,j,k;
	real dot,c,z;
	/*
	 * (1) A = M-N ,M=daig(daig(A)),N=-(tril(A,-1)+triu(A,1))
	 * (2) A = M-N ,M=daig(daig(A))+trilA,-1),N=-triu(A,1) Kenny Erleben的方法
	 * 直接将A=M'N,b=M'b
	 */
	for(i=0;i<n;i++){
		dot = A[i*n+i];
		if(dot==0){
			return 0;
		}
		b[i] /= dot;
		for(j=0;j<n;j++){
			if(i==j){
				A[i*n+j] = 0;
			}else{
				A[i*n+j]/=-dot;
			}
		}
	}
	k = 0;
	do{
		/*
		 * c是abs(x[k]-x[k+1]),中最大的那个
		 * 如果c越来越小，并且低于一个值就认为收敛。
		 * 因为x>=0是解的条件
		 */
		c = 0;
		for(i=0;i<n;i++){
			dot = 0;
			for(j=0;j<n;j++){
				dot += (A[i*n+j]*x[j]);
			}
			if(i<nub){
				z = dot-b[i];
				c = fmax(c,fabs(z));
			}else{
				z = fmax(0,dot-b[i]);
				c = fmax(c,fabs(x[i]-z));
			}
			x[i] = z;
		}
		/*
		 * 判断收敛
		 */
		if(c<acc)break;
	}while(k++<nMax);

	return k;
}

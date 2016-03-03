 #include "linear.h"
 #include "lcp.h"
 #include <stdlib.h>
 #include <memory.h>
 #include <stdio.h>
 
 /*
  * 对A的左上角(nub)x(nub)的区域进行消元处理，将该区块处理成单元矩阵
  * 
  */
 static int pivotNub(real *A,real *b,int *P,int nub,int n)
 {
	int skip = n;
	int i,j;
	for(i=0;i<nub;i++){
		if( A[i*skip+i] == 0 )return 0; //FIXBUG,使用交换行来解决
		real d = 1.0/A[i*skip+i];
		multiply_line(A,d,i,n,skip);
		A[i*skip+i] = 1;
		for(j=0;j<n;j++){
			if(j!=i&&A[j*skip+i]!=0){
				elimination(A,j,i,i,n,skip);
			}
		}
	}
	return 1;
 }
 
/*
 * M表示如下(M|y,x|'=0,'表示转置)，nub为自由区。
 * 因为在nub区x=R(任意实数)y=0，因此y的nub区被排除
 * x,y的base区中选择n-nub个基列和x的nub区组成方程得到解
 * 也就是y base区和x base区互补选择n-nub个列与x nub区组成最终的求解矩阵
 * 先将M弄成下面的样式
 * |			y				|			x				|d|b|
 * |	nub	|	base	|	nub	|	base	|d|b|
 * |			|0			0	|1			|				|
 * |			|0			0	|		1	|				|
 * |			|1			0	|			|				|
 * |			|	1			|			|				|
 * |			|		1		|			|				|
 * |			|0			1	|			|				|
 * 当把M转化为这样以后，y的base区和x的base区转化为一个标准线性互补问题
 */
int mlcpSolver(real * A,real *b,real *x,int nub,int n,lcpSolver lcpfunc)
{
	int i,j,skip;
	int * P;
	int result = 0;

	if(nub<0||nub>n)return result;
	P = (int *)malloc(nub*sizeof(int));
	/*
	 * nub = n问题完全转化为线性方程问题
	 * nub = 0问题是一个完全的线性互补问题
	 */
	for(i = 0;i<nub;i++)
		P[i] = i;
	if( pivotNub(A,b,P,nub,n) ){
		real *AA = (real *)malloc((n-nub)*(n-nub)*sizeof(real));
		real * bb = (real *)malloc((n-nub)*sizeof(real));
		real * xx = (real *)malloc(2*(n-nub)*sizeof(real));
		result = lcpfunc(AA,bb,xx,n-nub);
		/*
		 * 收集解
		 */
		for(i=0;i<nub;i++){
			x[i] = b[i]; //nub区x
			x[n+i] = 0; //nub区y
		}
		for(i=0;i<n-nub;i++){
			x[nub+i] = xx[i]; //互补区x
			x[n+nub+i] = xx[n-nub+i]; //互补区y
		}
		free(xx);
		free(bb);
		free(AA);
	}
	free(P);
	return result;
}
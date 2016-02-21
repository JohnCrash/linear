#include "lcp.h"
#include <math.h>

#define max(a,b) (a>b?a:b)
#define abs(a) (a>0?a:-a)
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
	int i,j,k;
	real d,z,c,e;
	/*
	 * A = M+N ,M=daig(A),N=tril(A,-1)+triu(A,1)
	 * 直接将A=M'N,b=M'b
	 */
	for(i=0;i<n;i++){
		d = A[i*n+i];
		if(d==0)return 0;
		b[i] /= d;
		for(j=0;j<n;j++){
			if(i==j){
				A[i*n+j] = 0;
			}else{
				A[i*n+j]/=d;
			}
		}
	}
	/*
	 * 开始迭代
	 */
	while(true){
		/*
		 * c是abs(x[k]-x[k+1]),中最大的那个
		 * 如果c越来越小，并且低于一个值就认为收敛。
		 * e是所有x中负数的和，如果e>=0表示在正常范围
		 * 因为x>=0是解的条件
		 */
		c = 0;
		e = 0;
		for(i=0;i<n;i++){
			d = 0;
			for(j=0;j<n;j++){
				d += (A[i*n+j]*x[j]-b[j]);
			}
			z = max(0,d);
			c = max(c,abs(x[i]-z));
			x[i] = z;
			if(x[i]<0)
				e+=x[i];
		}
		/*
		 * 判断收敛
		 */
		if(e>=0 && c<0.01)
			return 1;
		/*
		 * 判断发散,暂时不知道如何实现
		 */
		if(k>100)
			return 0;
		k++;
	}
	return 0;
}

#include "lcp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "misc.h"

#define max(a,b) (a>b?a:b)
#define abs(a) ((a)>0?(a):-(a))

static void printX(real * x,int n)
{
	for(int i=0;i<n;i++){
		printf("%.4f\t",x[i]);
	}
	printf("\n");	
}

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
	real * y = (real *)malloc(n*sizeof(real));
	/*
	 * A = M+N ,M=daig(A),N=tril(A,-1)+triu(A,1)
	 * 直接将A=M'N,b=M'b
	 */
	for(i=0;i<n;i++){
		d = A[i*n+i];
		if(d==0){
			free(y);
			return 0;
		}
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
	 printf("========================\n");
	  printMat("solve lcp A=",A);
	  printVec("lcp b=",b);
	 printVec("lcp x=",x);
	 printf("---------------------------------------------\n");
	k = 0;
	while(true){
		printf("pgs[%d]:",k);
		printX(x,n);
		/*
		 * c是abs(x[k]-x[k+1]),中最大的那个
		 * 如果c越来越小，并且低于一个值就认为收敛。
		 * e是所有x中负数的和，如果e>=0表示在正常范围
		 * 因为x>=0是解的条件
		 */
		c = 0;
		e = 0;
		memcpy(y,x,n*sizeof(real));
		for(i=0;i<n;i++){
			d = 0;
			for(j=0;j<n;j++){
				d += (A[i*n+j]*x[j]);
			}
			x[i] = d-b[i];
			//x[i] = b[i]-d;
		}
		printX(x,n);
		for(i=0;i<n;i++){
			x[i] = max(0,x[i]);
			c = max(c,abs(x[i]-y[i]));
			if(x[i]<0)
				e+=x[i];
		}		
		/*
		 * 判断收敛
		 */
		if(e>=0 && c<0.01){
			printf("pgs[%d]:",k);
			printX(x,n);			
			printf("---------------------------------------------\n");
			printf("pgs k=%d,e=%f,c=%f\n",k,e,c);
			printX(x,n);
			printf("---------------------------------------------\n");
			free(y);
			return k+1;
		}
		/*
		 * 判断发散,暂时不知道如何实现
		 */
		if(k>100)
			break;
		k++;
	}
	free(y);
	return 0;
}

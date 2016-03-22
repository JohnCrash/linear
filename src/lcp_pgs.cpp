#include "lcp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "misc.h"

static void printX(real * x,int n)
{
	for(int i=0;i<n;i++){
		printf("%.4f\t",x[i]);
	}
	printf("\n");	
}

static void printCol(real s,int j,int n)
{
	if(j==0){
		if(s>=0)
			printf("| %.2f\t",s);
		else
			printf("|%.2f\t",s);
	}else if(j==n-1){
		if(s>=0)
			printf(" %.2f\t|",s);
		else
			printf("%.2f\t|",s);
	}else{
		if(s>=0)
			printf(" %.2f\t",s);	
		else
			printf("%.2f\t",s);	
	}
}

static void printGussSolve2(real * A,real *b,int n)
{
	real  x[4][2];
	int c = 0;
	x[c][0] = 0; //[0]
	x[c][1] = 0; //[0]
	c++;
	if(A[0*n+0]!=0){
		x[c][0] = -b[0]/A[0*n+0];//[1]
		x[c][1] = 0;//[0]
		c++;
	}
	if(A[1*n+1]!=0){
		x[c][0] = 0;//[0]
		x[c][1] = -b[1]/A[1*n+1];//[1]	
		c++;
	}
	real d = A[0*n+0]*A[1*n+1]-A[0*n+1]*A[1*n+0];
	if(d!=0){
		x[3][0] = (b[1]*A[0*n+1]-b[0]*A[1*n+1])/d;//[1]
		x[3][1] = (b[0]*A[1*n+0]-b[1]*A[0*n+0])/d;//[1]
		c++;
	}
	for(int j=0;j<n;j++){
		real y[2];
		
		for(int i=0;i<c;i++){
			y[0] = A[0*n+0]*x[i][0]+A[0*n+1]*x[i][1]+b[0];
			y[1] = A[1*n+0]*x[i][0]+A[1*n+1]*x[i][1]+b[1];
			if(x[i][j]>0)
				printf("[ %.2f,\t",x[i][j]);
			else
				printf("[%.2f,\t",x[i][j]);			
			if(y[j]>0)
				printf(" %.2f]",y[j]);
			else
				printf("%.2f]",y[j]);			
		}
		printf("\n");
	}
}

/*
 *使用Gauss-Seidel迭代法解LCP
 * Gauss-Seidel迭代成立的条件是:
 * A是一个对称正定矩阵或者A是一个严格主元占优的
 * x'(Ax+b)=0 x'表示x转置
 * x>=0 and Ax+b>=0
 * x是迭代初始值
 * 解法不停的迭代下面公式直到x(k+1)收敛
 * max(0,(M'(Nx(k)-b))) = x(k+1) ,M'是M的逆，x(k)第k次迭代
 *下面使用的是Kenny Erleben的PGS
 */
int lcp_pgs2(real * A,real *b,real *x,int n)
{
	int i,j,k;
	real d,c;
	real * y = (real *)malloc(n*sizeof(real));
	real * AA = (real *)malloc(n*n*sizeof(real));
	real * bb = (real *)malloc(n*sizeof(real));
	memcpy(AA,A,n*n*sizeof(real));
	memcpy(bb,b,n*sizeof(real));
	
	/*
	 * (1) A = M-N ,M=daig(daig(A)),N=-(tril(A,-1)+triu(A,1))
	 * (2) A = M-N ,M=daig(daig(A))+trilA,-1),N=-triu(A,1) Kenny Erleben的方法
	 * 直接将A=M'N,b=M'b
	 */
	for(i=0;i<n;i++){
		d = A[i*n+i];
		if(d==0){
			free(y);
			free(AA);
			free(bb);
			return 0;
		}
		b[i] /= d;
		for(j=0;j<n;j++){
			if(i==j){
				A[i*n+j] = 0;
			}else{
				A[i*n+j]/=-d;
			}
		}
	}
	k = 0;
	while(true){
	//	printf("pgs[%d]:",k);
	//	printX(x,n);
		/*
		 * c是abs(x[k]-x[k+1]),中最大的那个
		 * 如果c越来越小，并且低于一个值就认为收敛。
		 * 因为x>=0是解的条件
		 */
		c = 0;
		//memcpy(y,x,n*sizeof(real));
		for(i=0;i<n;i++){
			d = 0;
			for(j=0;j<n;j++){
				d += (A[i*n+j]*x[j]);
			}
			y[i] = d-b[i];
			//y[i] = b[i]-d;
		}
		//printX(y,n);
		for(i=0;i<n;i++){
			y[i] = (real)fmax(0,y[i]);
			c = fmax(c,fabs(x[i]-y[i]));
			x[i] = y[i];
		}		
		/*
		 * 判断收敛
		 */
		if(c<0.01){
			free(y);
			free(AA);
			free(bb);				
			return k+1;
		}
		/*
		 * 判断发散,暂时不知道如何实现
		 */
		if(k>15)
			break;
		k++;
	}
	free(y);
	free(AA);
	free(bb);
	return 0;
}

/*
 *使用Gauss-Seidel迭代法解LCP
 * Gauss-Seidel迭代成立的条件是:
 * A是一个对称正定矩阵或者A是一个严格主元占优的
 * x'(Ax+b)=0 x'表示x转置
 * x>=0 and Ax+b>=0
 * x是迭代初始值
 * 解法不停的迭代下面公式直到x(k+1)收敛
 * max(0,(M'(Nx(k)-b))) = x(k+1) ,M'是M的逆，x(k)第k次迭代
 *下面使用的是Kenny Erleben的PGS
 */
int lcp_pgs(real * A,real *b,real *x,int n,int nMax,real acc)
{
	int i,j,k;
	real d,c,z;
	/*
	 * (1) A = M-N ,M=daig(daig(A)),N=-(tril(A,-1)+triu(A,1))
	 * (2) A = M-N ,M=daig(daig(A))+trilA,-1),N=-triu(A,1) Kenny Erleben的方法
	 * 直接将A=M'N,b=M'b
	 */
	for(i=0;i<n;i++){
		d = A[i*n+i];
		if(d==0){
			return 0;
		}
		b[i] /= d;
		for(j=0;j<n;j++){
			if(i==j){
				A[i*n+j] = 0;
			}else{
				A[i*n+j]/=-d;
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
			d = 0;
			for(j=0;j<n;j++){
				d += (A[i*n+j]*x[j]);
			}
			z = (real)fmax(0,d-b[i]);
			c = fmax(c,fabs(x[i]-z));
			x[i] = z;
		}
		/*
		 * 判断收敛
		 */
		if(c<acc)break;
	}while(k++<nMax);
	return k;
}

int Solve_GaussSeidel(real * A, real * b, real *x,int n,int kMax)
{
  int i, j;
  real delta;
	if(n==2)
		printGussSolve2(A,b,n);
  // Gauss-Seidel Solver
  for (int k = 0; k < kMax; ++k)
  {
    for (i = 0; i < n; ++i)
    {
      delta = 0.0f;

      for (j = 0;   j < i; ++j) delta += A[i*n+j]*x[j];
      for (j = i+1; j < n; ++j) delta += A[i*n+j]*x[j];

      x[i] = (b[i] - delta) / A[i*n+i];
    }
  }

  return 1;
}

/*
 *用来配接mlcpSolver,让lcp_pgs可以使用
 */
int lcp_pgsAdapter(real * A,real *b,real *x,int n)
{
	return lcp_pgs(A,b,x,n,100,0.01);
}
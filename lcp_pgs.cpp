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
 * x'(Ax+b)=0 x'表示x转置
 * x>=0 and Ax+b>=0
 * x是迭代初始值
 * 解法不停的迭代下面公式直到x(k+1)收敛
 * max(0,(M'(Nx(k)-b))) = x(k+1) ,M'是M的逆，x(k)第k次迭代
 */
int lcp_pgs(real * A,real *b,real *x,int n)
{
	int i,j,k;
	real d,z,c;
	real * y = (real *)malloc(n*sizeof(real));
	real * AA = (real *)malloc(n*n*sizeof(real));
	real * bb = (real *)malloc(n*sizeof(real));
	memcpy(AA,A,n*n*sizeof(real));
	memcpy(bb,b,n*sizeof(real));
	
	/*
	 * (1) A = M-N ,M=daig(daig(A)),N=-(tril(A,-1)+triu(A,1))
	 * 直接将A=M'N,b=M'b
	 */
	 /*
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
	}*/
	/*
	 * (2) A = M-N ,M=daig(daig(A))+trilA,-1),N=-triu(A,1)
	 * 直接将A=M'N,b=M'b	
	 */
	real * M = (real*)malloc(n*n*sizeof(real));
	real * NN = (real*)malloc(n*n*sizeof(real));
	real * invM = (real*)malloc(n*n*sizeof(real));
	real * v = (real*)malloc(n*sizeof(real));
	memset(M,0,n*n*sizeof(real));
	memset(NN,0,n*n*sizeof(real));
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if(i>=j)
				M[i*n+j] = A[i*n+j];
			else
				NN[i*n+j] = -A[i*n+j];
		}
	}
	
	inverse(M,invM,n);
	
	multiply0(A,invM,NN,n,n,n);
	multiply0(v,invM,b,n,n,1);
	memcpy(b,v,n*sizeof(real));
	printf("1\n");
	free(v);
	printf("2\n");
	free(invM);
	printf("3\n");
	free(M);
	printf("4\n");
	free(NN);
	printf("5\n");
	/*
	 * 开始迭代
	 */
	printf("------------------------------------------------\n");
	printf("|\tA\t|b\t|\tM'N\t|M'b\t|\n");
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			printCol(AA[i*n+j],j,n);
		}
		printCol(bb[i],n-1,n);
		for(j=0;j<n;j++){
			printCol(A[i*n+j],j,n);
		}		
		printCol(b[i],n-1,n);
		printf("\n");
	}
	printf("------------------------------------------------\n");
	if(n==2)
		printGussSolve2(AA,bb,n);
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
		//memcpy(y,x,n*sizeof(real));
		for(i=0;i<n;i++){
			d = 0;
			for(j=0;j<n;j++){
				d += (A[i*n+j]*x[j]);
			}
			y[i] = d-b[i];
			//y[i] = b[i]-d;
		}
		printX(y,n);
		for(i=0;i<n;i++){
			y[i] = max(0,y[i]);
			c = max(c,abs(x[i]-y[i]));
			x[i] = y[i];
		}		
		/*
		 * 判断收敛
		 */
		if(k>10&&c<0.01){
			printf("pgs[%d]:",k);
			printX(x,n);			
			printf("---------------------------------------------\n");
			printf("pgs k=%d,c=%f\n",k,c);
			printX(x,n);
			printf("---------------------------------------------\n");
			free(y);
			free(AA);
			free(bb);				
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
	free(AA);
	free(bb);
	return 0;
}

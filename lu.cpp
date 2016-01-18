#include "linear.h"
#include <stdio.h>

void printMat(real * A, int n, int m)
{
	int i, j;
	for (i = 0; i<n; i++){
		for (j = 0; j<m; j++){
			printf("%f ", A[i*n + j]);
		}
		printf("\n");
	}
}

void zero(real * A,int n,int m)
{
	real *a,*aa;
	aa = A+n*m;
	for(a=A;a<aa;a++)
		*a=0;
}

void identity(real * A,int n)
{
	int i;
	zero(A,n,n);
	for(i=0;i<n;i++){
		A[i*(n+1)] = 1;
	}
}

/* 交互i和j行上的全部元素 */
static void xchangeRaw(real * A,int n,int i,int j)
{
	int x;
	real temp;
	real *a,*b;
	a = A+i*n;
	b = A+j*n;
	for(x=0;x<n;x++){
		temp = a[x];
		a[x]=b[x];
		b[x] = temp;
	}
}

/* 查找raw,col下面(同一列)全部元素中绝对值最大的行 */
static int absMaxLeading(real * A,int n,int raw,int col)
{
	int i;
	real v;
	int maxcol=-1;
	real mv=0;
	for(i=raw;i<n;i++){
		v = fabs(A[i*n+col]);
		if(v>mv){
			maxcol=i;
			mv=v;
		}
	}
	return maxcol;
}

/* 消元，将低i行对应元素乘以d加到j上。 */
static void decol(real * A,int n,int i,int j,real d)
{
	int k;
	for(k=0;k<n;k++){
		A[j*n+k] += A[i*n+k]*d;
	}
}

/* 一个n*n低三角矩阵乘一个单位矩阵，该矩阵的i,j(i>j)位置为d */
static void multiplyL(real * L, int n, int i, int j, real d)
{
	for (int k = i; k < n; k++){
		L[k*n + j] += d*L[k*n+i];
	}
}

/*
 * 进行lu分解，将U存入A中，将L存入到L中
 *
 */
int lu(real * A,real * L,int n)
{
	int m,i,j;
	real mr,v,d;
	identity(L,n);
	for(i=0;i<n-1;i++){
		printf("A:\n");
		printMat(A, n, n);
		m=absMaxLeading(A,n,i,i);
		printf("absMaxLeading(A,%d,%d,%d)\n", n, i, i);
		printMat(A, n, n);
		if(m!=i)
			xchangeRaw(A,n,i,m);
		printf("xchangeRaw(A,%d,%d,%d)\n", n, i, m);
		printMat(A, n, n);
		mr = A[i*n+i];
		if(mr==0)
			return 0;
		for(j=i+1;j<n;j++){
			v = A[j*n+i];
			if(v!=0){
				d = -v/mr;
				printf("A: d = %f ,v = %f,mr = %f \n",d,v,mr);
				printMat(A, n, n);
				decol(A,n,i,j,d);
				printf("decol(A,%d,%d,%d,%f)\n", n, j, i, d);
				printMat(A, n, n);
				A[j*n+i]=0;
				printf("L=\n");
				printMat(L, 3, 3);
				multiplyL(L, n, j, i, -d);
				printf("decol(L,%d,%d,%d,%f)\n", n, j, i, -d);
				printMat(L, 3, 3);
			}
		}
	}
	return 1;
}

/* */
int pldu(real * A, real * P,real * D,real * L, int n)
{
	int m, i, j;
	real mr, v, d;
	identity(L, n);
	for (i = 0; i<n - 1; i++){
		printf("A:\n");
		printMat(A, n, n);
		m = absMaxLeading(A, n, i, i);
		printf("absMaxLeading(A,%d,%d,%d)\n", n, i, i);
		printMat(A, n, n);
		if (m != i)
			xchangeRaw(A, n, i, m);
		printf("xchangeRaw(A,%d,%d,%d)\n", n, i, m);
		printMat(A, n, n);
		mr = A[i*n + i];
		if (mr == 0)
			return 0;
		for (j = i + 1; j<n; j++){
			v = A[j*n + i];
			if (v != 0){
				d = -v / mr;
				printf("A: d = %f ,v = %f,mr = %f \n", d, v, mr);
				printMat(A, n, n);
				decol(A, n, i, j, d);
				printf("decol(A,%d,%d,%d,%f)\n", n, j, i, d);
				printMat(A, n, n);
				A[j*n + i] = 0;
				printf("L=\n");
				printMat(L, 3, 3);
				multiplyL(L, n, j, i, -d);
				printf("decol(L,%d,%d,%d,%f)\n", n, j, i, -d);
				printMat(L, 3, 3);
			}
		}
	}
	return 1;
}

int main(int argn,const char *argv[])
{
	real A[]={9,2.5,3,2,5,7,3,5,3};
	real L[9],C[9];
	printf("A=\n");
	printMat(A,3,3);
	lu(A,L,3);
	printf("U=\n");
	printMat(A,3,3);
	printf("L=\n");
	printMat(L,3,3);
	multiply0(C,L, A,3, 3, 3);
	printf("C=\n");
	printMat(C, 3, 3);
	return 0;
}


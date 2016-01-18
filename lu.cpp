#include "linear.h"
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>

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

/* 转置矩阵 */
void transpose(real * A, int n)
{
	real temp;
	for (int i = 1; i < n; i++){
		for (int j = 0; j < i; j++){
			temp = A[i*n + j];
			A[i*n + j] = A[j*n+i];
			A[j*n+i] = temp;
		}
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
		m=absMaxLeading(A,n,i,i);
		if(m!=i)
			xchangeRaw(A,n,i,m);
		mr = A[i*n+i];
		if(mr==0)
			return 0;
		for(j=i+1;j<n;j++){
			v = A[j*n+i];
			if(v!=0){
				d = -v/mr;
				decol(A,n,i,j,d);
				A[j*n+i]=0;
				multiplyL(L, n, j, i, -d);
			}
		}
	}
	return 1;
}

/* 
 * 进行lu分解，将U存入A中，将L存入到L中。P是交换矩阵。
 * D是一个对角线矩阵，D就是为了使U的对角线都是1。
 * A = P*L*D*U
 */
int pldu(real * A, real * P,real * D,real * L, int n)
{
	int m, i, j;
	real mr, v, d;
	identity(L, n);
	identity(P, n);
	identity(D, n);
	for (i = 0; i<n - 1; i++){
		m = absMaxLeading(A, n, i, i);
		if (m != i){
			xchangeRaw(A, n, i, m);
			xchangeRaw(P, n, i, m);
		}
		mr = A[i*n + i];
		if (mr == 0)
			return 0;
		for (j = i + 1; j<n; j++){
			v = A[j*n + i];
			if (v != 0){
				d = -v/mr;
				decol(A, n, i, j, d);
				A[j*n + i] = 0;
				multiplyL(L, n, j, i, -d);
			}
		}
	}
	/* 处理U的对角线 */
	for (i = 0; i < n; i++){
		mr = A[i*n + i];
		A[i*n + i] = 1;
		for (j = i + 1; j < n; j++){
			A[i*n + j] /= mr;
		}
		D[i*n + i] = mr;
	}
	return 1;
}

/* 求下三角矩阵的逆矩阵 */
static void inverse_low_triangle(real * L, int n)
{
	real * A, *B,*C;
	A = (real *)malloc(n*n*sizeof(real));
	B = (real *)malloc(n*n*sizeof(real));
	C = (real *)malloc(n*n*sizeof(real));
	identity(A,n);
	identity(B,n);
	for (int i = 0; i < n; i++){
		if (i == 0){
			for (int j = i + 1; j < n; j++){
				B[j*n + i] = -L[j*n + i];
			}
			multiply0(C, A, B, n, n, n);
			for (int j = i + 1; j < n; j++){
				B[j*n + i] = 0;
			}
			memcpy(A, C, n*n*sizeof(real));
		}
		else{
			for (int j = i + 1; j < n; j++){
				A[j*n + i] = -L[j*n + i];
			}
		}
	}
	memcpy(L, A, n*n*sizeof(real));
	free(A);
	free(B);
	free(C);
}

/* 求上三角逆矩阵 */
static void inverse_upper_triangle(real * L, int n)
{

}

/* 求交换矩阵的逆矩阵 */
static void inverse_pivoting(real * L, int n)
{
	transpose(L, n);
}

static void inverse_diagonal(real * D, int n)
{
	for (int i = 0; i < n; i++){
		D[i*n + i] = 1 / D[i*n + i];
	}
}
/*
 * A的逆矩阵，A=P*L*D*U就是上面pldu分解的结果。
 * 最后将结果放入到P矩阵中。
 */
void inverse(real * P, real * L, real * D, real * U, int n)
{
	real * T;

	/*
	 * A' = (P*L*D*U)' = P'*L'*D'*U'  '代表逆
	 */
	T = (real *)malloc(n*n*sizeof(real));
	/* 对角矩阵的逆矩阵 */
	inverse_diagonal(D, n);
	/* 三角的逆矩阵 */
	inverse_low_triangle(L, n);
	inverse_upper_triangle(L, n);
	/* 交换矩阵的逆矩阵 */
	inverse_pivoting(P, n);

	multiply0(T, P, L, n, n, n);
	multiply0(P, T, D, n, n, n);
	multiply0(T, P, U, n, n, n);
	memcpy(P, T, n*n*sizeof(real));
	free(T);
}

static void test_lu_1()
{
	real A[] = { 1, 2, 3, 2, 5, 7, 3, 5, 3 };
	real L[9], C[9];
	printf("A=\n");
	printMat(A, 3, 3);
	lu(A, L, 3);
	printf("U=\n");
	printMat(A, 3, 3);
	printf("L=\n");
	printMat(L, 3, 3);
	multiply0(C, L, A, 3, 3, 3);
	printf("C=\n");
	printMat(C, 3, 3);
}

static void test_pldu_1()
{
	real A[] = { 1, 2, 3, 2, 5, 7, 3, 5, 3 };
	real L[9], C[9];
	real P[9], D[9], T[9];
	printf("A=\n");
	printMat(A, 3, 3);
	pldu(A, P,D,L, 3);
	printf("P=\n");
	printMat(P, 3, 3);
	printf("D=\n");
	printMat(D, 3, 3);
	printf("U=\n");
	printMat(A, 3, 3);
	printf("L=\n");
	printMat(L, 3, 3);
	multiply0(C, P, L, 3, 3, 3);
	multiply0(T, C, D, 3, 3, 3);
	multiply0(C, T, A, 3, 3, 3);
	printf("C=\n");
	printMat(C, 3, 3);
}

static void test_inverse_1()
{
	real A[] = { 1, 2, 3, 2, 5, 7, 3, 5, 3 };
	real L[9], C[9];
	real P[9], D[9], T[9];
	printf("A=\n");
	memcpy(T, A, 9 * sizeof(real));
	printMat(A, 3, 3);
	pldu(A, P, D, L, 3);
	printf("P=\n");
	printMat(P, 3, 3);
	printf("D=\n");
	printMat(D, 3, 3);
	printf("U=\n");
	printMat(A, 3, 3);
	printf("L=\n");
	printMat(L, 3, 3);
	inverse(P, L, D, A, 3);
	printf("inverse=\n");
	printMat(P, 3, 3);
	multiply0(C, T, P, 3, 3, 3);
	printf("A*A'=\n");
	printMat(C, 3, 3);
}

int main(int argn,const char *argv[])
{
	test_inverse_1();
	return 0;
}


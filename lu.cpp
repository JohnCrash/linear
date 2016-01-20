#include "linear.h"
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>

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
void xchangeRaw(real * A,int n,int i,int j)
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
int absMaxLeading(real * A,int n,int raw,int col)
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

/* 消元，将第i行对应元素乘以d加到j上。仅仅做i右边的 */
static void decol(real * A,int n,int i,int j,real d)
{
	int k;
	for(k=i+1;k<n;k++){
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

static void printMat(const char * s, real * A,int N)
{
	int i, j;
	if (s)
		printf("%s\n", s);
	for (i = 0; i<N; i++){
		for (j = 0; j<N; j++){
			printf("%f ", A[i*N + j]);
		}
		printf("\n");
	}
}

/* 将一列的绝对值最大元素作为主元 */
static void sortPovit(real * A,real * P,int n)
{
	int m;
	for (int i = 0; i < n; i++){
		m = absMaxLeading(A, n, i, i);
		if (m != i&&m!=-1){
			xchangeRaw(A, n, i, m);
			if (P)
				xchangeRaw(P, n, i, m);
		}
	}
}

/*
 * 进行lu分解，将U存入A中，将L存入到L中
 */
int lu(real * A,real * P,real * L,int n)
{
	int i,j;
	real mr,v,d;
	identity(L,n);
	identity(P, n);
	sortPovit(A, P, n);
	for(i=0;i<n-1;i++){//column
		mr = A[i*n+i];
		if(mr==0)
			return 0;
		for(j=i+1;j<n;j++){
			v = A[j*n+i];
			if(v!=0){
				d = v/mr;
				decol(A,n,i,j,-d);
				A[j*n+i]=0;
				multiplyL(L, n, j, i, d);
			}
		}
	}
	transpose(P, n);
	return 1;
}

/* 
 * 进行lu分解，将U存入A中，将L存入到L中。P是交换矩阵。
 * D是一个对角线矩阵，D就是为了使U的对角线都是1。
 * A = P*L*D*U
 */
int pldu(real * A, real * P,real * D,real * L, int n)
{
	int i, j;
	real mr, v, d;
	identity(L, n);
	identity(P, n);
	identity(D, n);
	sortPovit(A, P, n);
	for (i = 0; i<n - 1; i++){
		mr = A[i*n + i];
		if (mr == 0)
			return 0;
		for (j = i + 1; j<n; j++){
			v = A[j*n + i];
			if (v != 0){
				d = v/mr;
				decol(A, n, i, j, -d);
				A[j*n + i] = 0;
				multiplyL(L, n, j, i, d);
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
	transpose(P, n);
	return 1;
}

/* crout algothim */
int crout_lu(real * A,real * L,int n)
{
	int i,j;
	real mr,v;
	identity(L,n);
//	sortPovit(A, NULL, n);
	for(i=0;i<n;i++){
		mr = A[i*n+i];
		if(mr==0)
			return 0;
		L[i*n + i] = mr;
		A[i*n+i] = 1;
		for(j=i+1;j<n;j++){
			A[i*n+j]/=mr;
		}
		for(j=i+1;j<n;j++){
			v = A[j*n + i];
			if(v!=0){
				decol(A,n,i,j,-v);
				A[j*n+i]=0;
				multiplyL(L, n, j, i, v);
			}
		}
	}
	return 1;
}

int crout_plu(real * A,real * P,real * L,int n)
{
	int i,j;
	real mr,v;
	identity(L,n);
	identity(P, n);
	sortPovit(A, P, n);
	for(i=0;i<n;i++){
		mr = A[i*n+i];
		if(mr==0)
			return 0;
		L[i*n + i] = mr;
		A[i*n+i] = 1;
		for(j=i+1;j<n;j++){
			A[i*n+j]/=mr;
		}
		for(j=i+1;j<n;j++){
			v = A[j*n + i];
			if(v!=0){
				decol(A,n,i,j,-v);
				A[j*n+i]=0;
				multiplyL(L, n, j, i, v);
			}
		}
	}
	transpose(P, n);
	return 1;	
}

/* 求下三角矩阵的逆矩阵，算法要求对角线为1 */
void inverse_low_triangle(real * L, int n)
{
	real * A, *B,*C;
	A = (real *)malloc(n*n*sizeof(real));
	B = (real *)malloc(n*n*sizeof(real));
	C = (real *)malloc(n*n*sizeof(real));
	identity(A,n);
	identity(B,n);
	for (int i = n-2; i >=0; i--){
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

/* 求上三角逆矩阵，算法要求对角线为1 */
void inverse_upper_triangle(real * L, int n)
{
	transpose(L, n);
	inverse_low_triangle(L, n);
	transpose(L, n);
}

/* 求交换矩阵的逆矩阵 */
void inverse_pivoting(real * L, int n)
{
	transpose(L, n);
}

void inverse_diagonal(real * D, int n)
{
	for (int i = 0; i < n; i++){
		D[i*n + i] = 1 / D[i*n + i];
	}
}

/*
 * A的逆矩阵，A=P*L*D*U就是上面pldu分解的结果。
 * 最后将结果放入到P矩阵中。
 */
int inverse(real * P, real * L, real * D, real * U, int n)
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
	inverse_upper_triangle(U, n);
	/* 交换矩阵的逆矩阵 */
	inverse_pivoting(P, n);

	multiply0(T, U, D, n, n, n);
	multiply0(U, T, L, n, n, n);
	multiply0(T, U, P, n, n, n);
	memcpy(P, T, n*n*sizeof(real));
	free(T);
	return 0;
}

/*
 *计算A的逆矩阵，将结果放入到B中 
 *将AB组成增广矩阵，将A部分通过初等变换转化为单位矩阵I。
 *这时候B即为A的逆矩阵 inverse0该方法比inverse速度快
*/
int inverse0(real * A,real * B,int n)
{
	int m,i,j,k;
	real mr,v;	
	real * P = (real *)malloc(n*n*sizeof(real));
	identity(B, n);

	/* 对行进行交换，确保主元位置是绝对值最大的。*/
	for (int i = 0; i < n; i++){
		m = absMaxLeading(A, n, i, i);
		if (m != i&&m!=-1){
			xchangeRaw(A, n, i, m);
			xchangeRaw(P, n, i, m);
			xchangeRaw(B, n, i, m);
		}
	}
	/* 使用crout法消下三角 */
	for(i=0;i<n;i++){
		mr = A[i*n+i];
		if(mr==0)
			return 0;
		A[i*n+i] = 1;
		for(j=i+1;j<n;j++){
			A[i*n+j]/=mr;
		}
		for(j=0;j<n;j++){
			B[i*n+j]/=mr;
		}
		for(j=i+1;j<n;j++){
			v = A[j*n + i];
			if(v!=0){
				decol(A,n,i,j,-v);
				A[j*n+i]=0;
				/* 将第i行乘v减第j行 */
				for(k=0;k<n;k++){
					B[j*n+k] -= B[i*n+k]*v;
				}
			}
		}
	}
	/* 消上三角，经过上面crout方法处理后A对角线都为1 */
	for(i=n-1;i>0;i--){
		for(j=i;j>=0;j--){
			/* 
			 *消除A[j][i]，主元A[i][i]=1，A[j][i]-=row(i)*A[j][i] 
			 *考虑到A的主元为1，且j行除了主元其他全是0
			 *因此可以不对A进行任何计算了
			 */
			v = A[j*n+i];
			for(k=0;k<n;k++){
				B[j*n+k] -= B[i*n+k]*v;
			}			
		}
	}
	/* 乘交换矩阵恢复位置 */
	memcpy(A,B,n*n*sizeof(real));
	multiply0(B,P,A,n,n,n);
	free(P);
}
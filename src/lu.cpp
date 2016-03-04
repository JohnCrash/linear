#include "linear.h"
#include "misc.h"
#include <math.h>

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
	a = &A[i*n];
	b = &A[j*n];
	for(x=0;x<n;x++){
		temp = a[x];
		a[x]=b[x];
		b[x] = temp;
	}
}

/* 
 * A是一个下三角矩阵，交互i和j行上的在i列左侧的元素 注意：j>i
 * [a 0 0]					[a 0 0]
 * [b 1 0] 比如交互1,2行	[c 1 0]
 * [c 0 1]					[b 0 1]
 */
void xchangeRawLowerTriangle(real * A,int n,int i,int j)
{
	int x;
	real temp;
	real *a,*b;
	a = &A[i*n];
	b = &A[j*n];
	for(x=0;x<i;x++){
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
	int i,j,m;
	real mr,v,d;
	identity(L,n);
	identity(P, n);
	
	for(i=0;i<n-1;i++){
		m = absMaxLeading(A, n, i, i);
		if (m != i&&m!=-1){
			xchangeRaw(A, n, i, m);
			xchangeRawLowerTriangle(L,n,i,m);
			xchangeRaw(P, n, i, m);
		}	
		mr = A[i*n+i];
		if (FTEQ(mr,0)){
			transpose(P, n);
			return 0;
		}
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
 * crout algothim 
 * http://www.physics.utah.edu/~detar/phys6720/handouts/crout.txt
 * 将n*n矩阵A分解为LU,返回U被写入到A中
  * 成功返回1，失败返回0
  * 下面的算法并不交互行，如果一旦主元为0将停止分解
 */
int crout_lu(real * A,real * L,int n)
{
	int i,j;
	real mr,v;
	identity(L,n);
	for(i=0;i<n;i++){
		mr = A[i*n+i];
		if (FTEQ(mr,0))
			return 0;
		/* 使第i行的主元为1 */
		L[i*n + i] = mr;
		A[i*n+i] = 1;
		for(j=i+1;j<n;j++){
			A[i*n+j]/=mr;
		}
		for(j=i+1;j<n;j++){
			v = A[j*n + i];
			if(v!=0){
				/* 消元操作，同时施加到L与U上 */
				decol(A,n,i,j,-v);
				A[j*n+i]=0;
				multiplyL(L, n, j, i, v);
			}
		}
	}
	return 1;
}

/*
 * A=PLU,其中P为一个交换矩阵,返回U被写入到A中
 * 成功返回1，失败返回0
 */
int crout_plu(real * A,real * P,real * L,int n)
{
	int i,j,m;
	real mr,v;
	identity(L,n);
	identity(P, n);

	for(i=0;i<n;i++){
		m = absMaxLeading(A, n, i, i);
		if (m != i&&m!=-1){
			xchangeRaw(A, n, i, m);
			/* 
			 * 为什么不能使用 xchangeRaw(L, n, i, m); ?
			 * 参见xchangeRawLowerTriangle的说明
			 */
			xchangeRawLowerTriangle(L,n,i,m);
			xchangeRaw(P, n, i, m);
		}
		mr = A[i*n+i];
		if (FTEQ(mr,0)){
			/* 确保失败的时候也保证等式成立 */
			transpose(P, n);
			return 0;
		}
		/* 使第i行的主元为1 */
		L[i*n + i] = mr;
		A[i*n+i] = 1;
		for(j=i+1;j<n;j++){
			A[i*n+j]/=mr;
		}
		for(j=i+1;j<n;j++){
			v = A[j*n + i];
			if(v!=0){
				/* 消元操作，同时施加到L与U上 */
				decol(A,n,i,j,-v);
				A[j*n+i]=0;
				multiplyL(L, n, j, i, v);
			}
		}
	}
	/*
	 * 为什么要加转置操作？
	 * 首先矩阵A1 = P1A ,P1 是第一次交换操作，这时候P=P1
	 * 而A1是第一个用来进行消元的矩阵,依次进行。
	 * A2 = P2A1 = P2P1A ，同时P=P2P1 (在前一个P的基础上做P2)操作
	 * 最后有An = LU = Pn...P2P1A ,同时P=Pn...P2P1
	 * 因为我们需要一个P(外部的)，似的A = PLU 这样P=P‘（逆矩阵）
	 */
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
 * 试验阶段
 * A的逆矩阵，A=P*L*D*U就是上面pldu分解的结果。
 * 最后将结果放入到P矩阵中。
 */
int inverse0(real * P, real * L, real * D, real * U, int n)
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
 *这时候B即为A的逆矩阵 inverse该方法比inverse0速度快
*/
int inverse(real * A,real * B,int n)
{
	int m,i,j,k;
	real mr,v;	
	real * P = (real *)malloc(n*n*sizeof(real));
	identity(B, n);
	identity(P,n);
	
	/* 使用crout法消下三角 */
	for(i=0;i<n;i++){
		m = absMaxLeading(A, n, i, i);
		if (m != i&&m!=-1){
			xchangeRaw(A, n, i, m);
			xchangeRawLowerTriangle(B,n,i,m);
			xchangeRaw(P, n, i, m);
		}		
		mr = A[i*n+i];
		if (FTEQ(mr,0))
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
		for(j=i-1;j>=0;j--){
			/* 
			 *消除A[j][i]，主元A[i][i]=1，A[j][i]-=row(i)*A[j][i] 
			 *考虑到A的主元为1，且j行除了主元其他全是0
			 *因此可以不对A进行任何计算了
			 */
			v = A[j*n+i];
			//A[j*n+i] = 0;
			for(k=0;k<n;k++){
				B[j*n+k] -= B[i*n+k]*v;
			}			
		}
	}
	//transpose(P,n);
	/* 乘交换矩阵恢复位置
	 * 求A的逆矩阵A'，首先将A进行变换B=P'A，通过增广矩阵计算出B'
	 * B'=(P'A)'=PA' ,B'P'=PA'P' , B'P' = A'
	 */
	memcpy(A,B,n*n*sizeof(real));
	multiply0(B,A,P,n,n,n);
	free(P);
	return 1;
}

/*
 * cholesky分解，A=L*L'，L'表示L的转置矩阵，另L是一个n*n下三角矩阵
 * 已知A求L
 */
int cholesky(real * A,real *L,int n)
{
	int i,j,k;
	real sum,s;
	for(k<0;k<n;k++){
		sum=0;
		for(i=0;i<k;i++){
			sum += L[k*n+i]*L[k*n+i];
		}
		if(A[k*n+k]<=sum)
			return 0;
		s = sqrt(A[k*n+k]-sum);
		L[k*n+k] = s;
		for(i=k+1;i<n;i++){
			sum = 0;
			for(j=0;j<k;j++){
				sum += L[i*n+j]*L[k*n+j];
			}
			L[i*n+k] = (A[i*n+k]-sum)/s;
		}
	}
	return 1;
}
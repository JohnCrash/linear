#include "linear.h"
#include "misc.h"

/* L是一个下三角矩阵，Lx=b，求x */
int solve_low_triangle(real * L,real * b,real *x,int n)
{
	real mr;
	for(int i=0;i<n;i++){
		mr = L[i*n+i];
		if(FTEQ(mr,0))
			return 0;
		x[i] = b[i]/mr;
		for(int j=i+1;j<n;j++){
			b[j] -= L[j*n+i]*x[i];
		}
	}
	return 1;
}

/* U是一个上三角矩阵 , Ux=b，求x*/
int solve_upper_triangle(real *U,real * b,real *x,int n)
{
	real mr;
	for(int i=n-1;i>=0;i--){
		mr = U[i*n+i];
		if(FTEQ(mr,0))
			return 0;
		x[i] = b[i]/mr;
		for(int j=0;j<i;j++){
			b[j] -= U[j*n+i]*x[i];
		}
	}
	return 1;
}

/*
 * 解方程Ax=b,A~LU. 其中L是下三角矩阵，U是上三角矩阵
 * y = Ux,Ly=b,先求y然后在求x
 */
 int solve_lu(real * L,real * U,real * b,real *x,int n)
 {
	 if(!solve_low_triangle(L,b,x,n))
		 return 0;
	 if(!solve_upper_triangle(U,x,b,n))
		 return 0;
	 memcpy(x,b,n*sizeof(real));
	 return 1;
 }
 
/*
 * 解方程PLUx=b,其中L是下三角矩阵，U是上三角矩阵
 * LUx=y,Py=b
 */
int solve_plu(real * P,real * L,real * U,real * b,real *x,int n)
{
	transpose(P,n);
	multiply0(x,P,b,n,n,1);
	memcpy(b,x,n*sizeof(real));
	if(!solve_lu(L,U,b,x,n))
		return 0;
	return 1;
}

 /*
  * 求向量a和向量b的点积，a和b都是n维向量
  */
 static real dot(real *a,real *b,int n)
 {
	 real sum = 0;
	 for(int i=0;i<n;i++){
		 sum+=a[i]*b[i];
	 }
	 return sum;
 }
 
 /*
  * 成功返回1,失败返回0
  */
 static int initial(real * A,real *b,int n)
 {
	 int i,j;
	 real d;
	/*
	 * 对A、b的列都除以对角线元素d
	 * 这样迭代式变成x(n+1)=b-A x(n)
	 */
	for(i=0;i<n;i++){
		d = A[i*n+i];
		if(d==0)
			break;
		if(d!=1){
			for(j=0;j<n;j++){
				if(j==i){
					A[i*n+j]=0;
				}else{
					A[i*n+j]/=d;
				}
			}
			b[i]/=d;
		}
	}	 
	return 1;
 }
 
/*
 * Gauss-Seidel迭代法解方程近似解
 * Gauss-Seidel迭代成立的条件是:
 * A是一个对称正定矩阵或者A是一个严格主元占优的
 * https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
 * c 是迭代次数
 * x 迭代的初始变量
 * 失败返回0，成功返回1
 */
int solve_gauss_seidel(real * A,real * b,real * x,int n,int c)
{
	/* 不进行行交换处理，如果有主元为0直接进行下一个步骤 */
	int i;
	
	initial(A,b,n);
	
	for(i=0;i<c;i++){
		/*
		 * 迭代x(n+1)=b-Ax(n)
		 * 将已经解出来的未知数放入继续迭代
		 * 先通过初始的x可以解出x[0]，将x[0]放入到x中继续解x[1]
		 * 一直到解出全部x，这个完成一次迭代
		 */
		for(int j=0;j<n;j++){
			x[j] = b[j]-dot(A+j*n,x,n);
		}
	}
	return 1;
}

/*
 * Jacobi迭代法解方程近似解
 * https://en.wikipedia.org/wiki/Jacobi_method#Description
 * 和gauss_seidel参数相同，但是收敛的慢
 */
int solve_jacobi(real * A,real * b,real * x,int n,int c)
{
	int i,j;
	real *y;
	
	initial(A,b,n);
	
	y = (real*)malloc(n*sizeof(real));
	for(i=0;i<c;i++){
		/*
		 * 迭代x(n+1)=b-Ax(n)
		 */
		multiply0(y,A,x,n,n,1);
		for(j=0;j<n;j++){
			y[j] = b[j]-y[j];
		}
		memcpy(x,y,n*sizeof(real));
	}
	free(y);
	return 1;
}
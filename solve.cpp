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
  * 
  */
 int solve_plu(real * P,real * L,real * U,real * b,real *x,int n)
 {
	 transpose(P,N);
	 multiply0(x,P,b,n,n,1);
	 memcpy(b,x,n*sizeof(real));
	 if(!solve_lu(L,U,b,x,n))
		 return 0;
	 
	// multiply0(b,P,x,n,n,1);
	// memcpy(x,b,n*sizeof(real));
	 return 1;
 }

#include "linear.h"
#include "lcp.h"
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
 
#define SKIP(n) (2*(n)+1)
 
 /*
  *
  */
static int pivot(real *M,int row,int col,int n)
{
	int i,skip;
	skip = SKIP(n);
	real d = M[row*skip+col];
	printf("pivot[%d,%d]\n",row,col);
	printM(M,NULL,n,skip);	
	if( d == 0 )return 0;
	multiply_line(M,1.0/d,row,n,skip);
	M[row*skip+col] = 1;
	for(i=0;i<n;i++){
		if( i != row&&M[i*skip+col]!=0 )
			elimination(M,i,row,col,n,skip);
	}
	return 1;
}

/*
 *
 */
static int check_get_result_and_free_base(real * M,int *base,real *x,int n)
{
	int i,row,result,skip;
	skip = SKIP(n);
	memset(x,0,2*n*sizeof(real));
	result = 1;
	for(i=0;i<skip;i++){
		row = i%n;
		if(base[i]){
			real d = M[row*skip+skip-1];
			if(i<n){
				x[i+n] = d;
			}else{
				x[i-n] = d;
			}
			if( d < 0 )
				result = 0;
		}
	}
	free(base);
	return result;
}

/*
 *
 */
static int last_negative(real * M,int n,int skip)
{
	for(int i=n-1;i>=0;i--){
		if(M[i*skip+skip-1] < 0 ){
			printM(M,NULL,n,skip);
			printf("last_negative row = %d\n",i);			
			return i;
		}
	}
	return -1;
}

/*
 *
 */
static int sovle_principalPivot(real * M,real * x,int n,int m)
{
	int i,row,prev,skip;
	int *base = (int *)malloc(2*n*sizeof(int));
	skip = SKIP(n);
	for(i=0;i<n;i++){
		base[i] = 1;
		base[n+i] = 0;
	}
	prev = -1;
	while( (row=last_negative(M,n,skip)) != prev && row!=-1 ){
		prev = row;
		if( M[row*skip+row] != 1){
			base[row] = 1;
			base[row+n] = 0;
			if(!pivot(M,row,row,n,m))
				break;
		}else if( M[row*skip+row+n] != 1 ){
			base[row] = 0;
			base[row+n] = 1;
			if(!pivot(M,row,row+n,n,m))
				break;
		}		
	}

	printf("final result:\n");
	printM(M,NULL,n,skip);
	return check_get_result_and_free_base(M,base,x,n);
}

int lcp_pivot( real *A,real *b,real *x,int n)
{
	return lcp_pivot2(A,b,x,n,0);
}

/*
 * 算法可行条件是A是P-matrix
 *(主子式都大于0，主子式principal minor，子式minor是选择一个元素i,j
 * 去掉i行和j列剩下的矩阵的行列式的值，主子式是i=j)
 * 表示在A矩阵下面还有m行，在进行消元处理的时候参与变换
 * 而b元素下面还有m个元素,也在消元处理的时候参与变换
 * 将m引入主要配接mlcp的求解
 */
int lcp_pivot2( real *A,real *b,real *x,int n,int m)
{
	int i,j,k,skip;
	skip = SKIP(n);
	real * M = (real *)malloc((n+m)*skip*sizeof(real));
	/* 
	 * 构造一个 [I,-A,b] 增广矩阵
	 */
	for(i=0;i<n;i++){
		k = i*skip;
		for(j=0;j<n;j++){
			if(i!=j)
				M[k+j] = 0;
			else
				M[k+j] = 1;
			M[k+n+j] = -A[i*n+j];
		}
		M[k+skip-1] = b[i];
	}
	/*
	 * 将附加行复制到M中去
	 */
	for(i=0;i<m;i++){
		k = (i+n)*skip;
		for(j=0;j<n;j++){
			if(j<n)
				M[k+j] = 0;
			else
				M[k+n+j] = -A[(n+i)*n+j];
		}
		M[k+skip-1] = b[n+i];
	}
	k = sovle_principalPivot(M,x,n,m);
	free(M);
	return k;	
}
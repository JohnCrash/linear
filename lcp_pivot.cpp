 #include "linear.h"
#include "lcp.h"
 #include <stdlib.h>
 #include <memory.h>
 #include <stdio.h>
 
#define SKIP(n) (2*(n)+1)
 
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

static int sovle_principalPivot(real * M,real * x,int n)
{
	int i,row,skip;
	bool infea;
	int *base = (int *)malloc(2*n*sizeof(int));
	skip = SKIP(n);
	for(i=0;i<n;i++){
		base[i] = 1;
		base[n+i] = 0;
	}
	do{
		infea = false;
		row = -1;
		for(i=n-1;i>=0;i--){
			if( M[i*skip+skip-1] < 0 ){
				row = i;
				printM(M,NULL,n,skip);
				printf("principal row = %d\n",row);
				break;
			}
		}
		if( row >= 0 ){
			if( M[row*skip+row] != 1){
				base[row] = 1;
				base[row+n] = 0;
				if(!pivot(M,row,row,n))
					break;
				infea = true;
			}else if( M[row*skip+row+n] != 1 ){
				base[row] = 0;
				base[row+n] = 1;
				if(!pivot(M,row,row+n,n))
					break;
				infea = true;
			}		
		}
	}while(infea);

	printM(M,NULL,n,skip);
	return check_get_result_and_free_base(M,base,x,n);
}

int lcp_pivot( real *A,real *b,real *x,int n)
{
	int i,j,k,skip;
	skip = SKIP(n);
	real * M = (real *)malloc(n*skip*sizeof(real));
	/* 
	 * 构造一个 [I,-M,b] 增广矩阵
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
	k = sovle_principalPivot(M,x,n);
	free(M);
	return k;	
}
#include "lcp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "misc.h"

void mlcp_sor1(real *A,real *b,real *x,real *hi,real *lo,real w,int n,int k_max)
{
	int i,j,k,m,skip;
	skip = n;
	
	for(k=0;k<k_max;k++){
		for(i=0;i<n;i++){
			real delta = 0;
			m = i*skip;
			for(j=0;j<i;j++){
				delta += A[m+j]*x[j];
			}
			for(j=i+1;j<n;j++){
				delta += A[m+j]*x[j];
			}
			delta = (b[i]-delta)/A[m+i];
			x[i] += w*(delta-x[i]);
			if(x[i]>hi[i])x[i] = hi[i];
			if(x[i]<lo[i])x[i] = lo[i];
		}
	}
}
 #include "linear.h"
 #include "lcp.h"
 #include <stdlib.h>
 #include <memory.h>
 #include <stdio.h>
 
 #define SKIP(n) (2*(n)+2)
 
 static void printM(real * M,int * N,int n)
 {
	 int skip = SKIP(n);
	 printf("-----------------------------------------------------------------\n");
	 for(int i=0;i<n;i++){
		 if(N){
			 if(N[i]!=-1)
				printf("[%d]",N[i]);
			else
				printf("[X]");
		 }
		for(int j=0;j<skip;j++){
			real d = M[i*skip+j];
			if(d>0)
				printf(" %.2f ",d);
			else if(d<0)
				printf("%.2f ",d);
			else
				printf(" 0.00 ");
		 }
		 printf("\n");
	 }
 }
 
static void sub_line(real * M,int i,int j,int n)
{
	 int k,skip = SKIP(n);
	 int src,des;
	 des = i*skip;
	 src = j*skip;
	 for(k=0;k<skip;k++){
		 M[des+k] -= M[src+k];
	 }
}
 
static void add_line(real * M,int i,int j,int n)
{
	 int k,skip = SKIP(n);
	 int src,des;
	 des = i*skip;
	 src = j*skip;
	 for(k=0;k<skip;k++){
		 M[des+k] += M[src+k];
	 }
}
 
static void multiply_line(real * M,real d,int i,int n)
 {
	 int k,skip = SKIP(n);
	 int des;
	 des = i*skip;
	 for(k=0;k<skip;k++){
		 M[des+k] *= d;
	 }
 }
 
 static void negative_line(real * M,int j,int n)
 {
	 int i,skip = SKIP(n);
	 int k = j*skip;
	 for(i=0;i<skip;i++){
		 M[k+i] = -M[k+i];
	 }
 }
 
 static void elimination(real *M,int eli,int row,int col,int n)
 {
	 int i,skip = SKIP(n);
	 real d = M[eli*skip+col];
	 for(i=0;i<skip;i++){
		 M[eli*skip+i] -= d*M[row*skip+i];
	 }
 }
 
static void pivot(real *M,int *N,int n,int row,int col)
{
	int skip = SKIP(n);
	real d = 1.0/M[row*skip+col];
	N[row] = col;
	printf("pivot[%d,%d]\n",row,col);
	printM(M,N,n);	
	multiply_line(M,d,row,n);
	for(int k=0;k<n;k++){
		if(k!=row&&M[k*skip+col]!=0){
			elimination(M,k,row,col,n);
		}
	}
}

 static int init_probrem(real * M,int *N,int n,int *enter)
 {
	 int i,j;
	 real d,r;
	 int skip = SKIP(n);
	 d = 0;
	 j = -1;
	 for(i=0;i<n;i++){
		 r = M[i*skip+skip-1];
		 if( r < d ){
			 d = r;
			 j = i;
		 }
	 }
	 if(j!=-1){
		 *enter = N[j];
		 pivot(M,N,n,j,2*n);
		 return 1;
	 }
	 return 0;
 }

static int argmin_element(real * M,int *N,int n,int* enter,int * prow,int *pcol)
{
	int i,j,k,skip = SKIP(n);
	real ratios,d,r;
	ratios = FLT_MAX;
	j = -1;
	k = *enter;
	printf("argmin_element enter=%d\n",k);
	printM(M,N,n);
	if(k==2*n)return 0;
	for(i=0;i<n;i++){
		d = M[i*skip+n+k];
		if(d>0){
			r = M[i*skip+skip-1]/d;
			if(r<ratios){
				j = i;
				ratios = r;
			}
		}
	}
	if(j!=-1){
		*prow = j;
		*pcol = n+k;
		*enter = N[j];
		printf("argmin_element return row=%d,col=%d,ratios=%.2f,enter=%d\n",
		j,k,ratios,N[j]);
		return 1;
	}else return 0;
}

/*
 *
 */
static int check_get_result_and_free_N(real *M,int * N,real * x,int n)
{
	int skip = SKIP(n);
	real d;
	int result = 1;
	memset(x,0,2*n*sizeof(real));
	
	for(int i=0;i<n;i++){
		d = M[i*skip+skip-1];
		if(N[i]<n){
			x[n+N[i]] = d; //获取y
		}else if(N[i]<2*n){
			x[N[i]-n] = d; //获取x
		}else{
			result = 0; //未将辅助变量2*n消去，lemke算法没有成功完成
		}
		if(d<0)result = 0; //q中还存在小于0的值，算法没有成功
	}
	free(N);
	return result;
} 

int solve_lemke(real * M,real * x,int n)
{
	int i,enter,skip;
	int * N = (int *)malloc(n*sizeof(int));
	skip = SKIP(n);
	for(i=0;i<n;i++)N[i]=i;
	printM(M,N,n);
	printf("init_probrem\n");
	if(init_probrem(M,N,n,&enter)){
		int row,col;
		printM(M,N,n);
		while( argmin_element(M,N,n,&enter,&row,&col) ){
			pivot(M,N,n,row,col);
		}
	}	 
	return check_get_result_and_free_N(M,N,x,n);
}
 
 /*
  *
  */
 int lcp_lemke(real * A,real *b,real *x,int n)
 {
	int i,j,k,skip;
	skip = SKIP(n);
	real * M = (real *)malloc(n*skip*sizeof(real));
	/* 
	 * 构造一个 [I,-M,-1,b] 增广矩阵
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
		M[k+skip-2] = -1;
		M[k+skip-1] = b[i];
	}
	k = solve_lemke(M,x,n);
	free(M);
	return k;
 }
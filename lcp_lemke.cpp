 #include "linear.h"
 #include "lcp.h"
 #include <stdlib.h>
 #include <memory.h>
 #include <stdio.h>
 
 #define SKIP(n) (2*(n)+1)
 
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
 
 static int positive_b(real * M,int n)
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
		 for(i=0;i<n;i++){
			 if(i!=j){
				 sub_line(M,i,j,n);
			 }
		 }
		 negative_line(M,j,n);
		 return 1;
	 }
	 return 0;
 }
 
 static bool isused(int *N,int j,int n)
 {
	 for(int i=0;i<n;i++){
		 if( N[i]!=-1 && N[i]%n==j%n )return true;
	 }
	 return false;
 }
 
 static int argmin_element(real * M,int *N,int n,int * prow,int *pcol)
 {
	 int i,j,k,s,r,skip = SKIP(n);
	 real line_max,ratios,d;
	 ratios = FLT_MAX;
	 k = s = -1;
	 for(i=0;i<n;i++){
		 if( N[i]==-1 ){
			 line_max = -FLT_MAX;
			 for(j=0;j<skip-1;j++){
				 if( M[i*skip+j] > line_max&&!isused(N,j,n) ){
					 r = j;
					 line_max = M[i*skip+j];
				 }
			 }
			 if(line_max>0){
				d = M[i*skip+skip-1]/line_max;
				if( d < ratios ){
					ratios = d;
					s = i;
					k = r;
				}
			 }
		 }
	 }
	 if(k!=-1&&s!=-1){
		 *prow=s;
		 *pcol = k;
		 return 1;
	 }else return 0;
 }
 
 int solve_lemke(real * M,real * x,int n)
 {
	 int i,j,k,skip;
	 real d;
	 int * N;
	 skip = SKIP(n);
	 printM(M,NULL,n);
	 printf("positive_b\n");
	 if(!positive_b(M,n)){
		 for(i=0;i<n;i++){
			 x[i] = 0;
			 x[n+i] = M[i*skip+skip-1];
		 }
		 return 1;
	 }
	 N = (int *)malloc(n*sizeof(int));
	 for(i=0;i<n;i++)N[i]=-1;
	 i = 0;
	 while(i<n){
		 int row,col;
		 printM(M,N,n);
		 if( argmin_element(M,N,n,&row,&col) ){
			N[row] = col;
			d = 1.0/M[row*skip+col];
			printf("argmin_element[%d,%d]\n",row,col);
			multiply_line(M,d,row,n);
			printM(M,N,n);
			for(k=0;k<n;k++){
				if(k!=row&&M[k*skip+col]!=0){
					elimination(M,k,row,col,n);
				}
			}
			i++;
		 }else{
			 printf("argmin stop loop\n");
			 goto failexit;
		 }
	 }
	 for(i=0;i<n;i++){
		 if(N[i]==-1||M[i*skip+skip-1]<0){
			 printf("q<0 stop loop\n");
			 goto failexit;
		 }
	 }
	 memset(x,0,2*n*sizeof(real));
	 for(i=0;i<n;i++){
		x[N[i]] = M[i*skip+skip-1];
	 }
	 return 1;
failexit:
	 free(N);
	 return 0;
 }
 
 /*
  * 1,2,3,4
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
		M[k+skip-1] = b[i];
	}
	k = solve_lemke(M,x,n);
	free(M);
	return k;
 }
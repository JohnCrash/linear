 #include "linear.h"
 #include "lcp.h"
 
static void sub_line(real * M,int i,int j,int n)
{
	 int k,skip = 2*n+2;
	 int src,des;
	 des = i*skip;
	 src = j*skip;
	 for(k=0;k<skip;k++){
		 M[des+k] -= M[src+k];
	 }
}
 
static void add_line(real * M,int i,int j,int n)
{
	 int k,skip = 2*n+2;
	 int src,des;
	 des = i*skip;
	 src = j*skip;
	 for(k=0;k<skip;k++){
		 M[des+k] += M[src+k];
	 }
}
 
static void multiply_line(real * M,real d,int i,int n)
 {
	 int k,skip = 2*n+2;
	 int des;
	 des = i*skip;
	 for(k=0;k<skip;k++){
		 M[des+k] *= d;
	 }
 }
 
 static void negative_line(real * M,int j,int n)
 {
	 int i,skip = 2*n+2;
	 int k = j*skip;
	 for(i=0;i<skip;i++){
		 M[k+i] = -M[k+i];
	 }
 }
 
 static int positive_b(real * M,int n)
 {
	 int i,j;
	 real d,r;
	 int skip = 2*n+2;
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
 
 static int argmin_element(real * M,int n,int * prow,int *pcol)
 {
	 int i,j,k,s,skip = 2*n+2;
	 real line_max,ratios,d;

	 ratios = FLT_MAX
	 for(i=0;i<n;i++){
		 line_max = -FLT_MAX;
		 for(j=0;j<skip-1;j++){
			 if( M[i*skip+j] > line_max ){
				 k = j;
				 line_max = M[i*skip+j];
			 }
		 }
		 if(line_max!=0){
			d = M[i*skip+skip-1]/line_max;
			if( d < ratios ){
				ratios = d;
				s = i;
			}
		 }else{
			return 0;
		 }
	 }
	 *prow=s;
	 *col = k;
	 return 1;
 }
 
 int solve_lemke(real * M,real * x,int n)
 {
	 int i,j,skip;
	 real d;
	 skip = 2*n+2;
	 if(!positive_b(M,n)){
		 for(i=0;i<n;i++){
			 x[i] = 0;
			 x[n+i] = M[i*skip+skip-1+i];
		 }
		 return 1;
	 }
	 
	 if( argmin_element(M,n,&i,&j) ){
		 d = M[i*skip+j];
		 if( d != 0 )
			multiply_line(M,1/d,i,n);
	 }
	 return 0;
 }
 
 /*
  * 1,2,3,4
  */
 int lcp_lemke(real * A,real *b,real *x,int n)
 {
	int i,j,k,skip;
	skip = 2*n+2;
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
		M[k+2*n] = -1;
		M[k+2*n+1] = b;
	}
	k = solve_lemke(M,x,n);
	free(M);
	return k;
 }
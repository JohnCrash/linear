 #include "linear.h"
 #include "lcp.h"
 #include <stdlib.h>
 #include <memory.h>
 #include <stdio.h>
 
void printM(real * M,int * N,int n,int skip)
 {
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
 
void sub_line(real * M,int i,int j,int n,int skip)
{
	 int k;
	 int src,des;
	 des = i*skip;
	 src = j*skip;
	 for(k=0;k<skip;k++){
		 M[des+k] -= M[src+k];
	 }
}
 
void add_line(real * M,int i,int j,int n,int skip)
{
	 int k;
	 int src,des;
	 des = i*skip;
	 src = j*skip;
	 for(k=0;k<skip;k++){
		 M[des+k] += M[src+k];
	 }
}
 
void multiply_line(real * M,real d,int i,int n,int skip)
 {
	 int k;
	 int des;
	 des = i*skip;
	 for(k=0;k<skip;k++){
		 M[des+k] *= d;
	 }
 }
 
 void negative_line(real * M,int j,int n,int skip)
 {
	 int i;
	 int k = j*skip;
	 for(i=0;i<skip;i++){
		 M[k+i] = -M[k+i];
	 }
 }
 
 /*
  * row,col为主元，(eli,col)是要消为0的元素。
  */
 void elimination(real *M,int eli,int row,int col,int n,int skip)
 {
	 int i;
	 real d = M[eli*skip+col];
	 for(i=0;i<skip;i++){
		 M[eli*skip+i] -= d*M[row*skip+i];
	 }
 }
 
 /*
  * pivot将row,col对应的元素作为主元，将col列的其他元素消为0
  * 同时将该行的基置为col(N[row] = col)
  * 条件pivot前(row,col)对应的元素不能为0
  */
static void pivot(real *M,int *N,int n,int row,int col,int m,int skip)
{
	int k;
	real d = 1.0/M[row*skip+col];
	N[row] = col;
	printf("pivot[%d,%d]\n",row,col);
	printM(M,N,n,skip);	
	multiply_line(M,d,row,n,skip);
	for(k=0;k<n;k++){
		if(k!=row&&M[k*skip+col]!=0){
			elimination(M,k,row,col,n,skip);
		}
	}
	if(m>0){
		for(k=0;k<m;k++){
			if(k!=row&&M[(n+k)*skip+col]!=0){
				elimination(M,n+k,row,col,n,skip);
			}			
		}
	}
}

/*
 * 初始化步骤
 * 1.如果b都>=0则直接返回0，solve_lemke会调用check_get_result_and_free_N获取解
 * 2.在b列中找到最负的数
 */
 static int init_probrem(real * M,int *N,int n,int *enter,int skip)
 {
	 int i,j;
	 real d,r;
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

 /*
  * enter是进基列索引，将在进基列中搜索ratios最小的行。
  * 然后确定下一个进基列。
  */
static int argmin_element(real * M,int *N,int n,int* enter,int * prow,int *pcol,int skip)
{
	int i,j,k,skip = SKIP(n);
	real ratios,d,r;
	ratios = FLT_MAX;
	j = -1;
	k = *enter;
	printf("argmin_element enter=%d\n",k);
	printM(M,N,n,skip);
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
 * 检查结果如果有小于0的则失败，如果基列中还存在辅助列也失败。
 * 如果成功根据基变量将结果取出。同时释放N
 */
static int check_get_result_and_free_N(real *M,int * N,real * x,int n,int skip)
{
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


int lcp_lemkeBlock(real *M,real *x,int n,int m,int skip)
{
	int i,enter,skip;
	int * N = (int *)malloc(n*sizeof(int));
	skip = SKIP(n);
	for(i=0;i<n;i++)N[i]=i;
	printM(M,N,n,skip);
	printf("init_probrem\n");
	if(init_probrem(M,N,n,&enter,skip)){
		int row,col;
		printM(M,N,n,skip);
		while( argmin_element(M,N,n,&enter,&row,&col,skip) ){
			pivot(M,N,n,row,col,m,skip);
		}
	}	 
	return check_get_result_and_free_N(M,N,x,n,skip);	
}

 /*
  * 使用lemke法求解线性互补问题
  * y = Ax+b
  * yx'>=0,x>=0,y>=0 其中x'表示x的转置,A是一个nxn矩阵
  * 构造一个增广矩阵
  * [ I][y]	
  * [-A][x] = b 	
  * 如果矩阵A是下面的一种，lemke可以用于求解LCP
  * A是copositive plus matrix,
  * 
  */
 int lcp_lemke(real * A,real *b,real *x,int n)
 {
	int skip;
	real * M = mallocLemkeMatrix(A,b,n,0,&skip);
	int result = lcp_lemkeBlock(M,x,n,0,skip);
	free(M);
	return result;
 }
 
 real *mallocLemkeMatrix(real * A,real *b,int n,int m,int *pskip)
 {
	int i,j,k,skip;
	skip = 2*n+2;
	*pskip = skip;
	real * M = (real *)malloc(n*skip*sizeof(real));
	/* 
	 * 构造一个 [I,-M,-1,b] 增广矩阵,-1是辅助变量
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
	return M;
 }
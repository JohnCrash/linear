#include "linear.h"
#include "lcp.h"
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>

/*
 * 将方阵src的第i行变负复制到des的第j行
 */
static void copyColumAndNeg(real *des,real *src,int j,int i,int n)
{
	for(int k=0;k<n;k++){
		des[k*n+j] = -src[k*n+i];
	}
}

/*
 * 将单位矩阵的第i行复制到des的第j行
 */
static void copyIdentity(real *des,int j,int i,int n)
{
	for(int k=0;k<n;k++){
		if(k==i)
			des[k*n+j] = 1;
		else
			des[k*n+j] = 0;
	}
}

static void printMat(const char * s,real * A,int N)
{
	int i, j;
	if(s)
		printf("%s\n",s);
	for (i = 0; i<N; i++){
		for (j = 0; j<N; j++){
			printf("%.4f\t", A[i*N + j]);
		}
		printf("\n");
	}
}

/*
 * N是一个互补集合，其中N[i]=0表示对于的x[i]=0 或者N[i]=1表示对应的x[i]>=0
 */
static int doSolveN(real *A,real *b,int *N,real * M,real *P,real *L,real *x,int nub,int n)
{
	int i,j;
	for(i=0,j=0;i<n;i++,j++){
		if(N[j]){
			copyColumAndNeg(M,A,j,i,n);
		}else{
			copyIdentity(M,j,i,n);
		}
	}
	if( crout_plu(M,P,L,n) ){
		if( solve_plu(P,L,M,b,x,n) ){
			for(i=0;i<n;i++){
				if(i>=nub&&x[i]<0){
					return 0;
				}
				x[n+i] = N[i];
			}
			return 1;
		}
	}
	return 0;
}

static void printN(int *N,int n)
{
	for(int i=0;i<n;i++){
		printf("%d ",N[i]);
	}
	printf("\n");
}

/*
 * 将线性互补问题全部解穷举出来（因此算法复杂度为O(pow(2,n)*pow(n,3))）
 * y = Ax+b;x,y>=0;x'y=0
 * 成功返回1并将全部的解都放入到xs中,失败返回0
 * 注意返回的解x,y是混合放置的，后面接互补集合，互补集合为1对应x,0对应y
 * 如a b c d 1 0 1 1,那么x=a 0 c d,y=0 c 0 0
 *
 * nub 是前nub个x变量是自由变量，而y=0
 * 其他部分和lcp函数相同，仅仅是前nub个变量不参与互补条件判断
 */
int mlcp(real * A,real *b,std::vector<real *>& xs,int nub,int n)
{
	if(nub<0||nub>n)return 0;
	int i;
	int * N = (int *)malloc(n*sizeof(int));
	real * M = (real *)malloc(n*n*sizeof(real));
	real * P = (real *)malloc(n*n*sizeof(real));
	real * L = (real *)malloc(n*n*sizeof(real));
	real * bb = (real *)malloc(n*sizeof(real));
	real * x = NULL;
	memset(N,0,n);
	for(i=0;i<nub;i++)
		N[i]=1;
	do{
		if(!x)
			x = (real*)malloc(2*n*sizeof(real));
		memcpy(bb,b,n*sizeof(real));
		/*
		 * 求解已知互补条件N的方程,M,P,L是中间变量防止多次分配
		 * x返回解,如果有解并且x>=0,y>=0返回1
		 */
		if( doSolveN(A,bb,N,M,P,L,x,nub,n) ){
			xs.push_back(x);
			x = NULL;
		}
		/* 
		 * 这里产生一个互补集N,覆盖全部可能性，每一个位上都可能为0或者1
		 * 这个正好和二进制的进一样，模拟二进制进位。
		 * 确保小于nub的都为0
		 */
		for(i=nub;i<n;i++){
			if( N[i]==0 ){
				N[i] = 1;
				break;
			}else
				N[i] = 0;
		}
	}while(i<n);
	free(N);
	free(M);
	free(P);
	free(L);
	free(bb);
	return 0;
}

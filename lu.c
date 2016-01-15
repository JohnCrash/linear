#include "linear.h"

void zero(real * A,int n)
{
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			A[i][j]=0;
		}
	}
}

void identity(real * A,int n)
{
	zero(A);
	for(int i=0;i<n;i++)
		A[i][i] = 1;
}

//交互i和j行上的全部元素
static void xchangeRaw(real * A,int n,int i,int j)
{
	int temp;
	for(int x=0;x<n;x++){
		temp = A[i][x];
		A[i][x] = A[j][x];
		A[j][x] = temp;
	}
}

//查找raw,col下面(同一列)全部元素中绝对值最大的行
static int absMaxLeading(real * A,int n,int raw,int col)
{
	int maxcol=-1;
	real mv=0;
	for(int i=raw;i<n;i++){
		int v = fabs(A[i][col]);
		if(v>mv){
			maxcol=i;
			mv=v;
		}
	}
	return maxcol;
}

//消元，将低i行对应元素乘以d加到j上。
static void decol(real * A,int n,int i,int j,real d)
{
	for(int k=0;k<n;k++)
		A[j][k] += A[i][k]*d;
}

/*
 * 进行lu分解，将U存入A中，将L存入到L中
 *
 */
int lu(real * A,real * L,int n)
{
	int m,i,j,k;
	real mr,v,d;
	identity(L,n);
	for(i=0;i<n-1;i++){ //columns
		m=absMaxLeading(A,n,i,i);
		if(m!=i)
			xchangeRaw(A,n,i,m);
		mr = A[i][i];
		if(mr==0)
			return i;
		for(j=i+1;j<n;j++){ //raws
			v = A[j][i];
			if(v!=0){
				d = -v/mr;
				decol(A,n,j,i,v);
				A[j][i]=0; 
			}
			
		}
	}
}
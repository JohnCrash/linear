#include "linear.h"
#include "lcp.h"
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>

/*
 * Principad pivot method
 * 问题：w = A*z+q , w'z=0 , z>=0 , w>=0
 * 将方程从下下层下面的
 *       [w]
 * [I,-M][ ]=[q] ,其中I是单位矩阵.
 *       [z]
 * 算法举例：
 * 基变量 | w1 w2 w3 z1 z2 z3 | q
 * -------+-------------------+---
 *	w1	  |  1  0  0 -1  0  0 |-1
 *	w2	  |  0  1  0 -2 -1  0 |-1
 *	w3	  |  0  0  1 -2 -2[-1]|-1 <--从下向上第一个负数
 * -------+-------------------+---
 *	w1	  |  1  0  0 -1  0  0 |-1
 *	w2	  |  0  1  0 -2[-1] 0 |-1 <--继续从下向上找到的第一个负数
 *	z3	  |  0  0 -1  2  2  1 | 1 <--交换基向量，进行pivot
 * -------+-------------------+---
 *	w1	  |  1  0  0 -1  0  0 |-1
 *	z2	  |  0 -1  0  2  1  0 | 1
 *	z3	  |  0  2[-1]-2  0  1 |-1
 * -------+-------------------+---
 *	w1	  |  1  0  0[-1] 0  0 |-1
 *	z2	  |  0 -1  0  2  1  0 | 1
 *	w3	  |  0 -2  1  2  0 -1 | 1
 * -------+-------------------+---
 *	z1	  | -1  0  0  1  0  0 | 1
 *	z2	  |  2 -1  0  0  1  0 |-1
 *	w3	  |  2 -2  1  0  0[-1]|-1
 * -------+-------------------+---
 *	z1	  | -1  0  0  1  0  0 | 1
 *	z2	  |  2[-1] 0  0  1  0 |-1
 *	z3	  | -2  2 -1  0  0  1 | 1
 * -------+-------------------+--- 
 *	z1	  | -1  0  0  1  0  0 | 1
 *	w2	  | -2  1  0  0 -1  0 | 1
 *	z3	  |  2  0[-1] 0  2  1 |-1
 * -------+-------------------+---
 *	z1	  | -1  0  0  1  0  0 | 1
 *	w2	  | -2  1  0  0 -1  0 | 1
 *	w3	  | -2  0  1  0 -2 -1 | 1
 * -------+-------------------+---
 *              ^  ^  ^
 *              |  |  |
 *    +---------+  |  |
 *    |  +---------+  |
 *    |  |  +---------+
 *    |  |  |
 *(0 ,1 ,1 ,1 ,0 ,0 )
 *(w1,w2,w3,z1,z2,z3)
 */ 
 /*
  * 将row,col位置变成1，列的其他元素消元为0
  * 运算作用于附加区
  */
static int pivot(real *M,int row,int col,int n,int m,int skip)
{
	int i;
	real d = M[row*skip+col];
	//printf("pivot[%d,%d]\n",row,col);
	//printM(M,NULL,n+m,skip);	
	if( d == 0 )return 0;
	multiply_line(M,1.0/d,row,n,skip);
	M[row*skip+col] = 1;
	for(i=0;i<n;i++){
		if( i != row&&M[i*skip+col]!=0 )
			elimination(M,i,row,col,n,skip);
	}
	/* 对附加区进行变换，应用于mlcp求解程序 */
	for(i=0;i<m;i++){
		if( M[(n+i)*skip+col]!=0 )
			elimination(M,n+i,row,col,n,skip);
	}	
	return 1;
}

/*
 * 根据基向量从b中取出x,y
 * base是一个基向量标识数组
 * |   y   |   x   |
 * |<--n-->|<--n-->|
 * |0 1 0 1|1 0 1 0|
 *    ^
 *    |
 *   表示在对应位置的b(M的最后一列)中取出y2
 */
static int check_get_result_and_free_base(real * M,int *base,real *x,int n,int skip)
{
	int i,row,result;
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
 * 从下向上查找b行的第一个负的元素，如果发现就返回该行的行数。否则返回-1
 */
static int last_negative(real * M,int n,int skip)
{
	for(int i=n-1;i>=0;i--){
		if(M[i*skip+skip-1] < 0 ){
			//printM(M,NULL,n,skip);
			//printf("last_negative row = %d\n",i);			
			return i;
		}
	}
	return -1;
}

/*
 * 矩阵M的结构是
 * |	y 	|	x	|b|
 * |1	   0|		|
 * |   1   0|		|
 * |	   1|		|
 * |下面接m行		|
 * principad pivoting methiod不停的交换y区的和x区的主元，直到b都为正
 * 或者算法无法继续推进（返回的列和上次的一样）。
 * x 是一个2n向量，前面的n个元素是解x,后面y。skip是M的行距。
 * 成功返回1，失败0.
 */
int lcp_pivotBlock(real *M,real *x,int n,int m,int skip)
{
	int i,row,prev,s;
	int *base = (int *)malloc(2*n*sizeof(int));
	/*
	 * base用来标识基向量
	 * |<--n-->|<--n-->|
	 * |   y   |   x   |
	 * |1111111|0000000|
	 * 初始化的时候将基向量都设置为y区
	 */
	for(i=0;i<n;i++){
		base[i] = 1;
		base[n+i] = 0;
	}
	prev = -1;
	s = 0;
	/*
	 *FIXBUG: 算法会陷入死循环，目前简单的加一个简单限制s++<2n
	 */
	while( (row=last_negative(M,n,skip)) != prev && row!=-1 && s++<=2*n){
		prev = row;
		/*
		 * 算法会确保基向量总是等于1,在x区和y区选择一个基向量
		 */
		if( M[row*skip+row] != 1){
			base[row] = 1;
			base[row+n] = 0;
			if(!pivot(M,row,row,n,m,skip))
				break;
		}else if( M[row*skip+row+n] != 1 ){
			base[row] = 0;
			base[row+n] = 1;
			if(!pivot(M,row,row+n,n,m,skip))
				break;
		}
	}

	//printf("final result:\n");
	//printM(M,NULL,n+m,skip);
	return check_get_result_and_free_base(M,base,x,n,skip);
}

/*
 * 算法可行条件是A是P-matrix
 *(主子式都大于0，主子式principal minor，子式minor是选择一个元素i,j
 * 去掉i行和j列剩下的矩阵的行列式的值，主子式是i=j)
 * 表示在A矩阵下面还有m行，在进行消元处理的时候参与变换
 * 而b元素下面还有m个元素,也在消元处理的时候参与变换
 * 将m引入主要配接mlcp的求解
 */
int lcp_pivot( real *A,real *b,real *x,int n)
{
	int skip;
	real * M = mallocPivotMatrix( A,b,n,0,&skip );
	int result = lcp_pivotBlock(M,x,n,0,skip);
	free(M);
	return result;	
}

/*
 * 构造一个求解矩阵M，A是一个nxn方阵，m附加在A和b下面的附加数据(参与pivot运算)
 * 附加行数由m制定，没有设置为0，附加的行包括b也要附加m个元素。
 * 通过pskip返回构造的矩阵的行距(一行有多少个real)
 * 成功返回一个构造好的矩阵。
 */
real *mallocPivotMatrix(real * A,real *b,int n,int m,int *pskip)
{
	int i,j,k,skip;
	skip = 2*n+1;
	*pskip = skip;
	real * M = (real *)malloc((n+m)*skip*sizeof(real));
	/* 
	 * 构造一个 [I,-A,b] 增广矩阵
	 */
	for(i=0;i<(n+m);i++){
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
	return M;
}
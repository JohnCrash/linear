 #include "linear.h"
 #include "lcp.h"
 #include <stdlib.h>
 #include <memory.h>
 #include <stdio.h>
 
 /*
  * 对A的左上角(nub)x(nub)的区域进行消元处理，将该区块处理成单元矩阵
  */
 static int pivotNub(real *A,int *P,int nub,int n,int skip)
 {
	int i,j;
	for(i=0;i<nub;i++){
		if( A[i*skip+i] == 0 )return 0; //FIXBUG,使用交换行来解决
		real d = 1.0/A[i*skip+i];
		multiply_line(A,d,i,n,skip);
		A[i*skip+i] = 1;
		for(j=0;j<n;j++){
			if(j!=i&&A[j*skip+i]!=0){
				elimination(A,j,i,i,n,skip);
			}
		}
	}
	return 1;
 }
 
/*
 * M表示如下(M|y,x|'=0,'表示转置)，nub为自由区。
 * 因为在nub区x=R(任意实数)y=0，因此y的nub区被排除
 * x,y的base区中选择n-nub个基列和x的nub区组成方程得到解
 * 也就是y base区和x base区互补选择n-nub个列与x nub区组成最终的求解矩阵
 * 先将M弄成下面的样式
 * |	y	 |	 x	  |d|b|
 * |nub	|base|nub|base|d|b|
 * |	|0	0|1 0|	  |
 * |	|0	0|0 1|	  |
 * |	|1  0|0 0|	  |
 * |	| 1  |	 |	  |
 * |	|  1 |	 |	  |
 * |	|0  1|0 0|	  |
 * x base 和 y base 区互补，y nub区因为y=0而被屏蔽了,由x base和y base挑选
 * 出来的基和x nub 区组成nxn方阵。
 * 当把M转化为这样以后，y的base区和x的base区转化为一个标准线性互补问题
 */
int mlcpSolver(real * A,real *b,real *x,int nub,int n,LCPSolver solver)
{
	int i,j,skip;
	int * P;
	int result = 0;
	
	if(nub<0||nub>n)return result;
	P = (int *)malloc(nub*sizeof(int));
	/*
	 * nub = n问题完全转化为线性方程问题
	 * nub = 0问题是一个完全的线性互补问题
	 */
	for(i = 0;i<nub;i++)
		P[i] = i;
	real * M = (real*)malloc(n*(n+1)*sizeof(real));
	skip = n+1;
	
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			M[i*skip+j] = -A[i*n+j];
		}
		M[i*skip+skip-1] = b[i];
	}
	if( pivotNub(M,P,nub,n,skip) ){
		if(nub<n){
			/*
			 *|<---n+1--->|
			 *|	  x	   | b|
			 *|nub|base|
			 *|1 0|	   |
			 *|0 1|	   |
			 *|0 0|	   |
			 *|0 0|	   |
			 *如果将x nub左上角的对角区叫做A,x base右上角叫做B，等等..
			 *|nub|base| b|
			 *| A | B  | e|
			 *| C | D  | E|
			 */
			printf("pivotNub\n");
			printM(M,NULL,n,skip);
			
			real *AA = (real *)malloc((n-nub)*n*sizeof(real));
			real * bb = (real *)malloc(n*sizeof(real));
			real * xx = (real *)malloc(2*(n-nub)*sizeof(real));
			real * MM = NULL;
			int skip2 = n-nub;
			for(i=nub;i<n;i++){
				for(j=nub;j<n;j++){
					AA[(i-nub)*skip2+j-nub] = -M[i*skip+j];
				}
				bb[i-nub] = M[i*skip+skip-1];
			}
			for(i=0;i<nub;i++){
				for(j=0;j<(n-nub);j++){
					AA[(i+n-nub)*skip2+j] = -M[i*skip+j+nub];
				}
				bb[i+n-nub] = M[i*skip+skip-1];
			}
			printf("new matrix AA\n");
			printM(AA,NULL,n,skip2);
			printf("new matrix bb\n");
			printM(bb,NULL,n,1);			
			/* 
			 *上面的操作,将AA组合成下面的结构AA
			 *| -D | n-nub行
			 *| -B | nub行
			 *而bb是一个下面的结构
			 *| E |
			 *| e |
			 */
			if(solver==LEMKE){
				MM = mallocLemkeMatrix(AA,bb,n-nub,nub,&skip2);
				printf("mallocLemkeMatrix MM\n");
				printM(MM,NULL,n,skip2);
				result = lcp_lemkeBlock(MM,xx,n-nub,nub,skip2);
			}else if(solver==PIVOT){
				MM = mallocPivotMatrix(AA,bb,n-nub,nub,&skip2);
				printf("mallocPivotMatrix MM\n");
				printM(MM,NULL,n,skip2);				
				result = lcp_pivotBlock(MM,xx,n-nub,nub,skip2);			
			}
			/*
			 * 通过上面的lcp求解,xx是n-nub对互补解,skip2是矩阵MM的行距
			 */
			printf("lcp solve\n");
			printM(MM,NULL,n,skip2);
			/*
			 * 收集解
			 */
			for(i=0;i<nub;i++){
				x[i] = MM[(i+n-nub)*skip2+skip2-1]; //nub区x
				x[n+i] = 0; //nub区y
			}
			for(i=0;i<n-nub;i++){
				x[nub+i] = xx[i]; //互补区x
				x[n+nub+i] = xx[n-nub+i]; //互补区y
			}
			free(MM);
			free(xx);
			free(bb);
			free(AA);
		}else{
			/*
			 * 存的线性方程求解,下面收集解
			 */
			for(i=0;i<nub;i++){
				x[i] = M[i*skip+skip-1]; //nub区x
				x[n+i] = 0;
			}
		}
	}
	free(M);
	free(P);
	return result;
}
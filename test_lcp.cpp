#include "lcp.h"
#include "misc.h"

bool check_lcp_result(real * A,real *b,real *x,int n)
{
	int i;
	real * y = (real*)malloc(n*sizeof(real));
	multiply0(y,A,x,n,n,1);
	for(i =0;i<n;i++){
		y[i] += b[i];
	}
	printf("ceck lcp result:\nx=");
	for(i=0;i<n;i++){
		printf("%.4f\t",x[i]);
	}
	printf("\ny=");
	for(i=0;i<n;i++){
		printf("%.4f\t",y[i]);
	}
	printf("\n");
	free(y);
}

static void test_pgs()
{
	/*
	 * 2,4
	 */
	srand(4);
	real * A = makeRandMatrix2();
	real * b = makeRandVec2();
	real * x = makeRandVec2();
	real * AA = makeMatrix();
	real * bb = makeRandVec2();
	copyMatrix(AA,A);
	memcpy(bb,b,N*sizeof(real));
	memset(x,0,sizeof(real)*N);

	printMat("solve lcp A=",A);
	printVec("lcp b=",b);
	printVec("lcp x=",x);
	int result = lcp_pgs(A,b,x,N);
	
	printf("lcp solve %s (%d)\n",result?"true":"false",result);
	printVec("lcp solve ",x);
	
	printMat("solve lcp AA=",AA);
	printVec("lcp bb=",bb);
	printVec("lcp x=",x);

	check_lcp_result(AA,bb,x,N);
	freeMatrix(A);
	freeMatrix(b);
	freeMatrix(x);
	freeMatrix(AA);
	freeMatrix(bb);
}

int main(int argn,char * argv[])
{
	test_pgs();
	return 0;
}
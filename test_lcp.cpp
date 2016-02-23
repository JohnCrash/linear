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

void printLCPVx(std::vector<real *>& vx)
{
	for(auto i=vx.begin();i!=vx.end();++i){
		printf("[%d]\t",i-vx.begin());
		for(int j =0;j<2*N;j++){
			if(j==N)
				printf("|");
			printf("%.2f ",(*i)[j]);
		}
		printf("\n");
	}	
}

static void test_pgs()
{
	real * A = makeRandSPDMatrix();
	real * b = makeRandVec2();
	real * x = makeRandVec2();
	real * xx = makeRandVec2();
	real * AA = makeMatrix();
	real * bb = makeRandVec2();
	std::vector<real *> vx;
	copyMatrix(AA,A);
	memcpy(bb,b,N*sizeof(real));
	//memset(x,0,sizeof(real)*N);
	memcpy(x,b,N*sizeof(real));
	//memcpy(xx,b,N*sizeof(real));
	memset(xx,0,sizeof(real)*N);
	
	printMat("solve lcp A=",A);
	printVec("lcp b=",b);
	printVec("lcp x=",x);
	int result = lcp(A,b,vx,N);	
	printLCPVx(vx);
	result = lcp_pgs(A,b,x,N,15,0.01);
	//int result = Solve_GaussSeidel(A,b,x,N,15);
	//for(int i =0;i<N;i++)
	//	b[i] = -b[i];
	//solve_gauss_seidel(A,b,xx,N,15);
	printf("lcp solve %s (%d)\n",result?"true":"false",result);
	printVec("lcp solve ",x);
	printVec("gs solve ",xx);
	printMat("solve lcp AA=",AA);
	printVec("lcp bb=",bb);
	printVec("lcp x=",x);

	check_lcp_result(AA,bb,x,N);
	freeMatrix(A);
	freeMatrix(b);
	freeMatrix(x);
	freeMatrix(AA);
	freeMatrix(bb);
	freeLcpSolve(vx);
}

int main(int argn,char * argv[])
{
	if(argn>=1){
		/*
		 * 2,4
		 */
		 int s = atoi(argv[1]);
		 printf("random seed : %d\n",s);
		srand(s);
	}
		
	test_pgs();
	return 0;
}
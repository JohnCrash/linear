#include "lcp.h"
#include "misc.h"

void check_lcp_result(real * A,real *b,real *x,int n)
{
	int i;
	bool p = true;
	real * y = (real*)malloc(n*sizeof(real));
	multiply0(y,A,x,n,n,1);
	for(i =0;i<n;i++){
		y[i] += b[i];
		if( abs(y[i]*x[i])>0.01 )
			p = false;
	}
	printf("ceck lcp_pgs result: %s \nx=",p?"pass":"fail");
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
	memset(x,0,N*sizeof(real));
	//memcpy(xx,b,N*sizeof(real));
	memset(xx,0,sizeof(real)*N);
	
	printf("--------------------------------------------------------\n");
	printMat("solve lcp A=",A);
	printVec("lcp b=",b);
	printf("--------------------------------------------------------\n");
	int result1 = lcp(A,b,vx,N);	
	int result2 = lcp_pgs(A,b,x,N,15,0.001);
	int result3 = Solve_GaussSeidel(AA,bb,xx,N,15);
	printf("lcp solve:\n");
	printf("--------------------------------------------------------\n");
	printLCPVx(vx);
	printf("lcp_pgs solve %s (%d)\n",result2?"true":"false",result2);
	printVec("lcp_pgs=",x);
	check_lcp_result(AA,bb,x,N);
	printVec("gs solve :",xx);

	freeMatrix(A);
	freeMatrix(b);
	freeMatrix(x);
	freeMatrix(AA);
	freeMatrix(bb);
	freeLcpSolve(vx);
}

static real exaples_m[]= {
	1,-1,-1,-1,
	-1,1,-1,-1,
	1,1,2,0,
	1,1,0,2
};

static real exaples_q[] = {
	3,5,-9,-5,
};

static void test_lemke()
{
	real * A = makeRandSPDMatrix();
	real * b = makeRandVec2();
	real * x = (real *)malloc(2*N*sizeof(real));
	real * xx = (real *)malloc(2*N*sizeof(real));
	real * AA = makeMatrix();
	real * bb = makeRandVec2();
	std::vector<real *> vx;
	copyMatrix(A,exaples_m);
	for(int i =0;i<N;i++)
		x[i] = exaples_q[i];
	copyMatrix(AA,A);
	memcpy(bb,b,N*sizeof(real));
	//memset(x,0,sizeof(real)*N);
	memset(x,0,N*sizeof(real));
	//memcpy(xx,b,N*sizeof(real));
	memset(xx,0,sizeof(real)*N);	
	
	printf("--------------------------------------------------------\n");
	printMat("solve lcp A=",A);
	printVec("lcp b=",b);
	printf("--------------------------------------------------------\n");
	int result1 = lcp(A,b,vx,N);	
	int result2 = lcp_lemke(A,b,x,N);
	
	printf("lcp solve:\n");
	printf("--------------------------------------------------------\n");
	printLCPVx(vx);
	printf("lcp_lemke solve %s (%d)\n",result2?"true":"false",result2);
	printVec("lcp_lemke=",x);
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
		
	//test_pgs();
	test_lemke();
	return 0;
}
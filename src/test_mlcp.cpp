#include "lcp.h"
#include "misc.h"
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>

void initMLCProblem(int n,real **PA,real **pb,real **px,real **py)
{
	*PA = (real *)malloc(n*n*sizeof(real));
	*pb = (real *)malloc(n*sizeof(real));
	*px = (real *)malloc(2*n*sizeof(real));	
	*py = (real *)malloc(2*n*sizeof(real));
}

void freeMLCPProblem(real *A,real *b,real *x,real *y)
{
	free(A);
	free(b);
	free(x);
	free(y);
}

/*
 * 构造一个随机的混合互补问题
 */
void MLCProblem(int n,real *A,real *b,real *x,int *pnub)
{
	*pnub = (int)(randomReal()*(n+1));
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			A[i*n+j] = (randomReal()>0.5?1:-1)*randomReal();
		}
		b[i] = 2*(randomReal()>0.5?1:-1)*randomReal();
		x[i] = 0;
		x[n+i] = 0;
	}
}

/*
 * 校验MLCP的解
 */
int verifyMLCP(real *A,real *b,real *x,real *y,int nub,int n)
{
	int i;
	
	multiply0(y,A,x,n,n,1);
	for(i=0;i<nub;i++){
		if(!FTEQ(y[i],0))return 0;
	}
	for(i=nub;i<n;i++){
		if(x[i]<0||y[i]<0)return 0;
		if(!FTEQ(x[i]*y[i],0))return 0;
	}
	return 1;
}

static void printNumber(real a)
{
	if( a > 0 )
		printf(" %3.2f ",a);
	else
		printf("%3.2f ",a);
}
/*
 *打印MLCP问题和解
 */
void printMLCP(real *A,real *b,real *x,int nub,int n)
{
	int i,j;
	for(i=0;i<n;i++){
		if(i<nub)
			printf("*");
		else
			printf(" ");
		for(j=0;j<n;j++){
			printNumber(A[i*n+j]);
		}
		printf("|");
		printNumber(b[i]);
		printf("\n");
	}
	printf("[");
	for(i=0;i<2*n;i++){
		printNumber(x[i]);
	}
	printf("]\n");
}

/*
 * 正确性检测
 */
int testSovler(LCPSolver sovler,int n,int reps)
{
	real *A,*b,*x,*y;
	int nub;
	int count = 0;
  	initMLCProblem(n,&A,&b,&x,&y);
	for(int i=0;i<reps;i++){
		MLCProblem(n,A,b,x,&nub);
		mlcpSolver(A,b,x,nub,n,sovler);
		if(verifyMLCP(A,b,x,y,nub,n)){
			//printMLCP(A,b,x,nub,n);
			count++;
		}
	}
	freeMLCPProblem(A,b,x,y);
	return count;
}

/*
 * 使用遍历解，看看到底有多少有解
 */
int calcSolveCount(int n,int reps)
{
	real *A,*b,*x,*y;
	int nub;
	int count = 0;
	std::vector<real *> xs;
  	initMLCProblem(n,&A,&b,&x,&y);
	for(int i=0;i<reps;i++){
		MLCProblem(n,A,b,x,&nub);
		mlcp(A,b,xs,nub,n);
		if(!xs.empty() && verifyMLCP(A,b,xs[0],y,nub,n) ){
			//printMLCP(A,b,x,nub,n);
			count++;
		}
		freeLcpSolve(xs);
	}
	freeMLCPProblem(A,b,x,y);
	return count;
} 
/*
 * 检查算法死循环的问题
 */
int testSovlerCheck(LCPSolver sovler,int n,int reps)
{
	real *A,*b,*x,*y;
	int i,nub;
	int count = 0;
	bool isexit = false;
	int tc = 0;
	std::thread * pthread = new std::thread(
		[&](){
			while(!isexit){
				while(tc<5000){
					std::this_thread::sleep_for(std::chrono::microseconds(10));
					tc++;
					if(isexit)return;
				}
				printf("%d'th\n",i);
				printMLCP(A,b,x,nub,n);
				tc = 0;
			}
		}
	);
  	initMLCProblem(n,&A,&b,&x,&y);
	for(i=0;i<reps;i++){
		MLCProblem(n,A,b,x,&nub);
		tc = 0;
		mlcpSolver(A,b,x,nub,n,sovler);
		if(verifyMLCP(A,b,x,y,nub,n)){
			//printMLCP(A,b,x,nub,n);
			count++;
		}
	}
	freeMLCPProblem(A,b,x,y);
	isexit = true;
	pthread->join();
	delete pthread;	
	return count;
}

int main(int argn,char * argv[])
{
	int s;
	if(argn>=1){
		/*
		 * 2,4
		 */
		s = atoi(argv[1]);
		printf("random seed : %d\n",s);
		srand(s);
	}
	for( int i=0;i<1;i++ ){
		int n = i*10+3;
		//srand(s);
		//printf("PGS %d - %d/100\n",n,testSovler(PGS,n));
		srand(s);
		printf("ALL %d - %d/1000\n",n,calcSolveCount(n,1000));
		srand(s);
		printf("LEMKE %d - %d/1000\n",n,testSovler(LEMKE,n,1000));
		srand(s);
		printf("PIVOT %d - %d/1000\n",n,testSovler(PIVOT,n,1000));
	}
	//test_mlcp_pgs();
	//test_mlcpSolver();
	//test_sor();
	return 0;
}
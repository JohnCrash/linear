#include "lcp.h"
#include "misc.h"
#include <vector>
#include <algorithm>
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
	for(i=0;i<n;i++){
		y[i] += b[i];
	}
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
	//
	multiply0(x+n,A,x,n,n,1);
	for(i=0;i<n;i++){
		x[n+i] += b[i];
	}
	printf("[");
	for(i=0;i<2*n;i++){
		printNumber(x[i]);
	}
	printf("]\n");	
}

std::vector<int> _sovleInx;
/*
 * 正确性检测
 */
int testSovler(LCPSolver sovler,int n,int reps,int &solveC,int &verifyC)
{
	real *A,*b,*x,*y;
	int nub;
	int count = 0;
	verifyC = solveC = 0;
  	initMLCProblem(n,&A,&b,&x,&y);
	for(int i=0;i<reps;i++){
		MLCProblem(n,A,b,x,&nub);
		if(mlcpSolver(A,b,x,nub,n,sovler)){
			solveC++;
			count++;
		}
		if(verifyMLCP(A,b,x,y,nub,n)){
			//printMLCP(A,b,x,nub,n);
			verifyC++;
			//if(std::find(_sovleInx.begin(),_sovleInx.end(),i)!=_sovleInx.end()){
			//	printf("Sovler : %s %d\n",sovler==LEMKE?"LEMKE":"PIVOT",i);
			//	printMLCP(A,b,x,nub,n);
			//}
		}
	}
	freeMLCPProblem(A,b,x,y);
	return count;
}

/*
 * 使用遍历解，看看到底有多少有解
 */
int calcSolveCount(int n,int reps,int &solveC,int &verifyC)
{
	real *A,*b,*x,*y;
	int nub,k;
	int count = 0;
	verifyC = solveC = 0;
	std::vector<real *> xs;
  	initMLCProblem(n,&A,&b,&x,&y);
	for(int i=0;i<reps;i++){
		MLCProblem(n,A,b,x,&nub);
		mlcp(A,b,xs,nub,n);
		if(!xs.empty()){
			//printMLCP(A,b,x,nub,n);
			if( verifyMLCP(A,b,xs[0],y,nub,n) )
				verifyC++;
			else{
				printf("mlcp:\n");
				printMLCP(A,b,xs[0],nub,n);
			}
			_sovleInx.push_back(i);
			solveC++;
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

#define S 100

int main(int argn,char * argv[])
{
	int s;
	int n = 3;
	if(argn>1){
		/*
		 * 2,4
		 */
		s = atoi(argv[1]);
		printf("random seed : %d\n",s);
		srand(s);
		if(argn>2)
			n = atoi(argv[2]);
	}
	for( int i=0;i<1;i++ ){
		
		//srand(s);
		//printf("PGS %d - %d/100\n",n,testSovler(PGS,n));
		srand(s);
		int solveC,verifyC;
		double t0 = getClock();
		int c = calcSolveCount(n,S,solveC,verifyC);
		printf("ALL [%d] - %d-%d-%d/%d\n",n,solveC,verifyC,c,S);
		printf("use time : %fms\n",1000*(getClock()-t0));
		srand(s);
		t0 = getClock();
		c = testSovler(LEMKE,n,S,solveC,verifyC);
		printf("LEMKE [%d] - %d-%d-%d/%d\n",n,solveC,verifyC,c,S);
		printf("use time : %fms\n",1000*(getClock()-t0));
		srand(s);
		t0 = getClock();
		c = testSovler(PIVOT,n,S,solveC,verifyC);
		printf("PIVOT [%d] - %d-%d-%d/%d\n",n,solveC,verifyC,c,S);
		printf("use time : %fms\n",1000*(getClock()-t0));
	}
	//test_mlcp_pgs();
	//test_mlcpSolver();
	//test_sor();
	return 0;
}
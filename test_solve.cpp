#include "linear.h"
#include "misc.h"

void printColumn(real *c)
{
	for(int i = 0;i<N;i++){
		printf("[%d][%.2f]\n",i,c[i]);
	}
}

void test_solve_plu()
{
	real * A = makeRandMatrix();
	real * P = makeMatrix();
	real * L = makeMatrix(); 
	real * U = makeMatrix();
	real * b,*x,*bb;
	b = (real*)malloc(N*sizeof(real));
	x = (real*)malloc(N*sizeof(real));
	bb = (real*)malloc(N*sizeof(real));
	copyMatrix(U,A);
	for(int i = 0;i<N;i++){
		b[i] = (int)(randomReal()*20);
		bb[i] = b[i];
	}
	int ret = lu(U,P,L,N);
	if(ret){
		printMat("A:",A);
		printColumn(b);
		printMat3(P,L,U,N);
		ret = solve_plu(P,L,U,b,x,N);
	}
	if(ret){
		printf("solve:\n");
		printColumn(x);
		multiply0(b,A,x,N,N,1);
		printDiffent1("result:\n",b,bb);
	}
	freeMatrix(A);
	freeMatrix(P);
	freeMatrix(L);
	freeMatrix(U);
	free(b);
	free(x);
}

void test_solve_crout_plu()
{
	real * A = makeRandMatrix();
	real * P = makeMatrix();
	real * L = makeMatrix(); 
	real * U = makeMatrix();
	real * b,*x,*bb;
	b = (real*)malloc(N*sizeof(real));
	x = (real*)malloc(N*sizeof(real));
	bb = (real*)malloc(N*sizeof(real));
	copyMatrix(U,A);
	for(int i = 0;i<N;i++){
		b[i] = (int)(randomReal()*20);
		bb[i] = b[i];
	}
	int ret = crout_plu(U,P,L,N);
	if(ret){
		printMat("A:",A);
		printColumn(b);
		printMat3(P,L,U,N);
		ret = solve_plu(P,L,U,b,x,N);
	}
	if(ret){
		printf("solve:\n");
		printColumn(x);
		multiply0(b,A,x,N,N,1);
		printDiffent1("result:\n",b,bb);
	}
	freeMatrix(A);
	freeMatrix(P);
	freeMatrix(L);
	freeMatrix(U);
	free(b);
	free(x);
}

int main(int argn, char *argv[])
{
	readom_init();
	test_solve_plu();
	test_solve_crout_plu();
	return 0;
}
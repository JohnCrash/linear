#include "linear.h"
#include "misc.h"

void printColumn(real *c)
{
	for(int i = 0;i<NN;i++){
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
	b = (real*)malloc(NN*sizeof(real));
	x = (real*)malloc(NN*sizeof(real));
	bb = (real*)malloc(NN*sizeof(real));
	copyMatrix(U,A);
	for(int i = 0;i<NN;i++){
		b[i] = (int)(randomReal()*20);
		bb[i] = b[i];
	}
	int ret = lu(U,P,L,NN);
	if(ret){
		printMat("A:",A);
		printColumn(b);
		printMat3(P,L,U,NN);
		ret = solve_plu(P,L,U,b,x,NN);
	}
	if(ret){
		printf("solve:\n");
		printColumn(x);
		multiply0(b,A,x,NN,NN,1);
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
	b = (real*)malloc(NN*sizeof(real));
	x = (real*)malloc(NN*sizeof(real));
	bb = (real*)malloc(NN*sizeof(real));
	copyMatrix(U,A);
	for(int i = 0;i<NN;i++){
		b[i] = (int)(randomReal()*20);
		bb[i] = b[i];
	}
	int ret = crout_plu(U,P,L,NN);
	if(ret){
		printMat("A:",A);
		printColumn(b);
		printMat3(P,L,U,NN);
		ret = solve_plu(P,L,U,b,x,NN);
	}
	if(ret){
		printf("solve:\n");
		printColumn(x);
		multiply0(b,A,x,NN,NN,1);
		printDiffent1("result:\n",b,bb);
	}
	freeMatrix(A);
	freeMatrix(P);
	freeMatrix(L);
	freeMatrix(U);
	free(b);
	free(x);
}

void printM(real *A,int n)
{
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			printf("%.4f\t",A[i*n+j]);
		}
		printf("\n");
	}	
}

void printC(real *b,int n)
{
	for(int j=0;j<n;j++){
		printf("%.4f\t",b[j]);
	}
	printf("\n");	
}

/*
 * 数值例子来至于wiki
 * https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method
 * https://en.wikipedia.org/wiki/Jacobi_method#Description
 */
void test_gauss_seidel()
{
	const real aa[]={
		10,-1,2,0,
		-1,11,-1,3,
		2,-1,10,-1,
		0,3,-1,8
	};
	const real bb[]={
		6,25,-11,15
	};
	real * A = (real*)malloc(4*4*sizeof(real));
	real * b = (real*)malloc(4*sizeof(real));
	real * x = (real*)malloc(4*sizeof(real));
	memcpy(A,aa,4*4*sizeof(real));
	memcpy(b,bb,4*sizeof(real));
	memset(x,0,4*sizeof(real));
	printf("Gauss-Seidel solve Ax=b\nA=\n");
	printM(A,4);
	printf("=======================\nb=\n");
	printC(b,4);
	printf("=======================\nx=\n");
	printC(x,4);	
	for(int j=0;j<5;j++){
		solve_gauss_seidel(A,b,x,4,1);
		//solve_jacobi(A,b,x,4,1);
		//printM(A,4);
		printf("=======================\nx=\n");
		printC(x,4);
	}
	free(b);
	free(A);
}

int main(int argn, char *argv[])
{
	readom_init();
	test_solve_plu();
	test_solve_crout_plu();
	test_gauss_seidel();
	return 0;
}
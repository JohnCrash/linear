#include "linear.h"
#include "misc.h"

void test_inverse0()
{
	real * A = makeRandMatrix();
	real * B = makeMatrix();
	real * C = makeMatrix();
	real * X = makeMatrix();

	identity(X, NN);
	printf("test inverse0\n");
	copyMatrix(B, A);
	printMat("A=", A);
	inverse(A,C, NN);
	printMat("inverse=", C);
	multiply0(A, B, C, NN, NN, NN);
	printMat("A*A'=", A);
	printDiffent("A*A'=I", X, A);
	freeMatrix(A);
	freeMatrix(B);
	freeMatrix(C);
	freeMatrix(X);
}

void test_lup()
{
	real * A = makeRandMatrix();
	real * B = makeMatrix();
	real * L = makeMatrix();
	real * X = makeMatrix();
	real * P = makeMatrix();
	
	identity(X, NN);
	printf("test crout_plu\n");
	copyMatrix(B, A);
	printMat("A=", A);
	crout_plu(A,P,L,NN);
	multiply0(X, P, L, NN, NN, NN);
	multiply0(P, X, A, NN, NN, NN);
	printMat("A*A'=", P);
	printDiffent("A*A'=I", B, P);
	freeMatrix(A);
	freeMatrix(B);
	freeMatrix(L);
	freeMatrix(X);
	freeMatrix(P);
}

int main(int argn, char *argv[])
{
	readom_init();
	test_inverse0();
	test_lup();
	return 0;
}
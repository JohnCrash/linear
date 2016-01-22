#include "linear.h"
#include "misc.h"

void test_inverse0()
{
	real * A = makeRandMatrix();
	real * B = makeMatrix();
	real * C = makeMatrix();
	real * X = makeMatrix();

	identity(X, N);
	printf("test inverse0\n");
	copyMatrix(B, A);
	printMat("A=", A);
	inverse(A,C, N);
	printMat("inverse=", C);
	multiply0(A, B, C, N, N, N);
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
	
	identity(X, N);
	printf("test crout_plu\n");
	copyMatrix(B, A);
	printMat("A=", A);
	crout_plu(A,P,L,N);
	multiply0(X, P, L, N, N, N);
	multiply0(P, X, A, N, N, N);
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
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
	inverse0(A,C, N);
	printMat("inverse=", C);
	multiply0(A, B, C, N, N, N);
	printMat("A*A'=", A);
	printDiffent("A*A'=I", X, A);
	freeMatrix(A);
	freeMatrix(B);
	freeMatrix(C);
	freeMatrix(X);
}

int main(int argn, char *argv[])
{
	test_inverse0();
	return 0;
}
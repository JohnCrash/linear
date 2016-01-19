#include "linear.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

int N = 5;
real S = 100;

void printMat(const char * s,real * A)
{
	int i, j;
	if(s)
		printf("%s\n",s);
	for (i = 0; i<N; i++){
		for (j = 0; j<N; j++){
			printf("%f ", A[i*N + j]);
		}
		printf("\n");
	}
}

static void copyMatrix(real * des,real * src)
{
	memcpy(des,src,N*N*sizeof(real));
}

static real * makeMatrix()
{
	real * A = (real*)malloc(N*N*sizeof(real));
	return A;	
}

static real * makeRandMatrix()
{
	real * A = makeMatrix();
	random_matrix(A,N);
	return A;
}

static void freeMatrix(real * A)
{
	free(A);
}

static void test_lu_1()
{
	real * A = makeRandMatrix();
	real * L = makeMatrix(); 
	real * C = makeMatrix();
	printf("test lu\n");
	printMat("A=",A);
	lu(A, L, N);
	printMat("U=",A);
	printMat("L=",L);
	multiply0(C, L, A, N, N, N);
	printMat("L*U=",C);
	freeMatrix(A);
	freeMatrix(L);
	freeMatrix(C);
}

static void test_pldu_1()
{
	real * A = makeRandMatrix();
	real * L = makeMatrix(); 
	real * C = makeMatrix();
	real * P = makeMatrix();
	real * D = makeMatrix();
	real * T = makeMatrix();
	printf("test pldu\n");
	printMat("A=",A);
	pldu(A, P,D,L, N);
	printMat("P=",P);
	printMat("D=",D);
	printMat("U=",A);
	printMat("L=",L);
	multiply0(C, P, L, N, N, N);
	multiply0(T, C, D, N, N, N);
	multiply0(C, T, A, N, N, N);
	printMat("C=P*L*D*U",C);
	freeMatrix(A);
	freeMatrix(L);
	freeMatrix(C);
	freeMatrix(P);
	freeMatrix(D);
	freeMatrix(T);	
}

static void test_inverse_low_triangle()
{
	real * A = makeRandMatrix();
	real * B = makeMatrix();
	real * C = makeMatrix();
	printf("test inverse_low_triangle\n");
	clearUpperTriangle(A,N);
	printMat("A=",A);
	copyMatrix(B,A);
	inverse_low_triangle(A,N);
	printMat("inverse=",A);
	multiply0(C, B, A, N, N, N);
	printMat("A*A'=",C);	
	freeMatrix(A);
	freeMatrix(B);
	freeMatrix(C);	
}

static void test_inverse_upper_triangle()
{
	real * A = makeRandMatrix();
	real * B = makeMatrix();
	real * C = makeMatrix();	
	printf("test inverse_upper_triangle\n");
	clearLowerTriangle(A,N);
	printMat("A=",A);
	copyMatrix(B,A);
	inverse_upper_triangle(A,N);
	printMat("inverse=",A);
	multiply0(C, B, A, N, N, N);
	printMat("A*A'=",C);	
	freeMatrix(A);
	freeMatrix(B);
	freeMatrix(C);		
}

static void test_inverse_pivoting()
{
	real * A = makeRandMatrix();
	real * B = makeMatrix();
	real * C = makeMatrix();
	printf("test inverse_pivoting\n");
	identity(A, N);
	for(int i=0;i<N;i++){
		int m = (int)(N*(real)rand()/(real)RAND_MAX);
		xchangeRaw(A,N,i,m);
	}
	copyMatrix(B, A);
	printMat("A=",A);
	inverse_pivoting(A, 3);
	printMat("inverse=",A);
	multiply0(C, B, A, 3, 3, 3);
	printMat("A*A'=",C);
	freeMatrix(A);
	freeMatrix(B);
	freeMatrix(C);
}

static void test_inverse_diagonal()
{
	real * A = makeRandMatrix();
	real * B = makeMatrix();
	real * C = makeMatrix();	

	printf("test inverse_diagonal\n");
	clearLowerTriangle(A,N);
	clearUpperTriangle(A,N);
	copyMatrix(B,A);
	printMat("A=",A);
	inverse_diagonal(A, N);
	printMat("inverse=",A);
	multiply0(C, B, A, N, N, N);
	printMat("A*A'=",C);
	freeMatrix(A);
	freeMatrix(B);
	freeMatrix(C);		
}

static void test_inverse_1()
{
	real * A = makeRandMatrix();
	real * L = makeMatrix(); 
	real * C = makeMatrix();
	real * P = makeMatrix();
	real * D = makeMatrix();
	real * T = makeMatrix();

	printf("test inverse_1\n");
	copyMatrix(T,A);
	printMat("A=",A);
	pldu(A, P, D, L, N);
	printMat("P=",P);
	printMat("D=",D);
	printMat("U=",A);
	printMat("L=",L);
	inverse(P, L, D, A, N);
	printMat("inverse=",P);
	multiply0(C, T, P, N, N, N);
	printMat("A*A'=",C);
	freeMatrix(A);
	freeMatrix(L);
	freeMatrix(C);
	freeMatrix(P);
	freeMatrix(D);
	freeMatrix(T);		
}

static void test_crout_lu()
{
	real * A = makeRandMatrix();
	real * L = makeMatrix();
	real * C = makeMatrix();

	printf("test crout_lu\n");
	printMat("A=",A);
	crout_lu(A, L, N);
	printMat("U=",A);
	printMat("L=",L);
	multiply0(C, L, A, N, N, N);
	printMat("C=L*U",C);	
	freeMatrix(A);
	freeMatrix(L);
	freeMatrix(C);
}

int main(int argn,const char *argv[])
{
	readom_init();
	test_lu_1();
	test_pldu_1();
	test_inverse_low_triangle();
	test_inverse_upper_triangle();
	test_inverse_pivoting();
	test_inverse_diagonal();
	test_inverse_1();
	test_crout_lu();
	return 0;
}


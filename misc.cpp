#include "linear.h"
#include "misc.h"

void printDiffent(const char * s,real * A, real *B)
{
	real v = 0;
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			real s = fabs(A[i*N+j]-B[i*N+j]);
			if (s>v){
				v = s;
			}
		}
	}
	printf("%s --> %s\n", s, v>0.001 ? "failed" : "passed");
}

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

void copyMatrix(real * des,real * src)
{
	memcpy(des,src,N*N*sizeof(real));
}

real * makeMatrix()
{
	real * A = (real*)malloc(N*N*sizeof(real));
	return A;	
}

real * makeRandMatrix()
{
	real * A = makeMatrix();
	random_matrix(A,N);
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			A[i*N + j] = (int)(A[i*N + j]*10);
		}
	}
	return A;
}

void freeMatrix(real * A)
{
	free(A);
}
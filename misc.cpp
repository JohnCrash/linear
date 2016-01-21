#include "linear.h"
#include "misc.h"

int _disablePrint = 0;

void disablePrint(int b)
{
	_disablePrint = b;
}

int printDiffent(const char * s,real * A, real *B)
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
	if(!_disablePrint)
	printf("%s --> %s\n", s, v>FTACC ? "failed" : "passed");
	return v>FTACC?0:1;
}

void printMat(const char * s,real * A)
{
	if(!_disablePrint){
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

#if defined(WIN32)||defined(_WIN32)
#include <windows.h>
double getClock()
{
	double t;
	t = (double)GetTickCount();
	return t/1000.0;
}
#else
#include <sys/time.h>
double getClock()
{
	double t;
	timeval tv;
	gettimeofday(&tv,NULL);
	t  = (double)tv.tv_sec+(double)(tv.tv_usec)/1000000.0;
	return t;
}	
#endif
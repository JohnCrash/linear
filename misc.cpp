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

//打印三个并列的矩阵[P][A][B]
void printMat3(real * P,real * A,real *B,int n)
{
	if(_disablePrint)
		return;
	real * X;
	real * S;
	X = (real *)malloc(n*n*sizeof(real));
	S = (real *)malloc(n*n*sizeof(real));
	multiply0(X,P,A,n,n,n);
	multiply0(S,X,B,n,n,n);

	for(int i=0;i<n;i++){
		printf("[%d][",i);
		for(int j=0;j<n;j++){
			printf("%.2f ",P[i*n+j]);
		}		
		printf("][");
		for(int j=0;j<n;j++){
			printf("%.2f ",A[i*n+j]);
		}
		printf("][");
		for(int j=0;j<n;j++){
			printf("%.2f ",B[i*n+j]);
		}
		printf("][");
		for(int j=0;j<n;j++){
			printf("%.2f ",S[i*n+j]);
		}		
		printf("]\n");
	}
	free(X);
	free(S);
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
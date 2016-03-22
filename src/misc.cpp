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
	for (int i = 0; i < NN; i++){
		for (int j = 0; j < NN; j++){
			real s = fabs(A[i*NN+j]-B[i*NN+j]);
			if (s>v){
				v = s;
			}
		}
	}
	if(!_disablePrint)
		printf("%s --> %s\n", s, v>FTACC ? "failed" : "passed");
	return v>FTACC?0:1;
}

int printDiffent1(const char * s,real * b,real *x)
{
	real v = 0;
	for(int i=0;i<NN;i++){
		real s = fabs(b[i]-x[i]);
		if(s>v){
			v=s;
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
		for (i = 0; i<NN; i++){
			for (j = 0; j<NN; j++){
				printf("%.4f\t", A[i*NN + j]);
			}
			printf("\n");
		}
	}
}

void printVec(const char * s,real *v)
{
	if(_disablePrint)
		return;
	if(s)
		printf("%s\n",s);	
	for(int i=0;i<NN;i++){
		printf("%.4f\t",v[i]);
	}
	printf("\n");
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
	memcpy(des,src,NN*NN*sizeof(real));
}

real * makeMatrix()
{
	real * A = (real*)malloc(NN*NN*sizeof(real));
	return A;	
}

real * makeRandMatrix()
{
	real * A = makeMatrix();
	random_matrix(A,NN);
	for (int i = 0; i < NN; i++){
		for (int j = 0; j < NN; j++){
			A[i*NN + j] = (real)((int)(A[i*NN + j]*10));
		}
	}
	return A;
}

real * makeRandSPDMatrix()
{
	real * B = makeMatrix();
	real * A = makeRandMatrix();
	random_matrix(A,NN);
	multiply0(B,A,A,NN,NN,NN);
	free(A);
	return B;
}

real * makeRandMatrixN(int n)
{
	real * A = makeMatrix();
	random_matrix(A,n);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			A[i*n + j] = (real)((int)(A[i*n + j]*10));
		}
	}
	return A;
}

real * makeRandSPDMatrixNUB(int nub)
{
	real * B = makeRandMatrix();
	real * A = makeRandMatrixN(NN-nub);
	real * C = makeRandMatrixN(NN-nub);
	random_matrix(A,NN-nub);
	multiply0(C,A,A,NN-nub,NN-nub,NN-nub);
	for(int i=nub;i<NN;i++){
		for(int j=nub;j<NN;j++){
			B[i*NN+j] = C[(i-nub)*NN+j-nub];
		}
	}
	free(A);
	free(C);
	return B;
}

real * makeRandVec()
{
	real * V = (real *)malloc(NN*sizeof(real));
	random_vector(V,NN);
	for(int i = 0; i < NN; i++){
		V[i] = (real)((int)(V[i]*10));
	}
	return V;
}

real randNegative()
{
	return (real)(randomReal()>0.5?1.0:-1.0);
}

real * makeRandMatrix2()
{
	real * A = makeMatrix();
	random_matrix(A,NN);
	for (int i = 0; i < NN; i++){
		for (int j = 0; j < NN; j++){
			A[i*NN + j] = (int)(A[i*NN + j]*10)*randNegative();
		}
	}
	return A;
}

real * makeRandVec2()
{
	real * V = (real *)malloc(NN*sizeof(real));
	random_vector(V,NN);
	for(int i = 0; i < NN; i++){
		V[i] = (int)(V[i]*10)*randNegative();
	}
	return V;
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
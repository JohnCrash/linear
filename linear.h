#ifndef _LINEAR_H_
#define _LINEAR_H_

#define dSINGLE 1

#if defined(dSINGLE)
#define REAL(x) (x##f)
typedef float real;
#else
#define REAL(x) (x)
typedef double real;	
#endif

#define fabs(x) ((x)>0?(x):-(x))

void zero(real * A,int n);
void identity(real * A,int n);

void readom_init();
void random_matrix(real * A,int n);

int lu(real * A,real * L,int n);
int pldu(real * A, real * P, real * D, real * L, int n);
int crout_lu(real * A,real * L,int n);
int crout_plu(real * A,real * P,real * L,int n);
void inverse(real * P, real * L, real * D, real * U, int n);

void multiplyC(real * A,real c,int n);

void multiply0(real *A, const real *B, const real *C, int p, int q, int r);
void multiply1(real *A, const real *B, const real *C, int p, int q, int r);
void multiply2(real *A, const real *B, const real *C, int p, int q, int r);

void transpose(real * A,int n);
void inverse_low_triangle(real * L, int n);
void inverse_upper_triangle(real * L, int n);
void inverse_pivoting(real * L, int n);
void inverse_diagonal(real * D, int n);

#endif
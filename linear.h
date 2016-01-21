#ifndef _LINEAR_H_
#define _LINEAR_H_

#define dSINGLE 1

#if defined(dSINGLE)
#define REAL(x) (x##f)
typedef float real;
#define FTACC 0.01f
#else
#define REAL(x) (x)
typedef double real;
#define FTACC 0.000001
#endif

#define fabs(x) ((x)>0?(x):-(x))
#define FTEQ(x,y) (fabs(x-y)>FTACC/100?0:1)

typedef long long int64;

void zero(real * A,int n);
void identity(real * A,int n);

void readom_init();
void random_matrix(real * A,int n);

int lu(real * A,real * P,real * L,int n);
int pldu(real * A, real * P, real * D, real * L, int n);
int crout_lu(real * A,real * L,int n);
int crout_plu(real * A,real * P,real * L,int n);
int inverse(real * P, real * L, real * D, real * U, int n);

void clearUpperTriangle(real * A, int n); 
void clearLowerTriangle(real * A,int n);
void multiplyC(real * A,real c,int n);

void xchangeRaw(real * A,int n,int i,int j);
int absMaxLeading(real * A,int n,int raw,int col);

void multiply0(real *A, const real *B, const real *C, int p, int q, int r);
void multiply1(real *A, const real *B, const real *C, int p, int q, int r);
void multiply2(real *A, const real *B, const real *C, int p, int q, int r);

void transpose(real * A,int n);
void inverse_low_triangle(real * L, int n);
void inverse_upper_triangle(real * L, int n);
void inverse_pivoting(real * L, int n);
void inverse_diagonal(real * D, int n);
int inverse0(real * A,real * B,int n);

#endif
#ifndef _LINEAR_H_
#define _LINEAR_H_

#define dSINGLE 1

#if defined(dSINGLE)
#define REAL(x) (x##f)
typedef float real;
#define FTACC 0.01f
	#ifndef FLT_MIN
	#define FLT_MIN 1.1754943508222875e-38
	#endif
	#ifndef FLT_MAX
	#define FLT_MAX 3.4028234663852886e+38
	#endif	
#else
#define REAL(x) (x)
typedef double real;
#define FTACC 0.000001
	#ifndef FLT_MIN
	#define FLT_MIN 2.2250738585072014e-308
	#endif
	#ifndef FLT_MAX
	#define FLT_MAX 1.7976931348623157e+308
	#endif	
#endif

#define fmax(a,b) ((a)>(b)?(a):(b))
#define fmin(a,b) ((a)<(b)?(a):(b))
#define fabs(x) ((x)>0?(x):-(x))
#define FTEQ(x,y) (fabs(x-y)>FTACC/100?0:1)

typedef long long int64;

/* 零矩阵 */
void zero(real * A,int n);

/* 单位矩阵 */
void identity(real * A,int n);

void readom_init();
void random_matrix(real * A,int n);
void random_vector(real * V,int n);
real randomReal();

/*
 * 进行lu分解，将U存入A中，将L存入到L中
 */
int lu(real * A,real * P,real * L,int n);

/*
 * crout algothim 
 * http://www.physics.utah.edu/~detar/phys6720/handouts/crout.txt
 * 将n*n矩阵A分解为LU,返回U被写入到A中
  * 成功返回1，失败返回0
 */
int crout_lu(real * A,real * L,int n);

/*
 * A=PLU,其中P为一个交换矩阵,返回U被写入到A中
 * 成功返回1，失败返回0
 */
int crout_plu(real * A,real * P,real * L,int n);

int cholesky(real * A,real *L,int n);

/*
 * 试验阶段
 * A的逆矩阵，A=P*L*D*U就是上面pldu分解的结果。
 * 最后将结果放入到P矩阵中。
 */
int inverse0(real * P, real * L, real * D, real * U, int n);

void clearUpperTriangle(real * A, int n); 
void clearLowerTriangle(real * A,int n);
void multiplyC(real * A,real c,int n);

void xchangeRaw(real * A,int n,int i,int j);
int absMaxLeading(real * A,int n,int raw,int col);

/* 
 * 矩阵乘法运算
 * multiply0 : A=B*C （B是一个p*q矩阵，C是q*r矩阵
 * multiply1 : A=B*C' C'是C的转置矩阵
 * multiply1 : A=B'*C
 */
void multiply0(real *A, const real *B, const real *C, int p, int q, int r);
void multiply1(real *A, const real *B, const real *C, int p, int q, int r);
void multiply2(real *A, const real *B, const real *C, int p, int q, int r);

/* 转置矩阵 */
void transpose(real * A,int n);
void inverse_low_triangle(real * L, int n);
void inverse_upper_triangle(real * L, int n);
void inverse_pivoting(real * L, int n);
void inverse_diagonal(real * D, int n);

/*
 *计算A的逆矩阵，将结果放入到B中 
 *将AB组成增广矩阵，将A部分通过初等变换转化为单位矩阵I。
 *这时候B即为A的逆矩阵 inverse该方法比inverse0速度快
*/
int inverse(real * A,real * B,int n);

/*
 * 解方程Ax=b,不能解返回0，成功返回1
 */
 int solve_lu(real * L,real * U,real * b,real *x,int n);
 int solve_plu(real * P,real * L,real * U,real * b,real *x,int n);
 int solve_low_triangle(real * L,real * b,real *x,int n);
 int solve_upper_triangle(real *L,real * b,real *x,int n);
 /*
  * 使用Gauss-Seidel迭代法解Ax=b的近似解
  * Gauss-Seidel迭代收敛条件是A为正定矩阵或者是主元占优矩阵
  */
 int solve_gauss_seidel(real * A,real * b,real * x,int n,int c);
 /*
  * 使用jacobi迭代法解Ax=b的近似解，收敛速度比Gauss-Seidel法慢
  */
 int solve_jacobi(real * A,real * b,real * x,int n,int c);
 
#endif
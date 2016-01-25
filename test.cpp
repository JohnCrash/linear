#include "linear.h"
#include "misc.h"

static int test_lu_1()
{
	real * A = makeRandMatrix();
	real * P = makeMatrix();
	real * L = makeMatrix(); 
	real * C = makeMatrix();
	real * X = makeMatrix();
	
	copyMatrix(X, A);
	printMat("A=",A);
	int r = lu(A,P, L, N);
	if(!r && N<5 ){
		printf("lu failed\n");
		disablePrint(0);
		printMat3(X,L,A,N);
		disablePrint(1);
	}
	printMat("U=",A);
	printMat("L=",L);
	printMat("P=", P);
	multiply0(C, P, L, N, N, N);
	multiply0(P, C, A, N, N, N);
	printMat("L*U=",P);
	int ret = printDiffent("A=L*U",X,P);
	freeMatrix(A);
	freeMatrix(P);
	freeMatrix(L);
	freeMatrix(C);
	freeMatrix(X);
	return r?ret:-1;
}
/*
static int test_pldu_1()
{
	real * A = makeRandMatrix();
	real * L = makeMatrix(); 
	real * C = makeMatrix();
	real * P = makeMatrix();
	real * D = makeMatrix();
	real * T = makeMatrix();
	real * X = makeMatrix();

	copyMatrix(X, A);
	printMat("A=",A);
	int r = pldu(A, P,D,L, N);
	printMat("P=",P);
	printMat("D=",D);
	printMat("U=",A);
	printMat("L=",L);
	multiply0(C, P, L, N, N, N);
	multiply0(T, C, D, N, N, N);
	multiply0(C, T, A, N, N, N);
	printMat("P*L*D*U=",C);
	int ret = printDiffent("A=P*L*D*U", X, C);
	freeMatrix(A);
	freeMatrix(L);
	freeMatrix(C);
	freeMatrix(P);
	freeMatrix(D);
	freeMatrix(T);	
	freeMatrix(X);
	return r?ret:-1;
}
*/
static int test_inverse_low_triangle()
{
	real * A = makeRandMatrix();
	real * B = makeMatrix();
	real * C = makeMatrix();
	real * X = makeMatrix();

	identity(X, N);
	clearUpperTriangle(A,N);
	for (int i = 0; i < N; i++)
		A[i*N + i] = 1;
	printMat("A=",A);
	copyMatrix(B,A);
	inverse_low_triangle(A,N);
	printMat("inverse=",A);
	multiply0(C, A, B, N, N, N);
	printMat("A*A'=",C);	
	printMat("X=", X);
	int ret = printDiffent("A*A'=I",X, C);
	freeMatrix(A);
	freeMatrix(B);
	freeMatrix(C);	
	freeMatrix(X);
	return ret;
}

static int test_inverse_upper_triangle()
{
	real * A = makeRandMatrix();
	real * B = makeMatrix();
	real * C = makeMatrix();	
	real * X = makeMatrix();

	identity(X, N);
	clearLowerTriangle(A,N);
	for (int i = 0; i < N; i++)
		A[i*N + i] = 1;
	printMat("A=",A);
	copyMatrix(B,A);
	inverse_upper_triangle(A,N);
	printMat("inverse=",A);
	multiply0(C, B, A, N, N, N);
	printMat("A*A'=",C);	
	int ret = printDiffent("A*A'=I", X, C);
	freeMatrix(A);
	freeMatrix(B);
	freeMatrix(C);	
	freeMatrix(X);
	return ret;
}

static int test_inverse_pivoting()
{
	real * A = makeRandMatrix();
	real * B = makeMatrix();
	real * C = makeMatrix();
	real * X = makeMatrix();

	identity(X, N);
	identity(A, N);
	for(int i=0;i<N;i++){
		int m = (int)((N-1)*((real)rand()/(real)RAND_MAX));
		if(i!=m && m>=0 && m<N){
			xchangeRaw(A,N,i,m);
		}
	}
	copyMatrix(B, A);
	printMat("A=",A);
	inverse_pivoting(A, N);
	printMat("inverse=",A);
	multiply0(C, B, A, N, N, N);
	printMat("A*A'=",C);
	int ret = printDiffent("A*A'=I", X, C);
	freeMatrix(A);
	freeMatrix(B);
	freeMatrix(C);
	freeMatrix(X);
	return ret;
}

static int test_inverse_diagonal()
{
	real * A = makeRandMatrix();
	real * B = makeMatrix();
	real * C = makeMatrix();	
	real * X = makeMatrix();

	identity(X, N);
	clearLowerTriangle(A,N);
	clearUpperTriangle(A,N);
	copyMatrix(B,A);
	printMat("A=",A);
	inverse_diagonal(A, N);
	printMat("inverse=",A);
	multiply0(C, B, A, N, N, N);
	printMat("A*A'=",C);
	int ret = printDiffent("A*A'=I", X, C);
	freeMatrix(A);
	freeMatrix(B);
	freeMatrix(C);
	freeMatrix(X);
	return ret;
}
/*
static int test_inverse_1()
{
	real * A = makeRandMatrix();
	real * L = makeMatrix(); 
	real * C = makeMatrix();
	real * P = makeMatrix();
	real * D = makeMatrix();
	real * T = makeMatrix();
	real * X = makeMatrix();

	identity(X, N);
	copyMatrix(T,A);
	printMat("A=",A);
	pldu(A, P, D, L, N);
	printMat("P=",P);
	printMat("D=",D);
	printMat("U=",A);
	printMat("L=",L);
	int r = inverse0(P, L, D, A, N);
	printMat("inverse=",P);
	multiply0(C, T, P, N, N, N);
	printMat("A*A'=",C);
	int ret = printDiffent("A*A'=I", X, C);
	freeMatrix(A);
	freeMatrix(L);
	freeMatrix(C);
	freeMatrix(P);
	freeMatrix(D);
	freeMatrix(T);
	freeMatrix(X);
	return r?ret:-1;
}
*/

static int test_crout_lu()
{
	real * A = makeRandMatrix();
	real * L = makeMatrix();
	real * C = makeMatrix();
	real * X = makeMatrix();

	copyMatrix(X, A);
	printMat("A=",A);
	int r = crout_lu(A, L, N);
	if(!r && N<5 ){
		printf("crout_lu failed\n");
		disablePrint(0);
		printMat3(X,L,A,N);
		disablePrint(1);
	}	
	printMat("U=",A);
	printMat("L=",L);
	multiply0(C, L, A, N, N, N);
	printMat("L*U=",C);	
	int ret = printDiffent("A=L*U", X, C);
	freeMatrix(A);
	freeMatrix(L);
	freeMatrix(C);
	freeMatrix(X);
	return r?ret:-1;
}

static int test_crout_plu()
{
	real * A = makeRandMatrix();
	real * L = makeMatrix();
	real * C = makeMatrix();
	real * P = makeMatrix();
	real * X = makeMatrix();

	copyMatrix(X, A);
	printMat("A=",A);
	int r = crout_plu(A, P,L, N);
	if(!r && N<5 ){
		printf("crout_plu failed\n");
		disablePrint(0);
		printMat3(X,L,A,N);
		disablePrint(1);
	}	
	printMat("U=",A);
	printMat("L=",L);
	printMat("P=",P);
//	multiply0(C, L, A, N, N, N);
//	multiply0(A, C, P, N, N, N);
	
	multiply0(C, P, L, N, N, N);
	multiply0(P, C, A, N, N, N);	
	copyMatrix(A,P);
	
	printMat("L*U=",A);	
	int ret = printDiffent("A=L*U", X, A);
	freeMatrix(A);
	freeMatrix(L);
	freeMatrix(C);
	freeMatrix(P);
	freeMatrix(X);
	return r?ret:-1;	
}

void test( const char *s,int (*tf)() )
{
	int c = 0;
	int d = 0;
	double t0 = getClock();
	for(int i=0;i<100;i++){
		int r = tf();
		if( r==1 ){
			c++;
		}
		else if( r==-1 ){
			d++;
		}
	}
	printf("test %s [%d%%|%d%%|%d%%] use time : %.2fs\n",s,c,d,c+d,getClock()-t0);
}

/*
 * 关于运算次数和进度之间的关系
 */
 double random_double()
 {
	 return (double)rand()/(double)RAND_MAX;
 }
 
 void test_p()
 {
	 double t,oo;
	 double o = random_double();
	 oo = o;
	 float s = (float)o;
	 for( int i=0;i<100;i++){
		 t = random_double();
		 o *= t;
		 s *= (float)t;
	 }
	 double e = o - (double)s;
	 printf("%f , %f , %f\n",oo,e,e/oo);
	 float a = 0.65f;
	 float b = 0.6f;
	 float c = a-b;
	 printf("%f\n",c);
 }
 
int main(int argn,const char *argv[])
{
	disablePrint(1);
	readom_init();
	test("lu",test_lu_1);
	//test("pldu",test_pldu_1);
	test("inverse_low_triangle",test_inverse_low_triangle);
	test("inverse_upper_triangle",test_inverse_upper_triangle);
	test("inverse_pivoting",test_inverse_pivoting);
	test("inverse_diagonal",test_inverse_diagonal);
	//test("inverse",test_inverse_1);
	//因为不交互行，因此有些能分解也不能进行正常分解
	//test("crout_lu",test_crout_lu);
	test("test_crout_plu",test_crout_plu);
	//test_p();
	return 0;
}

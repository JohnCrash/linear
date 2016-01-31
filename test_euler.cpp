#include "linear.h"
#include <stdio.h>
#include <math.h>

real func_0(real t,real x)
{
	return x+1;
}

real func_0_solve(real t)
{
	return pow(M_E,t)-1;
}

/*
 * Euler method solve difference equation
 */
real eulerStep(real (*f)(real,real),real t0,real x0,real step)
{
	//x0+f(t0)* h
	return x0+f(t0,x0)*step;
}

real eulerPrint(real (*f)(real,real),real t0,real x0,real step,int n)
{
	real t,x,h;
	t = t0;
	x = x0;
	h = step;	
	printf("euler method,initial value t=%.4f x=%.4f step=%.4f\n",t0,x0,step);
	printf("%.4f\t|%.4f\n",t,x);
	for(int i=0;i<n;i++){
		x = eulerStep(func_0,t,x,h);
		t+=h;
		printf("%.4f\t|%.4f\n",t,x);
	}	
}

real solvePrint(real (*f)(real),real t0,real x0,real step,int n)
{
	real t,x,h;
	t = t0;
	x = x0;
	h = step;
	printf("formula solve,initial value t=%.4f x=%.4f step=%.4f\n",t0,x0,step);
	printf("%.4f\t|%.4f\n",t,x);	
	for(int i=0;i<n;i++){
		x = f(t);
		t+=h;
		printf("%.4f\t|%.4f\n",t,x);
	}	
}

//x'=1+x , x(0)=1
void test_1()
{
	eulerPrint(func_0,0,1,0.1,10);
	solvePrint(func_0_solve,0,1,0.1,10);
}

int main(int argn,char *argv[])
{
	test_1();
	return 0;
}

#include "linear.h"
#include <stdio.h>
#include <math.h>

/*
 * 形如：x'=1+1的ODE
 */
real func_0(real t,real x)
{
	return x+1;
}

/*
 * x'=1+1的公式解，A是系数，x'=1+x时A=2
 * 但是其实可以根据初始条件计算出A
 * b=false求A，这时候t和A给出初值，b=true做正常函数计算
 */
real func_0_solve(real t,real A,bool b)
{
	if(b)
		return A*pow(M_E,t)-1;
	else
		return (A+1)/pow(M_E,t);
}

/*
 * Euler method solve difference equation
 * 使用起点的斜率做计算
 */
real eulerStep(real (*f)(real,real),real t0,real x0,real step)
{
	//x0+f(t0)* h
	return x0+f(t0,x0)*step;
}

/*
 * 终点点法使用step两端的斜率的平均值作为斜率进行计算
 * 起点斜率是f(t0,x0)，终点斜率是f(t0+step,x0+f(t0,x0)*step)
 * 终点是通过euler法进行的近似(t0+step,x0+f(t0,x0)*step)
 */
real midpointStep(real (*f)(real,real),real t0,real x0,real step)
{
	return x0+(f(t0,x0)+f(t0+step,x0+f(t0,x0)*step))*step/2;
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

real solvePrint(real (*f)(real,real,bool),real t0,real x0,real step,int n)
{
	real t,x,h,A;
	t = t0;
	x = x0;
	h = step;
	A = f(t,x,false);
	printf("formula solve,initial value t=%.4f x=%.4f step=%.4f\n",t0,x0,step);
	printf("%.4f\t|%.4f\n",t,x);	
	for(int i=0;i<n;i++){
		t+=h;
		x = f(t,A,true);
		printf("%.4f\t|%.4f\n",t,x);
	}	
}

real midpointPrint(real (*f)(real,real),real t0,real x0,real step,int n)
{
	real t,x,h;
	t = t0;
	x = x0;
	h = step;
	printf("midpoint method,initial value t=%.4f x=%.4f step=%.4f\n",t0,x0,step);
	printf("%.4f\t|%.4f\n",t,x);	
	for(int i=0;i<n;i++){
		x=midpointStep(func_0,t,x,h);
		t+=h;
		printf("%.4f\t|%.4f\n",t,x);
	}		
}

//x'=1+x , x(0)=1
void test_1()
{
	printf("solve x'=1+x\n");
	eulerPrint(func_0,0,1,0.1,10);
	printf("x=2*exp(t)-1\n");
	solvePrint(func_0_solve,0,1,0.1,10);
	printf("midpoint\n");
	midpointPrint(func_0,0,1,0.1,10);	
}

void test_2()
{
}

int main(int argn,char *argv[])
{
	test_1();
	return 0;
}

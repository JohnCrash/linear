#include "linear.h"
#include <stdio.h>
#include <math.h>

enum FormulaParamter
{
	FORMULA,
	DERIVATIVE,
	SOLVE,
	INITIAL,
};

typedef real (*OdeFormula)(real,real,FormulaParamter);

/*
 * x'=1+x的公式解，A是系数，x'=1+x时A=2
 * 但是其实可以根据初始条件计算出A
 * DERIVATIVE给出微分式，SOLVE给出公式解
 * INITIAL给出初值求的系数
 */
real func_0(real t,real A,FormulaParamter s)
{
	switch(s){
		case DERIVATIVE:
			return A+1;
		case SOLVE:
			return A*pow(M_E,t)-1;
		case INITIAL:
			return (A+1)/pow(M_E,t);
	}
	//assert(0,"error func_0_solve FormulaParamter");
	return 0;
}

/*
 * y'=x+y
 */
real func_1(real t,real A,FormulaParamter s)
{
	switch(s){
		case DERIVATIVE:
			return A+t;
		case SOLVE:
			return A*pow(M_E,t)-1-t;
		case INITIAL:
			return (A+1+t)/pow(M_E,t);
	}
	//assert(0,"error func_0_solve FormulaParamter");
	return 0;
}

/*
 * y'=y*sin(x)
 */
real func_2(real t,real A,FormulaParamter s)
{
	switch(s){
		case DERIVATIVE:
			return A*sin(t);
		case SOLVE:
			return A*pow(M_E,-cos(t));
		case INITIAL:
			return A/pow(M_E,-cos(t));
	}
	//assert(0,"error func_0_solve FormulaParamter");
	return 0;
}

/*
 * y'=-x/y
 */
real func_3(real t,real A,FormulaParamter s)
{
	switch(s){
		case DERIVATIVE:
			return -t/A;
		case SOLVE:
			return sqrt(2*A-t*t);
		case INITIAL:
			return (A*A+t*t)/2;
	}
	//assert(0,"error func_0_solve FormulaParamter");
	return 0;
}

/*
 * Euler method solve difference equation
 * 使用起点的斜率做计算
 */
real eulerStep(OdeFormula f,real t0,real x0,real step)
{
	//x0+f(t0)* h
	return x0+f(t0,x0,DERIVATIVE)*step;
}

/*
 * 终点点法使用step两端的斜率的平均值作为斜率进行计算
 * 起点斜率是f(t0,x0)，终点斜率是f(t0+step,x0+f(t0,x0)*step)
 * 终点是通过euler法进行的近似(t0+step,x0+f(t0,x0)*step)
 */
real midpointStep(OdeFormula f,real t0,real x0,real step)
{
	//dx = f(t0,x0)*step
	//fmid=(f(t0,x0)+f(t0+step,x0+dx))/2
	//x0+step*fmid
	real fmid,dx;
	dx=f(t0,x0,DERIVATIVE)*step;
	fmid=(f(t0,x0,DERIVATIVE)+f(t0+step,x0+dx,DERIVATIVE))/2;
	return x0+step*fmid;
}

/*
 * 测试不同的方法,David Baraff的Differential Equation Basics A.pdf
 * 中提到的中点法，和中点法基本一样，精度在O3
 */
real midpointStep2(OdeFormula f,real t0,real x0,real step)
{
	//dx = step*f(t0,x0)
	//fmid=f((x0+dx)/2,(t0+step)/2) ，我认为这个公式是错误的
	//应该写成fmid=f(x0+dx/2,t0+step/2)
	//x0+h*fmid
	real fmid,dx;
	dx = step*f(t0,x0,DERIVATIVE);
	fmid = f(t0+step/2,x0+dx/2,DERIVATIVE); 
	return x0+step*fmid;
}

/*
 * Runge-Kutta
 * 精度在O5
 */
real testStep(OdeFormula f,real t0,real x0,real step)
{
	real k1,k2,k3,k4;
	k1 = step*f(t0,x0,DERIVATIVE); //欧拉法
	k2 = step*f(t0+step/2,x0+k1/2,DERIVATIVE); //中点法
	k3 = step*f(t0+step/2,x0+k2/2,DERIVATIVE);
	k4 = step*f(t0+step,x0+k3,DERIVATIVE);
	return x0+k1/6+k2/3+k3/3+k4/6;
}
/*
 * 打印通过欧拉法
 */
void eulerPrint(OdeFormula f,real t0,real x0,real step,int n)
{
	real t,x,h;
	t = t0;
	x = x0;
	h = step;	
	printf("euler method,initial value t=%.4f x=%.4f step=%.4f\n",t0,x0,step);
	printf("%.4f\t|%.4f\n",t,x);
	for(int i=0;i<n;i++){
		x = eulerStep(f,t,x,h);
		t+=h;
		printf("%.4f\t|%.4f\n",t,x);
	}	
}

void solvePrint(OdeFormula f,real t0,real x0,real step,int n)
{
	real t,x,h,A;
	t = t0;
	x = x0;
	h = step;
	A = f(t,x,INITIAL);
	printf("formula solve,initial value t=%.4f x=%.4f step=%.4f\n",t0,x0,step);
	printf("%.4f\t|%.4f\n",t,x);	
	for(int i=0;i<n;i++){
		t+=h;
		x = f(t,A,SOLVE);
		printf("%.4f\t|%.4f\n",t,x);
	}	
}

/*
 * 中点法
 */
void midpointPrint(OdeFormula f,real t0,real x0,real step,int n)
{
	real t,x,h;
	t = t0;
	x = x0;
	h = step;
	printf("midpoint method,initial value t=%.4f x=%.4f step=%.4f\n",t0,x0,step);
	printf("%.4f\t|%.4f\n",t,x);	
	for(int i=0;i<n;i++){
		x=midpointStep(f,t,x,h);
		t+=h;
		printf("%.4f\t|%.4f\n",t,x);
	}		
}

/*
 * 对三种方法比较进行比较
 */
void comparePrint(OdeFormula f,real t0,real x0,real step,int n)
{
	real t,x,h,mx,sx,A,tx;
	t = t0;
	tx = sx = mx=x = x0;
	h = step;	
	A = f(t,x,INITIAL);
	printf("变量\t|欧拉法\t误差\t|中点法\t误差\t误差比\t|测试\t误差|公式法|\n");
	for(int i=0;i<n;i++){
		x = eulerStep(f,t,x,h);
		mx = midpointStep(f,t,mx,h);
		tx = testStep(f,t,tx,h);
		t+=h;
		sx = f(t,A,SOLVE);
		printf("%.4f\t|%.4f\t%.4f\t|%.4f\t%.4f\t%.4f\t|%.4f\t%.7f|%.4f|\n",t,x,sx-x,mx,mx-sx,fabs(mx-sx)/fabs(sx-x),tx,tx-sx,sx);
	}		
}
//x'=1+x , x(0)=1
void test_1()
{
	printf("solve x'=1+x\n");
	eulerPrint(func_0,0,1,0.1,10);
	printf("x=2*exp(t)-1\n");
	solvePrint(func_0,0,1,0.1,10);
	printf("midpoint\n");
	midpointPrint(func_0,0,1,0.1,10);	
}

void test_2()
{
	printf("x'=1+x,x(0)=1\n");
	comparePrint(func_0,0,1,0.1,10);	
	printf("y'=y+x,y(0)=1\n");
	comparePrint(func_1,0,1,0.1,10);		
	printf("y'=y sin(x),y(0)=1\n");
	comparePrint(func_2,-1,1.0/M_E,0.3,10);		
	printf("y'=-x/y,y(0)=1\n");
	comparePrint(func_3,-1,2.0,0.3,10);	
}

int main(int argn,char *argv[])
{
	//test_1();
	test_2();
	return 0;
}

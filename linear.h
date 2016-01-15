#ifndef _LINEAR_H_
#define _LINEAR_H_

typedef float real;
#define fabs(x) ((x)>0?(x):-(x))

void zero(real * A,int n);
void identity(real * A,int n);
int lu(real * A,real * L,int n);

#endif
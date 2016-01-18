#include "linear.h"

void multiply0(real *A, const real *B, const real *C, int p, int q, int r)
{
    const int qskip = q;
    const int rskip = r;
    real *aa = A;
    const real *bb = B;
    for (int i=p; i; aa+=rskip, bb+=qskip, --i) {
        real *a = aa;
        const real *cc = C, *ccend = C + r;
        for (; cc != ccend; ++a, ++cc) {
            real sum = REAL(0.0);
            const real *c = cc;
            const real *b = bb, *bend = bb + q;
            for (; b != bend; c+=rskip, ++b) {
                sum += (*b)*(*c);
            }
            (*a) = sum; 
        }
    }
}

void multiply1(real *A, const real *B, const real *C, int p, int q, int r)
{
    const int pskip = p;
    const int rskip = r;
    real *aa = A;
    const real *bb = B, *bbend = B + p;
    for (; bb != bbend; aa += rskip, ++bb) {
        real *a = aa;
        const real *cc = C, *ccend = C + r;
        for (; cc != ccend; ++a, ++cc) {
            real sum = REAL(0.0);
            const real *b = bb, *c = cc;
            for (int k=q; k; b+=pskip, c+=rskip, --k) {
                sum += (*b)*(*c);
            }
            (*a) = sum;
        }
    }
}


void multiply2(real *A, const real *B, const real *C, int p, int q, int r)
{
    const int rskip = r;
    const int qskip = q;
    real *aa = A;
    const real *bb = B;
    for (int i=p; i; aa+=rskip, bb+=qskip, --i) {
        real *a = aa, *aend = aa + r;
        const real *cc = C;
        for (; a != aend; cc+=qskip, ++a) {
            real sum = REAL(0.0);
            const real *b = bb, *c = cc, *cend = cc + q;
            for (; c != cend; ++b, ++c) {
                sum += (*b)*(*c);
            }
            (*a) = sum; 
        }
    }
}
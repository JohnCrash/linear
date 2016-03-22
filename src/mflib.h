#ifndef _MFLIB_H_
#define _MFLIB_H_

#include "config.h"

#ifndef HAVE_FABS
#define fabs(x) ((x)>0?(x):-(x))
#endif

#ifndef HAVE_FMAX
#define fmax(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef HAVE_FMIN
#define fmin(a,b) ((a)<(b)?(a):(b))
#endif

#endif
#ifndef _MFLIB_H_
#define _MFLIB_H_

#include "config.h"

#define ODE_PURE_INLINE inline
#ifdef WIN32
	#define HAVE_FABS 1
	#define HAVE_FMAX 1
	#define HAVE_FMIN 1
	#define HAVE_FLT_MIN 1
	#define HAVE_FLT_MAX 1

	#ifdef _DLLEXPORT_
	#define DYNFUNC __declspec( dllexport )
	#else
	#define DYNFUNC __declspec( dllimport )
	#endif
#else
	#define DYNFUNC
#endif

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
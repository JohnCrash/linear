lib_LTLIBRARIES=liblinear.la
liblinear_la_SOURCES=lu.cpp misc.cpp matrix.cpp
liblinear_ladir=linear
liblinear_la_HEADERS=linear.h misc.h
liblinear_la_LDFLAGS=-no-undefined -avoid-version

test_lu_LDADD=liblinear.la
test_inverse_LDADD=liblinear.la

bin_PROGRAMS=test_lu test_inverse
test_lu_SOURCES=test_lu.cpp
test_inverse_SOURCES=test_inverse.cpp
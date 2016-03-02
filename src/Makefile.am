AM_CXXFLAGS=-std=c++11
lib_LTLIBRARIES=liblinear.la

liblinear_la_SOURCES=lu.cpp misc.cpp matrix.cpp \
solve.cpp lcp_pgs.cpp \
lcp.cpp lcp_lemke.cpp lcp_pivot.cpp \
mlcp.cpp mlcp_lemke.cpp mlcp_pivot.cpp mlcp_pgs.cpp

liblinear_ladir=linear
liblinear_la_HEADERS=linear.h misc.h lcp.h
liblinear_la_LDFLAGS=-no-undefined -avoid-version

test_mlcp_LDADD=liblinear.la
test_lcp_LDADD=liblinear.la
test_lu_LDADD=liblinear.la
test_inverse_LDADD=liblinear.la
test_LDADD=liblinear.la
test_solve_LDADD=liblinear.la

bin_PROGRAMS=test_lu test_inverse test_solve test test_euler test_prime test_lcp test_mlcp

test_prime_SOURCES=test_prime.cpp
test_euler_SOURCES=test_euler.cpp
test_solve_SOURCES=test_solve.cpp
test_lu_SOURCES=test_lu.cpp
test_inverse_SOURCES=test_inverse.cpp
test_SOURCES=test.cpp
test_lcp_SOURCES=test_lcp.cpp
test_mlcp_SOURCES=test_mlcp.cpp
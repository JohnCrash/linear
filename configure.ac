AC_PREREQ([2.69])
AC_INIT([linear], [0.1], [xxxxx@sina.com])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

AC_CONFIG_SRCDIR([src/linear.h])
AC_CONFIG_HEADERS([src/config.h])
AM_PROG_AR

LT_PREREQ([2.4.0])
LT_INIT([disable-static win32-dll dlopen])
AC_PROG_LIBTOOL

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([memory.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

# Process this file with autoconf to produce a configure script.
AC_INIT(harminv, 0.2)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR(harminv.c)
AM_CONFIG_HEADER(config.h)

# Checks for programs.
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_INSTALL

# Fortran libraries are (probably) required for BLAS and LAPACK:
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS

# Add lots of compiler warnings to check for if we are using gcc:
# (The variable $GCC is set to "yes" by AC_PROG_CC if we are using gcc.)
if test "$GCC" = "yes"; then
	CFLAGS="$CFLAGS -Wall -W -Wbad-function-cast -Wcast-qual -Wpointer-arith -Wcast-align -pedantic"
fi

HARMINV_CC="$CC"
HARMINV_CFLAGS="$CFLAGS"
AC_SUBST(HARMINV_CC)
AC_SUBST(HARMINV_CFLAGS)

###########################################################################

AC_ARG_WITH(cxx, [  --with-cxx=<dir>        force use of C++ and complex<double>], with_cxx=$withval, with_cxx=no)

if test "x$with_cxx" = "xno"; then

have_c_complex=yes
AC_CACHE_CHECK([for C complex keyword], acx_cv_c_complex,
 [acx_cv_c_complex=unsupported
  for acx_kw in complex _Complex __complex__; do
    AC_TRY_COMPILE([#include <complex.h>], [float $acx_kw foo;],
                   [acx_cv_c_complex=$acx_kw; break])
  done
 ])

if test "$acx_cv_c_complex" = "unsupported"; then
   AC_MSG_WARN([C doesn't support complex numbers; switching to C++.])
   have_c_complex=no
else
   if test "$acx_cv_c_complex" != "complex"; then
     AC_DEFINE_UNQUOTED(complex, $acx_cv_c_complex, [Define to equivalent of C99 complex keyword.  Do not define if complex is supported directly.])
   fi
fi

AC_CACHE_CHECK([for C complex I constant], acx_cv_c_complex_I,
 [acx_cv_c_complex_I=unsupported
 for acx_kw in I "(1.0i)"; do
   AC_TRY_COMPILE([#include <complex.h>],
                  [float $acx_cv_c_complex foo = $acx_kw;],
                  [acx_cv_c_complex_I=$acx_kw; break])
 done
])

if test "$acx_cv_c_complex_I" = "unsupported"; then
   AC_MSG_ERROR([C doesn't support complex numbers.])
   have_c_complex=no
else
   if test "$acx_cv_c_complex_I" != "I"; then
     AC_DEFINE_UNQUOTED(I, $acx_cv_c_complex_I, [Define to equivalent of C99 complex I constant.  Do not define if I is supported directly.])
   fi
fi

AC_CACHE_CHECK([for CX_LIMITED_RANGE_ON macro], acx_cv_c_cx_limited_range_on,
 [acx_cv_c_cx_limited_range_on=no
  AC_TRY_COMPILE([#include <complex.h>],
                 [CX_LIMITED_RANGE_ON;], [acx_cv_c_cx_limited_range_on=yes])])
if test "$acx_cv_c_cx_limited_range_on" = "no"; then
   AC_DEFINE(CX_LIMITED_RANGE_ON, 0, [Define to no-op if C99 macro is not available.])
fi

fi # $with_cxx = no

if test "$with_cxx" = "yes" -o "$have_c_complex" = "no"; then
	AC_PROG_CXX
	HARMINV_CC="$CXX"
	HARMINV_CFLAGS="$CXXFLAGS"
fi

if test "$have_c_complex" = "yes"; then
	AC_CHECK_FUNCS(carg)
fi

###########################################################################

# Checks for header files.
AC_CHECK_HEADERS(unistd.h getopt.h)

AC_CHECK_LIB(m, sqrt)

###########################################################################

# Checks for BLAS/LAPACK libraries:

sinclude(acx_blas.m4)
ACX_BLAS([], [AC_MSG_ERROR([BLAS library not found])])
sinclude(acx_lapack.m4)
ACX_LAPACK([], [AC_MSG_ERROR([LAPACK library not found])])

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

###########################################################################

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

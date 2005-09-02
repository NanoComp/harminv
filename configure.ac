# Process this file with autoconf to produce a configure script.
AC_INIT(harminv, 1.3, stevenj@alum.mit.edu)
AM_INIT_AUTOMAKE(1.6)
AC_CONFIG_SRCDIR(harminv.c)
AM_CONFIG_HEADER(config.h)
AM_MAINTAINER_MODE

# Shared-library version number; indicates api compatibility, and is
# *not* the same as the "public" version number.  CURRENT:REVISION:AGE
SHARED_VERSION_INFO="2:3:0"
AC_SUBST(SHARED_VERSION_INFO)
AM_ENABLE_SHARED(no) dnl shared libs cause too many headaches to be default

##############################################################################

AC_ARG_ENABLE(debug, [AC_HELP_STRING([--enable-debug,compile for debugging])],
	      enable_debug=$enableval, enable_debug=no)
if test "$enable_debug" = "yes"; then
	CFLAGS="-g"
	CXXFLAGS="-g"
	FFLAGS="-g"
	AC_DEFINE(DEBUG,1,[define to enable debugging code])
fi

##############################################################################

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
if test "$GCC" = "yes" && test "$USE_MAINTAINER_MODE" = yes; then
	CFLAGS="$CFLAGS -Wall -W -Wbad-function-cast -Wcast-qual -Wpointer-arith -Wcast-align -pedantic"
fi

# Checks for header files.
AC_CHECK_HEADERS(unistd.h getopt.h)

###########################################################################

AC_ARG_WITH(cxx, [  --with-cxx=<dir>        force use of C++ and complex<double>], with_cxx=$withval, with_cxx=no)

if test "x$with_cxx" = "xno"; then

have_c_complex=yes
AC_CHECK_HEADERS(complex.h)
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

fi # C has complex keyword

fi # $with_cxx = no

if test "$have_c_complex" = "yes"; then
	AC_CHECK_FUNCS(carg)
fi

###########################################################################

AC_PROG_CXX
if test "$with_cxx" = "yes" -o "$have_c_complex" = "no"; then
	CC="$CXX"
	CFLAGS="$CXXFLAGS"
fi

###########################################################################

AC_CHECK_TYPES([long double])

###########################################################################

AC_PROG_LIBTOOL                                                                                
###########################################################################

AC_CHECK_LIB(m, sqrt)

###########################################################################

# Checks for BLAS/LAPACK libraries:

ACX_BLAS([], [AC_MSG_ERROR([BLAS library not found])])
ACX_LAPACK([], [AC_MSG_ERROR([LAPACK library not found])])

LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

###########################################################################

AC_CONFIG_FILES([Makefile harminv.pc])
AC_OUTPUT

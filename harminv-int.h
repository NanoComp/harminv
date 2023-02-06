/* Copyright (C) 2017 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef HARMINV_INT_H
#define HARMINV_INT_H 1

#include "config.h"

#ifndef __cplusplus
/* Require C99 complex number support; this is just too painful
   without it.  Alternatively, use the complex<double> STL class in
   C++. */

#include <complex.h>
#endif

#include "harminv.h"

typedef harminv_complex cmplx; /* shortcut */

/* define extended precision (long double) cmplxl type */
#if defined(HAVE_LONG_DOUBLE) && !defined(__cplusplus)
#if defined(__cplusplus)
#include <complex>
typedef std::complex<long double> cmplxl; /* problematic w/some C++ comp's */

#elif defined(_Complex_I) && defined(complex) && defined(I)
/* C99 <complex.h> header was included before harminv.h */
typedef long double _Complex cmplxl;

#endif
#else
typedef cmplx cmplxl;
#endif

#ifdef __cplusplus
#define I cmplx(0, 1)
#define creal(c) real(c)
#define cimag(c) imag(c)
#define cabs(c) abs(c)
#define carg(c) arg(c)
#define cexp(c) exp(c)
#define csqrt(c) sqrt(c)
#define clog(c) log(c)
#else
#ifndef HAVE_CARG /* Cray doesn't have this for some reason */
#define carg(c) atan2(cimag(c), creal(c))
#endif
#endif

struct harminv_data_struct {
  const cmplx *c;
  int n, K, J, nfreqs;
  double fmin, fmax;
  cmplx *z;
  cmplx *U0, *U1;
  cmplxl *G0, *G0_M, *D0; /* cached G, G_M, and D arrays for U0 */
  cmplx *B, *u;           /* eigen-solutions of U1*B = u*U0*B */
  cmplx *amps;            /* mode amplitudes */
  double *errs;           /* relative "error" estimates */
};

#endif /* HARMINV_INT_H */

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

#ifndef HARMINV_H
#define HARMINV_H

/* the following need to be kept in sync with configure.ac: */
#define HARMINV_VERSION_MAJOR 1
#define HARMINV_VERSION_MINOR 4
#define HARMINV_VERSION_PATCH 1

/**************************************************************************/

#if defined(__cplusplus)
#include <complex>
extern "C" {
typedef std::complex<double> harminv_complex;

#elif defined(_Complex_I) && defined(complex) && defined(I)
/* C99 <complex.h> header was included before harminv.h */
typedef double _Complex harminv_complex;

#else
typedef double harminv_complex[2];
#endif

typedef struct harminv_data_struct *harminv_data;

/**************************************************************************/

typedef int (*harminv_mode_ok_func)(harminv_data d, int k, void *);

extern harminv_data harminv_data_create(int n, const harminv_complex *signal, double fmin,
                                        double fmax, int nf);
extern void harminv_data_destroy(harminv_data d);

extern void harminv_solve(harminv_data d);

extern int harminv_get_num_freqs(harminv_data d);
extern double harminv_get_freq(harminv_data d, int k);
extern double harminv_get_Q(harminv_data d, int k);
extern double harminv_get_decay(harminv_data d, int k);
extern void harminv_get_omega(harminv_complex *omega, harminv_data d, int k);
extern void harminv_get_amplitude(harminv_complex *amplitude, harminv_data d, int k);
extern double harminv_get_freq_error(harminv_data d, int k);

/* the following routines are undocumented and not recommended */
extern void harminv_solve_once(harminv_data d);
extern void harminv_solve_again(harminv_data d, harminv_mode_ok_func ok, void *ok_d);
extern void harminv_solve_ok_modes(harminv_data d, harminv_mode_ok_func ok, void *ok_d);
extern double *harminv_compute_freq_errors(harminv_data d);
extern harminv_complex *harminv_compute_amplitudes(harminv_data d);

/**************************************************************************/

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* HARMINV_H */

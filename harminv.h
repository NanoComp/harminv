#ifndef HARMINV_H
#define HARMINV_H

/* Require C99 complex number support; this is just too painful
   without it. */
#include <complex.h>

#include "config.h"

/**************************************************************************/

typedef double complex cmplx;

typedef struct harminv_data_struct {
     const cmplx *c;
     int n, K, J, nfreqs;
     double fmin, fmax;
     cmplx *z;
     cmplx *U0, *U1;
     cmplx *B, *u;  /* eigen-solutions of U1*B = u*U0*B */
} *harminv_data;

/**************************************************************************/

extern harminv_data harminv_data_create(int n,
					const cmplx *signal,
					double fmin, double fmax, int nf);
extern void harminv_data_destroy(harminv_data d);

extern void harminv_solve(harminv_data d);
extern void harminv_solve_again(harminv_data d);

extern int harminv_get_num_freqs(harminv_data d);
extern double harminv_get_freq(harminv_data d, int k);
extern double harminv_get_decay(harminv_data d, int k);

extern double *harminv_compute_frequency_errors(harminv_data d);
extern cmplx *harminv_compute_amplitudes(harminv_data d);

/**************************************************************************/

#endif /* HARMINV_H */

/* Copyright (C) 2004 Massachusetts Institute of Technology.
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

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "harminv-int.h"
#include "check.h"
#include "copyright.h"

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif

/* eat whitespace, including #... comments, from the file.  Returns the
   number of newlines read (so that a line count can be maintained).  If
   echo_comments != 0, then echo #... comments to stdout.  Commas count
   as whitespace, so that we can read comma-delimited text. */
static int eat_whitespace(FILE *f, int echo_comments)
{
     int c, newlines = 0;
     do {
	  do {
	       c = getc(f);
	       newlines += c == '\n';
	  } while (isspace(c) || c == ',');
	  
	  if (c == '#') { /* # begins comments that extend to the newline */
	       if (echo_comments)
		    putc(c, stdout);
	       do {
		    c = getc(f);
		    if (echo_comments) {
			 if (c != EOF)
			      putc(c, stdout);
			 else /* terminate line if we hit EOF */
			      putc('\n', stdout);
		    }
		    newlines += c == '\n';
	       } while (c != EOF && c != '\n');
	  }
     } while (isspace (c));
     ungetc(c, f); /* put back the last character read */
     newlines -= c == '\n';
     return newlines;
}

static int eat_plus(FILE *f)
{
     int c = getc(f);
     if (c != EOF && c != '+')
	  ungetc(c, f);
     return (c == '+' || c == '-');
}

static int eat_i(FILE *f)
{
     int c = getc(f);
     if (c != EOF && tolower(c) != 'i')
	  ungetc(c, f);
     return (tolower(c) == 'i');
}

static cmplx *read_input_data(FILE *f, int *n, int verbose)
{
    cmplx *data = NULL;
    int line = 1, n_alloc = 0;
    *n = 0;

     do {
	  double re, im;
	  int nread;

	  line += eat_whitespace(f, verbose);
	  nread = fscanf(f, "%lg", &re);
	  if (nread == 1 && eat_plus(f)) {
	       nread = fscanf(f, "%lg", &im);
	       if (nread == 1) nread = eat_i(f);
	  }
	  else
	       im = 0.0;
	  if (nread != EOF) {
	       if (nread < 1) {
		    fprintf(stderr, "harminv: invalid input on line %d.\n",
			    line);
		    free(data); 
		    *n = 0;
		    return NULL;
	       }
	       if (*n >= n_alloc) {
		    n_alloc = (n_alloc + 1) * 2;
		    data = (cmplx*) realloc(data, sizeof(cmplx) * n_alloc);
		    CHECK(data, "out of memory");
	       }
	       data[*n] = re + I*im;
	       ++*n;
	  }
     } while (!feof(stdin));

     data = (cmplx*) realloc(data, sizeof(cmplx) * *n);
     return data;
}

#ifdef INFINITY
const double inf = INFINITY;
#else
const double inf = 1.0 / 0.0;
#endif


#define NF 100
#define ERR_THRESH 0.1
#define REL_ERR_THRESH inf
#define AMP_THRESH 0.0
#define REL_AMP_THRESH -1.0
#define Q_THRESH 10.0

static void usage(FILE *f)
{
     fprintf(f, "Usage: harminv [options] <freq-min>-<freq-max>...\n"
	     "Options: \n"
	     "         -h : this help message\n"
             "         -V : print version number and copyright\n"
             "         -v : verbose output\n"
	     "         -T : specify periods instead of frequencies\n"
	     "         -w : specify/output angular frequency, not frequency\n"
	     "    -t <dt> : specify sampling interval dt [default: 1]\n"
	     "    -f <nf> : specify initial spectral density [default: %d]\n"
	     "  -s <sort> : sort by <sort> = freq/err/decay/amp [default: freq]\n"
	     "        -F : discard frequencies outside of specified range\n"
	     "    -a <a> : discard amplitudes < max * <a> [default: %e]\n"
	     "    -A <A> : discard amplitudes < <A> [default: %g]\n"
	     "    -e <e> : discard relative errors > min * <e> [default: %e]\n"
	     "    -E <E> : discard relative errors > <E> [default: %e]\n"
	     "    -Q <Q> : discard Q > <E> [default: %g]\n",
	     NF,
	     AMP_THRESH, REL_AMP_THRESH,
	     ERR_THRESH, REL_ERR_THRESH,
	     Q_THRESH);
}

#define TWOPI 6.2831853071795864769252867665590057683943388

harminv_data hd;

enum {
     SORT_FREQUENCY, SORT_DECAY, SORT_ERROR, SORT_AMPLITUDE, SORT_Q
} sortby = SORT_FREQUENCY;

static int cmp(double a, double b)
{
     return a > b ? 1 : (a < b ? -1 : 0);
}

static int compar(const void *a, const void *b)
{
     const int *ia = (const int *) a;
     const int *ib = (const int *) b;

     switch (sortby) {
	 case SORT_FREQUENCY:
	      return cmp(harminv_get_freq(hd,*ia), harminv_get_freq(hd,*ib));
	 case SORT_DECAY:
	      return cmp(harminv_get_decay(hd,*ia), harminv_get_decay(hd,*ib));
	 case SORT_ERROR:
	      return cmp(harminv_get_freq_error(hd, *ia), 
			 harminv_get_freq_error(hd, *ib));
	 case SORT_AMPLITUDE:
	      return cmp(cabs(harminv_get_amplitude(hd, *ia)), 
			 cabs(harminv_get_amplitude(hd, *ib)));
	 case SORT_Q:
	      return cmp(harminv_get_freq(hd,*ia) / harminv_get_decay(hd,*ia),
			 harminv_get_freq(hd,*ib) / harminv_get_decay(hd,*ib));
     }
     return 0;
}

typedef struct {
     int verbose;
     double fmin, fmax;
     int only_f_inrange;
     double err_thresh, rel_err_thresh, amp_thresh, rel_amp_thresh, Q_thresh;
     double min_err, max_amp;
     int num_ok;
} mode_ok_data;

static int mode_ok(harminv_data d, int k, void *ok_d_)
{
     mode_ok_data *ok_d = (mode_ok_data *) ok_d_;
     double errk, ampk, f;
     int ok;

     if (k == -1) { /* initialize */
	  int i;
	  ok_d->num_ok = 0;
	  if (!harminv_get_num_freqs(d))
	       return 0;
	  ok_d->min_err = harminv_get_freq_error(d, 0);;
	  ok_d->max_amp = cabs(harminv_get_amplitude(d, 0));
	  for (i = 1; i < harminv_get_num_freqs(d); ++i) {
	       double err, amp;
	       if ((err = harminv_get_freq_error(d, i)) < ok_d->min_err)
		    ok_d->min_err = err;
	       if ((amp = cabs(harminv_get_amplitude(d, i))) > ok_d->max_amp)
		    ok_d->max_amp = amp;
	  
	  }
	  return 0;
     }
     else if (k == -2) { /* finish */
	  if (ok_d->verbose && harminv_get_num_freqs(d))
	       printf("# harminv: %d/%d modes are ok: "
		      "errs <= %e and %e * %e\n, "
		      "amps >= %g, %e * %g, "
		      "|Q| >= %g\n", 
		      ok_d->num_ok, harminv_get_num_freqs(d),
		      ok_d->err_thresh, ok_d->rel_err_thresh, ok_d->min_err,
		      ok_d->amp_thresh, ok_d->rel_amp_thresh, ok_d->max_amp,
		      ok_d->Q_thresh);
	  return 0;
     }

     f = fabs(harminv_get_freq(d, k));
     errk = harminv_get_freq_error(d, k);
     ampk = cabs(harminv_get_amplitude(d, k));

     ok = ((!ok_d->only_f_inrange || (f >= ok_d->fmin && f <= ok_d->fmax))
	   && errk <= ok_d->err_thresh
	   && errk <= ok_d->min_err * ok_d->rel_err_thresh
	   && ampk >= ok_d->amp_thresh
	   && ampk >= ok_d->rel_amp_thresh * ok_d->max_amp
	   && fabs(harminv_get_Q(d,k)) >= ok_d->Q_thresh);

     ok_d->num_ok += ok;

     return ok;
}

#define SOLVE_ONCE_ONLY 0 /* 1 to use harminv_solve_once */
#define SOLVE_OK_ONLY 0 /* 1 for experimental solver */

int main(int argc, char **argv)
{
     int verbose = 0;
     int c;
     extern char *optarg;
     extern int optind;
     int specify_periods = 0;
     int specify_omega = 0;
     double dt = 1.0;
     mode_ok_data ok_d;
     int n, nf = NF;
     int iarg;
     cmplx *data;

     ok_d.only_f_inrange = 0;
     ok_d.err_thresh = ERR_THRESH;
     ok_d.rel_err_thresh = REL_ERR_THRESH;
     ok_d.amp_thresh = AMP_THRESH;
     ok_d.rel_amp_thresh = REL_AMP_THRESH;
     ok_d.Q_thresh = Q_THRESH;

     while ((c = getopt(argc, argv, "hvVTFwt:f:s:e:E:a:A:Q:")) != -1)
	  switch (c) {
	      case 'h':
		   usage(stdout);
		   return EXIT_SUCCESS;
	      case 'V':
		   printf("harminv " PACKAGE_VERSION " by Steven G. Johnson\n"
			  COPYRIGHT);
		   return EXIT_SUCCESS;
	      case 'v':
		   verbose = 1;
		   break;
	      case 'T':
		   specify_periods = 1;
		   break;
	      case 'w':
		   specify_omega = 1;
		   break;
	      case 'F':
		   ok_d.only_f_inrange = 1;
		   break;
	      case 'a':
		   ok_d.rel_amp_thresh = atof(optarg);
		   break;
	      case 'A':
		   ok_d.amp_thresh = atof(optarg);
		   break;
	      case 'E':
		   ok_d.err_thresh = atof(optarg);
		   break;
	      case 'e':
		   ok_d.rel_err_thresh = atof(optarg);
		   break;
	      case 'Q':
		   ok_d.Q_thresh = atof(optarg);
		   break;
	      case 't':
		   dt = atof(optarg);
		   break;
	      case 'f':
		   nf = atoi(optarg);
		   if (nf < 2) {
			fprintf(stderr, "harminv: -f argument must be > 1\n");
			return EXIT_FAILURE;
		   }
		   break;
	      case 's':
		   switch (tolower(optarg[0])) {
		       case 'f':
			    sortby = SORT_FREQUENCY;
			    break;
		       case 'd':
			    sortby = SORT_DECAY;
			    break;
		       case 'e':
			    sortby = SORT_ERROR;
			    break;
		       case 'a':
			    sortby = SORT_AMPLITUDE;
			    break;
		       case 'q':
			    sortby = SORT_Q;
			    break;
		       default:
			    fprintf(stderr, "harminv: invalid sort type -s %c\n", tolower(optarg[0]));
			    usage(stderr);
			    return EXIT_FAILURE;
		   }
		   break;
	      default:
		   fprintf(stderr, "harminv: invalid argument -%c\n", c);
		   usage(stderr);
		   return EXIT_FAILURE;
	  }
     if (optind == argc) {  /* no parameters left */
	  fprintf(stderr, "harminv: missing required frequency range(s)\n");
          usage(stderr);
          return EXIT_FAILURE;
     }	  

     data = read_input_data(stdin, &n, verbose);

     if (n < 1) {
	  fprintf(stderr, "harminv: no data read\n");
	  return EXIT_FAILURE;
     }

     if (verbose)
	  printf("# harminv: %d inputs, dt = %g, nf = %d\n", n, dt, nf);

     printf("frequency, decay constant, Q, amplitude, phase, error\n");

     ok_d.verbose = verbose;
	  
     for (iarg = optind; iarg < argc; ++iarg) {
	  double fmin, fmax;
	  int i;
	  int *isort = NULL;

	  if (sscanf(argv[iarg], "%lf-%lf", &fmin, &fmax) != 2) {
	       fprintf(stderr, "harminv: invalid argument \"%s\"\n",
		       argv[iarg]);
	       return EXIT_FAILURE;
	  }
	  if (specify_periods) {
	       if (fmin == 0 || fmax == 0) {
		    fprintf(stderr, "harminv: invalid argument \"%s\""
			    ": 0 not a valid period\n", argv[iarg]);
		    return EXIT_FAILURE;
	       }
	       fmin = 1/fmin;
	       fmax = 1/fmax;
	  }
	  if (specify_omega) {
	       fmin /= TWOPI;
	       fmax /= TWOPI;
	  }
	  if ((fmin > fmax && dt > 0) || (fmin < fmax && dt < 0)) {
	       double dummy = fmin;
	       fmin = fmax;
	       fmax = dummy;
	  }
	  if (verbose)
	       printf("# searching frequency range %g - %g\n", fmin, fmax);

	  ok_d.fmin = fmin*dt;
	  ok_d.fmax = fmax*dt;

	  hd = harminv_data_create(n, data, fmin*dt, fmax*dt, nf);
	  
#if SOLVE_OK_ONLY
	  harminv_solve_ok_modes(hd, mode_ok, &ok_d);
#elif SOLVE_ONCE_ONLY
	  harminv_solve_once(hd);
#else
	  harminv_solve(hd);
#endif

	  mode_ok(hd, -1, &ok_d); /* initialize ok_d */

	  CHK_MALLOC(isort, int, harminv_get_num_freqs(hd));
	  for (i = 0; i < harminv_get_num_freqs(hd); ++i) 
	       isort[i] = i;
	  qsort(isort, harminv_get_num_freqs(hd), sizeof(int), compar);

	  for (i = 0; i < harminv_get_num_freqs(hd); ++i) {
	       double freq, decay, err;
	       cmplx amp;
	       int j = isort[i];
#if SOLVE_OK_ONLY
	       CHECK(mode_ok(hd, j, &ok_d), "bug: invalid mode");
#else
	       if (!mode_ok(hd, j, &ok_d))
		    continue;
#endif
	       freq = harminv_get_freq(hd, j) / dt;
	       decay = harminv_get_decay(hd, j) / fabs(dt);
	       amp = harminv_get_amplitude(hd, j);
	       err = harminv_get_freq_error(hd, j);
	       printf("%g, %e, %g, %g, %g, %e\n",
		      freq * (specify_omega ? TWOPI : 1.0), decay,
		      harminv_get_Q(hd, j),
		      cabs(amp), -carg(amp), err);
	  }

#if !SOLVE_OK_ONLY
	  mode_ok(hd, -2, &ok_d);
#endif

	  harminv_data_destroy(hd);
     }

     free(data);
     return EXIT_SUCCESS;
}

#ifdef F77_DUMMY_MAIN
#  ifdef __cplusplus
extern "C"
#  endif
int F77_DUMMY_MAIN() { return 1; }
#endif

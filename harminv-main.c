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
   echo_comments != 0, then echo #... comments to stdout. */
int eat_whitespace(FILE *f, int echo_comments)
{
     int c, newlines = 0;
     do {
	  do {
	       c = getc(f);
	       newlines += c == '\n';
	  } while (isspace(c));
	  
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

cmplx *read_input_data(FILE *f, int *n, int verbose)
{
    cmplx *data = NULL;
    int line = 1, n_alloc = 0;
    *n = 0;

     do {
	  double re, im;
	  int nread;

	  line += eat_whitespace(f, verbose);
	  nread = fscanf(f, "%lg+%lgi", &re, &im);
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
	       data[*n] = re;
	       if (nread == 2)
		    data[*n] += I*im;
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
#define REL_AMP_THRESH 0.0


void usage(FILE *f)
{
     fprintf(f, "Usage: harminv [options] <freq-min>-<freq-max>...\n"
	     "Options: \n"
	     "         -h : this help message\n"
             "         -V : print version number and copyright\n"
             "         -v : verbose output\n"
	     "         -T : specify periods instead of frequencies\n"
	     "    -t <dt> : specify sampling interval dt [default: 1]\n"
	     "    -f <nf> : specify initial spectral density [default: %d]\n"
	     "  -s <sort> : sort by <sort> = freq/err/decay/amp [default: freq]\n"
	     "    -a <a> : discard amplitudes < max * <a> [default: %e]\n"
	     "    -A <A> : discard amplitudes < <A> [default: %g]\n"
	     "    -e <e> : discard relative errors > min * <e> [default: %e]\n"
	     "    -E <E> : discard relative errors > <E> [default: %e]\n"
,
	     NF,
	     AMP_THRESH, REL_AMP_THRESH,
	     ERR_THRESH, REL_ERR_THRESH);
}

#define TWOPI 6.2831853071795864769252867665590057683943388

harminv_data hd;
cmplx *amps = NULL;
double *errs = NULL;

enum {
     SORT_FREQUENCY, SORT_DECAY, SORT_ERROR, SORT_AMPLITUDE, SORT_Q
} sortby = SORT_FREQUENCY;

int cmp(double a, double b)
{
     return a > b ? 1 : (a < b ? -1 : 0);
}

int compar(const void *a, const void *b)
{
     const int *ia = (const int *) a;
     const int *ib = (const int *) b;

     switch (sortby) {
	 case SORT_FREQUENCY:
	      return cmp(harminv_get_freq(hd,*ia), harminv_get_freq(hd,*ib));
	 case SORT_DECAY:
	      return cmp(harminv_get_decay(hd,*ia), harminv_get_decay(hd,*ib));
	 case SORT_ERROR:
	      return cmp(errs[*ia], errs[*ib]);
	 case SORT_AMPLITUDE:
	      return cmp(cabs(amps[*ia]), cabs(amps[*ib]));
	 case SORT_Q:
	      return cmp(harminv_get_freq(hd,*ia) / harminv_get_decay(hd,*ia),
			 harminv_get_freq(hd,*ib) / harminv_get_decay(hd,*ib));
     }
     return 0;
}


int main(int argc, char **argv)
{
     int c;
     extern char *optarg;
     extern int optind;
     int verbose = 0;
     int specify_periods = 0;
     double dt = 1.0;
     double err_thresh = ERR_THRESH;
     double rel_err_thresh = REL_ERR_THRESH;
     double amp_thresh = AMP_THRESH;
     double rel_amp_thresh = REL_AMP_THRESH;
     int n, nf = NF;
     int iarg;
     cmplx *data;

     while ((c = getopt(argc, argv, "hvVTt:f:s:e:E:a:")) != -1)
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
	      case 'a':
		   rel_amp_thresh = atof(optarg);
		   break;
	      case 'A':
		   amp_thresh = atof(optarg);
		   break;
	      case 'E':
		   err_thresh = atof(optarg);
		   break;
	      case 'e':
		   rel_err_thresh = atof(optarg);
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

     for (iarg = optind; iarg < argc; ++iarg) {
	  double fmin, fmax;
	  double min_err = 1e20, max_amp = 0;
	  int i, cur_nf, prev_nf;
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
	  if ((fmin > fmax && dt > 0) || (fmin < fmax && dt < 0)) {
	       double dummy = fmin;
	       fmin = fmax;
	       fmax = dummy;
	  }
	  if (verbose)
	       printf("# searching frequency range %g - %g\n", fmin, fmax);

	  hd = harminv_data_create(n, data, fmin*dt, fmax*dt, nf);
	  
	  harminv_solve_once(hd);
	  prev_nf = cur_nf = harminv_get_num_freqs(hd);

	  /* keep re-solving as long as spurious solutions are eliminated */
	  if (verbose) {
	       printf("# harminv number of solution frequencies = %d", cur_nf);
	       fflush(stdout);
	  }
	  do {
	       prev_nf = cur_nf;
	       harminv_solve_again(hd);
	       cur_nf = harminv_get_num_freqs(hd);
	       if (verbose) {
		    printf(", then %d", cur_nf);
		    fflush(stdout);
	       }
	  } while (cur_nf < prev_nf);
	  if (verbose)
	       printf("\n");

	  if (cur_nf > prev_nf)
	       fprintf(stderr, "harminv: warning, number of solutions increased from %d to %d!\n", prev_nf, cur_nf);

	  errs = harminv_compute_frequency_errors(hd);
	  amps = harminv_compute_amplitudes(hd);

	  if (harminv_get_num_freqs(hd)) min_err = errs[0];
	  for (i = 0; i < harminv_get_num_freqs(hd); ++i) {
	       if (errs[i] < min_err) min_err = errs[i];
	       if (cabs(amps[i]) > max_amp) max_amp = cabs(amps[i]);
	  }
	  if (verbose) {
	       cur_nf = 0;
	       for (i = 0; i < harminv_get_num_freqs(hd); ++i)
		    if (errs[i] / fabs(dt) <= err_thresh
			&& errs[i] <= min_err * rel_err_thresh
			&& cabs(amps[i]) >= amp_thresh
			&& cabs(amps[i]) >= rel_amp_thresh * max_amp)
			 ++cur_nf;
	       printf("# harminv number of frequencies with err < %e and "
		      "< %e * %e, amp > %g and %e * %g: %d\n",
		      err_thresh, rel_err_thresh, min_err / fabs(dt),
		      amp_thresh, rel_amp_thresh, max_amp, cur_nf);
	  }
	       

	  CHK_MALLOC(isort, int, harminv_get_num_freqs(hd));
	  for (i = 0; i < harminv_get_num_freqs(hd); ++i) 
	       isort[i] = i;
	  qsort(isort, harminv_get_num_freqs(hd), sizeof(int), compar);
	  
	  for (i = 0; i < harminv_get_num_freqs(hd); ++i) {
	       double freq, decay;
	       int j = isort[i];
	       if (errs[j] / fabs(dt) > err_thresh
		   || errs[j] > min_err * rel_err_thresh
		   || cabs(amps[j]) < amp_thresh
		   || cabs(amps[j]) < rel_amp_thresh * max_amp)
		    continue;
	       freq = harminv_get_freq(hd, j) / dt;
	       decay = harminv_get_decay(hd, j) / fabs(dt);
	       printf("%g, %e, %g, %g, %g, %e\n",
		      freq, decay, TWOPI * fabs(freq) / (2 * decay),
		      cabs(amps[j]), carg(amps[j]), errs[j] / fabs(dt));
	  }

	  free(amps);
	  free(errs);

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

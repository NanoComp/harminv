/* Copyright (C) 2000 Massachusetts Institute of Technology.
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

#define NF 100

void usage(FILE *f)
{
     fprintf(f, "Usage: harminv [options] <freq-min>-<freq-max>...\n"
	     "Options: \n"
	     "         -h : this help message\n"
             "         -V : print version number and copyright\n"
             "         -v : verbose output\n"
	     "         -T : specify periods instead of frequencies\n"
	     "    -t <dt> : specify sampling interval dt [default: 1]\n"
	     "    -f <nf> : specify initial spectral density [default: %d]\n",
	     NF);
}

#define TWOPI 6.2831853071795864769252867665590057683943388

int main(int argc, char **argv)
{
     int c;
     extern char *optarg;
     extern int optind;
     int verbose = 0;
     int specify_periods = 0;
     double dt = 1.0;
     int n, nf = NF;
     int iarg;
     cmplx *data;

     while ((c = getopt(argc, argv, "hvVTt:f:")) != -1)
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
	  int i, cur_nf, prev_nf;
	  harminv_data hd;
	  cmplx *amps = NULL;
	  double *errs = NULL;

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
	  
	  for (i = 0; i < harminv_get_num_freqs(hd); ++i) {
	       double freq, decay;
	       freq = harminv_get_freq(hd, i) / dt;
	       decay = harminv_get_decay(hd, i) / fabs(dt);
	       printf("%g, %g, %g, %g, %g, %g\n",
		      freq, decay, TWOPI * fabs(freq) / (2 * decay),
		      cabs(amps[i]), carg(amps[i]), errs[i] / fabs(dt));
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

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
#include <time.h>
#include <math.h>

#include "config.h"
#include "check.h"
#include "harminv.h"
#include "copyright.h"

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif
#ifdef HAVE_GETOPT_H
#  include <getopt.h>
#endif

int verbose = 0;

#define TWOPI 6.2831853071795864769252867665590057683943388
#define NPERIODS 10

void usage(FILE *f)
{
     fprintf(f, "Usage: sines [options] <freq>...\n"
	     "Note that Im[freq] is the decay rate, or its inverse for -T.\n"
	     "Options: \n"
	     "         -h : this help message\n"
             "         -V : print version number and copyright\n"
             "         -v : verbose output\n"
	     "         -T : specify periods instead of frequencies\n"
	     "         -r : random amplitudes\n"
	     "     -n <n> : output <n> points (default %d * max period)\n"
	     "    -t <dt> : time step <dt> (default 1.0)\n",
	     NPERIODS);
}

typedef struct {
     double freq, decay, amplitude, phase;
} sinusoid;

#define MAX2(a,b) ((a) > (b) ? (a) : (b))

int main(int argc, char **argv)
{
     int c;
     extern char *optarg;
     extern int optind;
     int iarg;
     int specify_periods = 0;
     int random_amplitudes = 0;
     sinusoid *sines = NULL;
     int nsines = 0, nalloc = 0, n = 0;
     double max_period = 0;
     double dt = 1.0;
     int i, is;

     srand(time(NULL));

     while ((c = getopt(argc, argv, "hVvTrn:t:")) != -1)
	  switch (c) {
	      case 'h':
		   usage(stdout);
		   return EXIT_SUCCESS;
	      case 'V':
		   printf("sines " VERSION " by Steven G. Johnson.\n"
			  "Test program for harminv.\n"
			  COPYRIGHT);
		   return EXIT_SUCCESS;
	      case 'v':
		   verbose = 1;
		   break;
	      case 'T':
		   specify_periods = 1;
		   break;
	      case 'r':
		   random_amplitudes = 1;
		   break;
	      case 'n':
		   n = atoi(optarg);
		   if (n < 1) {
			fprintf(stderr, "sines: "
				"invalid non-positive -n argument %d\n", n);
			return EXIT_FAILURE;
		   }
		   break;
	      case 't':
		   dt = atof(optarg);
		   break;
	      default:
		   fprintf(stderr, "sines: invalid argument -%c\n", c);
		   usage(stderr);
		   return EXIT_FAILURE;
	  }

     if (optind == argc) {  /* no parameters left */
	  fprintf(stderr, "sines: missing required frequency(ies)\n");
          usage(stderr);
          return EXIT_FAILURE;
     }	  

     for (iarg = optind; iarg < argc; ++iarg) {
	  sinusoid s = { 0, 0, 0, 0 };

	  if (sscanf(argv[iarg], "%lf+%lfi", &s.freq, &s.decay) < 1) {
	       fprintf(stderr, "sines: invalid argument \"%s\"\n",
		       argv[iarg]);
	       return EXIT_FAILURE;
	  }
	  if (specify_periods) {
	       if (s.freq == 0) {
		    fprintf(stderr, "sines: invalid argument \"%s\""
			    ": 0 not a valid period\n", argv[iarg]);
		    return EXIT_FAILURE;
	       }
	       s.freq = 1/s.freq;
	       if (s.decay != 0)
		    s.decay = 1/s.decay;
	  }

	  if (s.decay == 0 || fabs(1/s.freq) > 1/s.decay)
	       max_period = MAX2(max_period, fabs(1/s.freq));
	  else
	       max_period = MAX2(max_period, 1/s.decay);

	  if (random_amplitudes) {
	       s.amplitude = rand() * 1.0/RAND_MAX;
	       s.phase = rand() * TWOPI/RAND_MAX - TWOPI/2;
	  }
	  else {
	       s.amplitude = (iarg - optind + 1);
	       s.phase = (iarg - optind + 1) * TWOPI / (argc-optind) - TWOPI/2;
	  }

	  if (verbose)
	       printf("# mode: frequency = %g (period %g), decay = %g (lifetime %g), amplitude = %g, phase = %g\n", s.freq, 1/s.freq, s.decay, s.decay != 0 ? 1/s.decay : 0, s.amplitude, s.phase);

	  if (nsines >= nalloc) {
	       nalloc = (nalloc + 1) * 2;
	       sines = (sinusoid*) realloc(sines, sizeof(sinusoid) * nalloc);
	       CHECK(sines, "out of memory");
	  }
	  sines[nsines++] = s;
     }

     if (n == 0 && max_period == 0) {
	  fprintf(stderr, "sines: must specify -n or non-zero period\n");
	  return EXIT_FAILURE;
     }

     if (n == 0)
	  n = (int) (max_period / fabs(dt) + 0.5) * NPERIODS;

     for (i = 0; i < n; ++i) {
	  cmplx output = 0;

	  for (is = 0; is < nsines; ++is)
	       output += sines[is].amplitude *
		    cexp(I * (sines[is].phase - TWOPI*sines[is].freq * i*dt)
			 - sines[is].decay * i*dt);
	  printf("%0.17e+%0.17ei\n", creal(output), cimag(output));
     }

     free(sines);
     return EXIT_SUCCESS;
}

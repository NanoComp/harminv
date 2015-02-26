## Harminv 1.4

26 February 2015

* Change `get_omega` and `get_amplitude` functions to return complex results
  via a pointer argument rather than a return value, since the latter
  is not generally binary compatible between C++ and C.

* Build script updates.

## Harminv 1.3.1

4 July 2006

* Fixed the phase output column to have the same sign as the documentation
  (the sign was flipped).  Thanks to Andrew Norton for the bug report.

* Added `-n` option to flip the sign of the frequency convention.

## Harminv 1.3:

18 October 2005

* Switch back to eigensolver technique used in Harminv 1.0.x, which first
  removes the singular null space as described by Wall and Neuhauser.
  This seems to make the solution much more stable and reliable.  The
  command-line tool once again defaults to 100 basis modes rather than
  a particular spectral density.

## Harminv 1.2.1:

20 May 2004

* Impose a maximum number of basis modes (300) to prevent the matrices
  from getting too large.

* Corrected typo in man page (for definition of `-d` density).

## Harminv 1.2:

19 May 2004

* Command line tool now defaults to a particular spectral "density"
  (set by the `-d` option), rather than a number of basis modes,
  since using too many basis modes leads to a singular eigenproblem
  and numerical instability.  (Based on defaults from M&T references.)

* Use `long double` precision, if available, to reduce accumulation
  of floating-point errors while computing Fourier/Z transforms.

## Harminv 1.1:

18 May 2004

* Corrected bug in frequency-error calculation; thanks to V. A.
  Mandelshtam for helpful discussions and for letting me look
  at his code to check against mine.

* Amplitude calculation is no longer unstable for strongly decaying modes.

* Used a more accurate eigensolver routine.

* Added `-e`/`-E`/`-a`/`-A`/`-Q` options to screen outputs with
  error/amplitude/Q too large/small/small, respectively.

* Added `-w` option to use angular frequency instead of frequency.

* More flexible input format: allow comma-delimited, a-bi as well as a+bi.

* API cleanups.

## Harminv 1.0.2:

15 May 2004

* Corrected inadvertent windowing of data that degraded accuracy
  in the case of very short signals.

## Harminv 1.0.1:

15 May 2004

* Corrected some minor release glitches.

## Harminv 1.0:

14 May 2004

* Initial release (after 4 years of private use).

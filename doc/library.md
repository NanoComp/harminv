[Harminv](../README.md) is mostly used via the stand-alone
`harminv` program, but it can also be called as a library
from C or C++, as described below.

## Library Usage

The usage of the library `-lharminv` is analogous to the program.  In C
or C++, you first `#include <harminv.h>`, then specify the data and the
frequency range by calling `harminv_data_create`, returning a
`harminv_data` data structure:
```c
harminv_data harminv_data_create(int n,
                                 const harminv_complex *signal,
                                 double fmin, double fmax, int nf);
```
Here, `signal` is a pointer to an array of `n` complex numbers.  In C++,
`harminv_complex` is `std::complex<double>`.  In C, `harminv_complex` is a
`double[2]` with the real parts in `signal[i][0]` and the imaginary parts
in `signal[i][1]`.  (For a real signal, set the imaginary parts to
zero.)  `fmin` and `fmax` are the frequency range to search, and `nf` is the
number of spectral basis functions (see below).  Frequencies are in
units corresponding to a sampling interval of 1 time unit; if your
actual sampling interval is dt, then you should rescale your
frequencies by multiplying them by dt.

A good default for `nf` is `min(300, (fmax - fmin) * n * 1.1)`,
corresponding to a spectral "density" of at most 1.1 (see also the `-d`
option of the command-line tool).  That is, this uses a number of
initial basis functions corresponding to the Fourier resolution of
`1/n`.  This does *not* determine the frequency resolution of the
outputs, which can be much greater than the Fourier resolution.  It
sets an upper bound on the number of modes to search for, and in some
sense is the density with which the bandwidth is initially "searched"
for modes.  Spectral densities much larger than 1 are not recommended,
as they lead to large and singular matrices and unstable results.
Note also that the computation time goes as O(n * nf) + O(nf^3).

Then, you solve for the frequencies by calling:
```c
void harminv_solve(harminv_data d);
```
Then, the frequencies and other data can be extracted from `d` by the
following routines.  The number N of frequencies found is returned by:
```c
int harminv_get_num_freqs(harminv_data d);
```
Then, for each index 0 <= k < N, the corresponding frequency and decay
constant (as defined in `man harminv`) are returned by:
```c
double harminv_get_freq(harminv_data d, int k);
double harminv_get_decay(harminv_data d, int k);
```
Alternative, you can get the complex angular frequency (omega =
2π freq - i decay) by:
```c
    void harminv_get_omega(harminv_complex *omega, harminv_data d, int k);
```
You can get the "quality factor" Q (pi |freq| / decay) by:
```c
double harminv_get_Q(harminv_data d, int k);
```
The complex amplitude (|amp| * exp(-I phase)) for each k is returned by:
```c
void harminv_get_amplitude(harminv_complex *amplitude, harminv_data d, int k);
```
A crude estimate of the relative error in the (complex) frequency is:
```c
double harminv_get_freq_error(harminv_data d, int k);
```
As described in `man harminv`, this is not really an error bar, and
should be treated more as a figure of merit (smaller is better).

### Linking

To link to the library, you need to not only link to `-lharminv`, but
also to the math library, the BLAS and LAPACK libraries (see below),
and any libraries that are required to link C with Fortran code (like
LAPACK).  If you have the `pkg-config` program installed (standard on most
GNU/Linux systems), you can simply do:
```
pkg-config --cflags harminv
pkg-config --libs harminv
```
to get the flags for compiling and linking, respectively.  You may
need to tell `pkg-config` where to find `harminv.pc` if `harminv` was
installed under `/usr/local` (the default)...in this case, you would
specify `/usr/local/lib/pkgconfig/harminv.pc` instead of `harminv`
above.

There is an additional wrinkle.  If you configured harminv with
`--with-cxx`, or if your C compiler did not support C99 complex numbers
and the configure script automatically switched to C++, then you will
need to link to harminv with the C++ linker, even if your program is
written in C, in order to link the C++ libraries.

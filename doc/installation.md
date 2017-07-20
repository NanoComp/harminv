[Harminv](../README.md) is designed to run on any Unix-like system
(GNU/Linux is fine), and uses a `configure` script to make it easy to
install. However, you do need a couple of prerequisites:

# Prerequisites

However, you do need a couple of prerequisites:

* [BLAS](http://www.netlib.org/blas/)

Basic Linear Algebra Subroutines (matrix-multiplications, etcetera),
following a standard interface, for use by LAPACK (see below).  There
are many optimized versions of BLAS available, e.g. a free library
called [ATLAS](http://www.netlib.org/atlas/) or [OpenBLAS](http://www.openblas.net/).

* [LAPACK](http://www.netlib.org/lapack/)

A standard, free, linear-algebra package.  Note that the default
configuration script looks for LAPACK by linking with `-llapack`.  This
means that the library must be called `liblapack.a` and be installed in
a standard directory like `/usr/local/lib` (alternatively, you can
specify another directory via the `LDFLAGS` environment variable; see
below).  

Often, you will install an optimized implementation of BLAS;
an excellent choice is the free [OpenBLAS library](http://www.openblas.net/), which also includes LAPACK.

# Compiling from Git

Most users should download an official Harminv release (a prepackaged `.tar.gz` file).  *If* you are installing from the raw `git` repository (rather than downloading a `.tar.gz` release), then you will also need

* [GNU Autotools](https://en.wikipedia.org/wiki/GNU_Build_System): `autoconf`, `automake`, and `libtool`.

and run
```
sh autogen.sh
make
```
where `autogen.sh` is a script included with Harminv that
re-generates the `configure` script and other necessary files.

# Compiling Harminv

Given the above, you can compile and install harminv.
Harminv comes with a [GNU-style](http://www.gnu.org/software/autoconf/)
`configure` script, so on Unix-like systems compilation is ideally just
a matter of:
```
./configure
make
```
and then switching to root and running:
```
make install
```
By default, this installs under `/usr/local`, i.e. in `/usr/local/bin`
etcetera. You can change this by passing the standard `--prefix=`*dir*
option to `configure`. Other `configure` options can be found by running
`./configure --help`.

In order to compile, harminv requires either:

-   An ANSI C compiler supporting complex numbers, as defined in the
    ANSI C99 standard (or a reasonable approximation thereof). For
    example, gcc-2.95 with GNU libc is fine.
-   A C++ compiler supporting the complex standard template class.

The `configure` script looks for a C compiler with complex numbers
first, and then, if that fails, for a C++ compiler. You can force it to
use C++ by passing `--with-cxx` to `configure`.

If you need to, you can have further control over the `configure`
script's behavior by setting enviroment variables, by passing
`VARIABLE=VALUE` arguments to `configure`. This can be useful
especially if you have libraries installed in nonstandard locations
(e.g. in your home directory, if you are not a system administrator), to
tell the compiler where to look. The most common variables to set are:

-   `CC`: the C compiler command
-   `CFLAGS`: the C compiler flags
-   `CXX`: the C++ compiler command
-   `CXXFLAGS`: the C++ compiler flags
-   `F77`: the Fortran 77 compiler command. **Important:** if you have
    more than one Fortran compiler, use the same compiler here as you
    used for BLAS/LAPACK.
-   `FFLAGS`: the Fortran 77 compiler flags
-   `CPPFLAGS`: `-I''dir''` flags to tell the C compiler additional
    places to look for header files.
-   `LDFLAGS`: `-L''dir''` flags to tell the linker additional places to
    look for libraries.
-   `LIBS`: additional libraries to link against.

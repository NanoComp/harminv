# harminv command-line program

[Harminv](../README.md) is installed as both a [library](library.md)
and as a command-line program `harminv`.   This page
describes the usage of the command-line program.

## Synopsis

    harminv [OPTION]... [freq-min-freq-max]...

## Description

`harminv` is a program designed to solve the problem of "harmonic inversion": given a time series consisting of a sum of sinusoids ("modes"), extract their frequencies and amplitudes. It can also handle the case of exponentially-decaying sinusoids, in which case it extracts their decay rates as well.

`harminv` is often able to achieve much greater accuracy and robustness than Fourier-transform methods, essentially because it assumes a specific form for the input.

It uses a low-storage "filter-diagonalization method" (FDM), as described in V. A. Mandelshtam and H. S. Taylor, "Harmonic inversion of time signals," *J. Chem. Phys.*, vol. 107, p. 6756 (1997). See also erratum, ibid `109`, 4128 (1998).

## Input

`harminv` reads in a sequence of whitespace-separated real or complex numbers from standard input, as well as command-line arguments indicating one or more frequency ranges to search, and outputs the modes that it extracts from the data. (It preferentially finds modes in the frequency range you specify, but may sometimes find additional modes outside of that range.) The data should correspond to equally-spaced time intervals, but there is no constraint on the number of points.

Complex numbers in the input should be expressed in the format `RE+IMi` (*no whitespace*). Otherwise, whitespace is ignored. Also, comments beginning with `#` and extending to the end of the line are ignored.

A typical invocation is something like

```sh
harminv -t 0.02 1-5 < input.dat
```

which reads a sequence of samples, spaced at 0.02 time intervals (in ms, say, corresponding to 50 kHz), and searches for modes in the frequency range 1–5 kHz. (See below on units.)

## Output

`harminv` writes four comma-delimited columns to standard output, one line for each mode: frequency, decay constant, Q, amplitude, phase, and error. Each mode corresponds to a function of the form:

    amplitude * exp[-i (2 pi frequency t - phase) - decay t]

Here, *i* is sqrt(-1), *t* is the time (see below for units), and the other parameters in the output columns are:

 * `frequency` — The frequency of the mode. If you don't recognize that from the expression above, you should recall Euler's formula: `exp(ix) = cos(x) + i sin(x)`. Note that for complex data, there is a distinction between positive and negative frequencies.
 * `decay constant` — The exponential decay constant, indicated by decay in the above formula. The inverse of this is often called the "lifetime" of the mode. The "half-life" is `ln(2)/decay`.
 * `Q` — A conventional, dimensionless expression of the decay lifetime: `Q = π |frequency| / decay`. Q, which stands for "quality factor", is the number of periods for the "energy" in the mode (the squared amplitude) to decay by `exp(-2π)`. Equivalently, if you look at the power spectrum (|Fourier transform|²), 1/Q is the fractional width of the peak at half maximum.
 * `amplitude` — The (real, positive) amplitude of the sinusoids. The amplitude (and phase) information generally seems to be less accurate than the frequency and decay constant.
 * `phase` — The phase shift (in radians) of the sinusoids, as given by the formula above.
 * `error` — A crude estimate of the relative error in the (complex) frequency. This is not really an error bar, however, so you should treat it more as a figure of merit (smaller is better) for each mode.

## Spurious Modes

Typically, harminv will find a number of spurious solutions in addition to the desired solutions, especially if your data are noisy. Such solutions are characterized by large errors, small amplitudes, and/or small Q (large decay rates / broad linewidths). You can omit these from the output by the error/Q/amplitude screening options defined below.

By default, modes with `error > 0.1` and `Q < 10` are automatically omitted, but it is likely that you will need to set stricter limits.

## Units

The frequency (and decay) values, both input and output, are specified in units of 1/time, where the units of time are determined by the sampling interval *dt* (the time between consecutive inputs). *dt* is by default 1, unless you specify it with the `-t dt` option.

In other words, pick some units (e.g. ms in the example above) and use them to express the time step. Then, be consistent and use the inverse of those units (e.g. kHz = 1/ms) for frequency.

Note that the frequency is the usual 1/period definition; it is not the angular frequency.

## Options

 * `-h` — Display help on the command-line options and usage.
 * `-V` — Print the version number and copyright info for `harminv`.
 * `-v` — Enable verbose output, printed to standard output as comment lines (starting with a `#` character). Also, any `#` comments in the input are echoed to the output.
 * `-T` — Specify period-ranges instead of frequency-ranges on the command line (in units of time corresponding to those specified by `-t`). The output is still frequency and not period, however.
 * `-w` — Specify angular frequencies instead of frequencies, and output angular frequency instead of frequency. (Angular frequency is frequency multiplied by 2π).
 * `-t dt` — Specify the sampling interval `dt`; this determines the units of time used throughout the input and output. Defaults to 1.0.
 * `-d d` — Specify the spectral "density" d to search for modes, where a density of 1 indicates the usual Fourier resolution. That is, the number of basis functions (which sets an upper bound on the number of modes) is given by `d × (freq-max - freq-min) × dt × (number of samples)`. A maximum of 300 is used, however, to prevent the matrices from getting too big (you can force a larger number with `-f`, below).
   - Note that the frequency resolution of the outputs is not limited by the spectral density, and can generally be much greater than the Fourier resolution. The density determines how many modes, at most, to search for, and in some sense is the density with which the bandwidth is initially "searched" for modes.
   - The default density is 1.1 (or lower, to keep within the 300 maximum), which is usually a good value for most applications. If you set the density too high, then the matrices become large and singular; if you set the density too low, then you risk missing modes.
 * `-f nf` — Specify a lower bound nf on the number of spectral basis functions (defaults to 2), setting a lower bound on the number of modes to search for. This option is sometimes a more convenient way to specify the number of basis functions than the `-d` option, above.
   - `-f` also allows you to employ more than 300 basis functions, but careful: the computation time scales as O(N nf) + O(nf^3), where N is the number of samples, and very large matrices can also have degraded accuracy.
 * `-s sort` — Specify how the outputs are sorted, where `sort` is one of `frequency`/`error`/`Q`/`decay`/`amplitude`. (Only the first character of `sort` matters, e.g. `-s a` sorts by amplitude.) All sorts are in ascending order. The default is to sort by frequency (`-s f`).
 * `-e err` — Omit any modes with error (see above) greater than err times the largest error among the computed modes. Defaults to no limit.
 * `-E err` — Omit any modes with error (see above) greater than err. Defaults to 0.1.
 * `-F` — Omit any modes with frequencies outside the specified range. (Such modes are not necessarily spurious, however.)
 * `-a amp` — Omit any modes with amplitude (see above) less than amp times the largest amplitude among the computed modes. Defaults to no limit.
 * `-A amp` — Omit any modes with amplitude (see above) less than amp. Defaults to no limit.
 * `-Q q` — Omit any modes with |Q| (see above) less than `q`. Defaults to 10.

## Bugs

Report bugs by [filing a github issue](https://github.com/stevengj/harminv).

## Authors

Written by Steven G. Johnson. Copyright © 2004–2017 by the Massachusetts Institute of Technology.

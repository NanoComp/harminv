Date: Sun, 16 May 2004 17:05:22 -0700 (PDT)
From: Vladimir Mandelshtam <mandelsh@uci.edu>
To: Steven G. Johnson <stevenj@fftw.org>
Subject: filter diagonalization code  (fwd)

Steven,
The message below is the standard message I send with these files.
This program does not include the latest developments (multi-window
with a multi-scale basis)
We have the corresponding code which I can also send you, but it was not
meant for distribution, so it may not be easy to figure out what it does.
Vladimir.

_________________________________________________


I am attaching to this message the filter-diagonalization code "fd_qz99.f"

to invert  c_n=sum_k d_k exp(-i*n*tau*w_k)

for unknown d_k and w_k

a typical input file "in_qz" and a fortran code to generate a test signal
with little noise, "generate.f". The signal is called "Jacob's ladder".
The noise is due to "integerization", i.e. some kind of roundoff.

Compile the file generate.f

        f77 -O generate.f -o generate.exe

then run
        generate.exe

which will create a file, "signal", with the model
signal and file "param", with spectral parameters.

Compile the FDM code

        f77 -O3 fd_qz.f -o fd_qz.exe

and then

        fd_qz.exe<in_qz

The absorption mode spectrum, i.e. Re \int C(t) exp(iwt) dt,
will be stored in file "ReSp" and the
spectral parameters in "parameters".

Let me know how it works. You can try to change the parameters such
as the size of the window, the number of basis functions in the window,
or introduce a time delay by using a non-zero Nskip and changing the t0
accordingly, etc.

Note, that this program was not written for distribution purpose.

Please, don't distribute it.

Regards,

Vladimir.


Vladimir Mandelshtam                           (949) 824-5509 (of)
Chemistry Department,                          (949) 824-8571 (FAX)
University of California at Irvine,            (949) 854-3243 (h)
Irvine, CA 92697                               e-mail: mandelsh@uci.edu
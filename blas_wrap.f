c     Wrappers around some functions from BLAS, changing them from
c     functions into subroutines.  This allows us to call them from C,
c     since function call protocols can sometimes differ between the two
c     languages.

c     This will have to be modified on a Cray vector machine, to use
c     'real' instead of 'double precision', 'complex' instead of
c     'double complex', and snrm2/cdotu, since on those machines
c     single precision *is* double precision.

      subroutine harminv_dnrm2(result, n, x, incx)
      implicit none
      integer n, incx
      double precision result, x(n)
      double precision dnrm2
      external dnrm2
      result = dnrm2(n, x, incx)
      end

      subroutine harminv_zdotu(result, n, x, incx, y, incy)
      implicit none
      integer n, incx, incy
      double complex result, x(n), y(n)
      double complex zdotu
      external zdotu
      result = zdotu(n, x, incx, y, incy)
      end

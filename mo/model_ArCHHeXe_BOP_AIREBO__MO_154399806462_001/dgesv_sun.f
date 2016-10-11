c-----------------------------------------------------------------------
c     dgesvw provides a wrapper for the dgesv routine.
c-----------------------------------------------------------------------
c     Sun supplies a version of dgesv with the Sun Performance Library.
c     ...sjs 11/27/07
c-----------------------------------------------------------------------
c
      subroutine dgesvw(n, nra, a, b, istatus)

      implicit none

c     n       = order of the matrix a
c     nra     = number of rows (leading dimension) of array a
c     a(i,j)  = i,j element of array a (on input, changed on output)
c     b(i)    = ith element of RHS vector on input; solution x on output
c     istatus = success (0 => successful, >0 => failure due to singular matrix)

      integer n
      real*8  a(nra,n)
      real*8  b(n)
      integer istatus

c     local variables

      integer nmax
      parameter (nmax = 10000)

      integer ipivot(nmax),
     .     nra, nrb, nrhs

      if (n .gt. nmax) then
         write(0, *) 'dgesvw: n = ', n, ' is larger than nmax = ', nmax
         write(0, *) '        increase nmax and recompile'
         stop
      endif

c     call Sun's dgesv

      nrhs = 1
      nrb = n

      call dgesv(n, nrhs, a, nra, ipivot, b, nrb, istatus)

      end
c

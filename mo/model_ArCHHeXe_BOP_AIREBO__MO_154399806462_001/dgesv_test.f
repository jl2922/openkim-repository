c-----------------------------------------------------------------------
c     wrapper program used to determine how to call dgesv().
c     ...sjs 11/27/07
c-----------------------------------------------------------------------
c
      program wrappr

      implicit none

c     local variables

      integer nvar
      parameter (nvar = 3)

      integer istatus, ivar, jvar, nra
      real*8 a(nvar,nvar),
     .     b(nvar), x(nvar),
     .     error

c     the correct answer, x

      do 100 ivar = 1, nvar
         x(ivar) = dble(ivar)
 100  continue

c     construct the matrix a and vector b = a x

      do 210 ivar = 1, nvar
         b(ivar) = 0.d0
         a(ivar,ivar) = 2.d0
         b(ivar) = b(ivar) + a(ivar,ivar) * x(ivar)
         do 200 jvar = 1, nvar
            if (jvar .ne. ivar) then
               a(ivar,jvar) = 1.d0
               b(ivar) = b(ivar) + a(ivar,jvar) * x(jvar)
            endif
 200     continue
 210  continue
               
c     solve x = a^-1 b by Gaussian elimination

      nra = nvar
      call dgesvw(nvar, nra, a, b, istatus)

      error = 0.d0
      do 300 ivar = 1, nvar
         error = error + dabs(b(ivar) - x(ivar))
 300  continue

      if (error .eq. 0.d0) then
         write(6, *) 'dgesv is working'
      else
         write(6, *) 'dgesv is not working. error = ', error
      endif

      end
c

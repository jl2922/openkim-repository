c-----------------------------------------------------------------------
c     flush provides a wrapper for whatever the native flush routine is,
c     if it is not called flush.
c-----------------------------------------------------------------------
c     With the Absoft compiler on Mac OS X,
c     the routine is called flush_, and is accessible via the libU77.a
c     library....sjs 11/14/05
c-----------------------------------------------------------------------
      subroutine flush(iunit)

      integer*4 iunit

      call flush_(iunit)

      return
      end
c

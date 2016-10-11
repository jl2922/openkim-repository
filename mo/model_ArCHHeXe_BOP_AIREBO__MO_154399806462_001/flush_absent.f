c-----------------------------------------------------------------------
c     flush provides a wrapper for whatever the native flush routine is,
c     if it is not called flush.
c-----------------------------------------------------------------------
c      For compilers than don't offer any version of the flush routine
c      (Digital Visual Fortran on Windows?) provide a dummy routine.
c      ...sjs 3/9/06
c-----------------------------------------------------------------------
      subroutine flush(iunit)

      integer*4 iunit

      return
      end
c

c-----------------------------------------------------------------------
c     getdir is a stub, designed to stand in for the wrapper-side 
c     routine, when called from codes that have no access to the wrapper
c     side.
c     ...sjs 10/2/09
c-----------------------------------------------------------------------
c
      subroutine getdir(dirnam)

      include 'common_files.inc'

      character dirnam*80

      dirnam = '.'

      return
      end
c

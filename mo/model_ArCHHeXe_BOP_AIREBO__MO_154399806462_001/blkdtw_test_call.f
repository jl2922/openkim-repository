c-----------------------------------------------------------------------
c     wrapper program used to determine whether an explicit call to the
c     block data subprogram is allowed....sjs 11/3/08
c-----------------------------------------------------------------------
c
      program wrappr

      real*8 x

      common/tstcom/x

      call blkdat

      if (x .eq. 1.234d0) then
         write(6, '(a5)') "works"
      else
         write(6, '(a5)') "fails"
      endif

      end
c

c-----------------------------------------------------------------------
c     wrapper program used to determine whether common blocks are 
c     initialized properly when the block data subprogram is not
c     explicitly called....sjs 11/3/08
c-----------------------------------------------------------------------
c
      program wrappr

      real*8 x

      common/tstcom/x

      if (x .eq. 1.234d0) then
         write(6, *) "works"
      else
         write(6, *) "fails"
      endif

      end
c

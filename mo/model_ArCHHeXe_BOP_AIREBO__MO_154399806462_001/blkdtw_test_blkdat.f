c-----------------------------------------------------------------------
c     block data subprogram used to test whether an explicit call is
c     needed or allowed....sjs 11/3/08
c-----------------------------------------------------------------------
c
      block data blkdat

      real*8 x

      common/tstcom/x

      data x/1.234d0/

      end
c

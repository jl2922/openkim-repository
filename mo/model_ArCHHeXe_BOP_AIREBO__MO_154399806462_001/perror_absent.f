c-----------------------------------------------------------------------
c      perror provides a dummy routine for the perror() intrinsic
c      routine, for machines that don't have it....sjs 3/14/05
c-----------------------------------------------------------------------
c
      subroutine perror(messg)

      character messg*160

      write(0, *) messg

      end

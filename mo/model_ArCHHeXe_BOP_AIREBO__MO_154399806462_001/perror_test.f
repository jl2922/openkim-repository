c-----------------------------------------------------------------------
c     wrapper program used to determine how to call perror().
c     ...sjs 7/4/07
c-----------------------------------------------------------------------
c
      program wrappr

      character messg*80

      messg="test"
      call perror(messg)

      end
c

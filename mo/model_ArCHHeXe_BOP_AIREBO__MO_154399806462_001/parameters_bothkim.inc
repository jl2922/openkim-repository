!     Change these parameters to alter the size of system that can be simulated
!     by the (statically allocated) potential routines. 
!     Suggested values work well together.
!     npmax  = maximum # of atoms (1000 or 7000 or 12500)
!     nlmax  = maximum # pairs in REBO neighbor list (10000 or 75000 or 150000)

      integer nlmax, npmax, nvmax
      parameter (npmax = 7000)
      parameter (nlmax = 75000)
      parameter (nvmax = 100000)

!     mxbcpa = maximum number of bond charges per atom. Allows for
!     4 covalent bonds per atom, plus a couple of extras for

!     conjugation
      integer mxbcpa
      parameter (mxbcpa = 6)

!     nbcmax = maximum number of bond charges. Should be no more than 
!     4 covalent bonds per atom, plus some nonbonded.

      integer nbcmax
      parameter (nbcmax = npmax * mxbcpa)

!     ndim = number of spatial dimensions

      integer ndim
      parameter (ndim   = 3)

!     potentials

      integer irebo
      parameter (irebo  = 1)

!     ntypes = maximum # of different LJ atom types

      integer ntypes
      parameter (ntypes = 9)

!     atom types

      integer iargon, icarb, ifluor, igerm, ihel, ihyd, ioxy, isili,    &
     &     ixenon
      parameter (ihyd   = 1)
      parameter (ihel   = 2)
      parameter (isili  = 3)
      parameter (igerm  = 4)
      parameter (iargon = 5)
      parameter (icarb  = 6)
      parameter (ixenon = 7)
      parameter (ioxy   = 8)
      parameter (ifluor = 9)

!     angular potential types

      integer ibop, isplin 
      parameter (isplin = 1)
      parameter (ibop   = 2)

!     bond charge types

      integer ipi, isigma
      parameter (isigma = 1)
      parameter (ipi    = 2)

!     maxc1p = maximum number of parameters in the singleton custom
!     parameter vector

      integer maxc1p
      parameter (maxc1p = 10)

!     maxc1t = maximum number of singleton custom parameters.
!     in principle it could by ntypes, but in practice it would only be
!     5 for CHNOF

      integer maxc1t
      parameter (maxc1t = 5)

!     maxc2p = maximum number of parameters in the pairwise custom parameter
!     vector

      integer maxc2p
      parameter (maxc2p = 13)

!     maxc2t = maximum number of pairwise custom parameters.
!     in principle it could be ntypes * (ntypes - 1) / 2, but in practice
!     it would only reach 20 for CHNOF

      integer maxc2t
      parameter (maxc2t = 20)

!     maxc3p = maximum number of parameters in the 3-body custom parameter
!     vector

      integer maxc3p
      parameter (maxc3p = 4)

!     maxc3t = maximum number of 3-body custom parameters.
!     in principle it could be ntypes ** 3, but in practice it would only
!     reach 125 for CHNOF

      integer maxc3t
      parameter (maxc3t = 125)


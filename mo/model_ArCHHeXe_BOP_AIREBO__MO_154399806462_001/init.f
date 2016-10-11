c-----------------------------------------------------------------------
c     ezinit is a wrapper for the init routine, that sets many of the 
c     initialization flags for the default cases of REBO and AIREBO
c     potentials....sjs 6/30/08
c-----------------------------------------------------------------------
c
      subroutine ezinit(np, r0, iatno, lpbc, cube, mypot, rbuffr, iverb)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np      = number of atoms
c     r0      = Cartesian coordinates of atoms
c     iatno   = atomic numbers
c     lpbc    = whether or not to use periodic boundaries
c     cube(d) = box length in dimension d (A)
c     mypot   = which potential to use (1 = REBO, 2 = AIREBO)
c     rbuffr  = buffer distance
c     iverb   = verbosity level (0 = silent ... 9 = debug)

      integer np
      real*8  r0(npmax, ndim)
      integer iatno(npmax)
      logical lpbc(ndim)
      real*8  cube(ndim)
      integer mypot
      real*8  rbuffr
      integer iverb

c     local variables

      integer ityc3p(maxc3t,3),
     .     ityc2p(maxc2t,2),
     .     ityc1p(maxc1t),
     .     igfunc, ipot, ncpair, nctrip, nctype
      logical usec3p(maxc3t,maxc3p),
     .     usec2p(maxc2t,maxc2p),
     .     usec1p(maxc1t,maxc1p),
     .     lbc, lckfil, les, lewald, ljdirw, 
     .     llj, lpdirw, lljdrw, lvor
      real*8  valc3p(maxc3t,maxc3p),
     .     valc2p(maxc2t,maxc2p),
     .     valc1p(maxc1t,maxc1p),
     .     ewkmxw, ewlkpl

      if (mypot .eq. 1) then
         write(istout, *) "Using REBO potential."
         write(istout, *) "Please cite as D.W. Brenner, O.A. ",
     .        "Shenderova, J.A. Harrison, S.J. Stuart, B. Ni, and S. ",
     .        "Sinnott, J. Phys.: Cond. Matt. 14, 783 (2002)."

c     set flags for REBO

         ipot = 1
         llj = .false.
         les = .false.
         lbc = .false.
         lewald = .false.
         ewlkpl = 6.d0
         ewkmxw = 1.6d0
         igfunc = isplin
         lvor = .false.
         ljdirw = .false.
         lpdirw = .false.
         lljdrw = .false.
         
      else if (mypot .eq. 2) then
         write(istout, *) "Using AIREBO potential."
         write(istout, *) "Please cite as S.J. Stuart, A.B. Tutein, ",
     .        "and J.A. Harrison, J. Chem. Phys. 112, 6472 (2000)."

c     set flags for AIREBO

         ipot = 1
         llj = .true.
         les = .false.
         lbc = .false.
         lewald = .false.
         ewlkpl = 6.d0
         ewkmxw = 1.6d0
         igfunc = isplin
         lvor = .false.
         ljdirw = .false.
         lpdirw = .false.
         lljdrw = .false.

      else
         write(isterr, *) "ezinit: Potential ", mypot, " unrecognized."
         write(isterr, *) "Please use 1 for REBO or 2 for AIREBO."
         call ioclos
         stop
      endif

      nctype = 0
      ncpair = 0
      nctrip = 0
      lckfil = .true.
      call preint(nctype, ityc1p, usec1p, valc1p,
     .     ncpair, ityc2p, usec2p, valc2p, 
     .     nctrip, ityc3p, usec3p, valc3p, lckfil)
      call runint(igfunc)
      call init(lbc, ipot, llj, rbuffr, les, lvor, ljdirw, lpdirw,
     .     lljdrw, iverb)
      call build(np, r0, iatno, cube, lpbc, lewald, ewlkpl, ewkmxw,
     .     lvor)

      return
      end
c
c-----------------------------------------------------------------------
c     preint initializes some stuff on the potential side, before the
c     run conditions are known or the system state is initialized. This
c     includes everything that MUST be initialized before the system
c     state is known, such as the mask identifying which parameters are
c     known. It can include anything that MAY be initialized before the
c     system state is known (like many potential parameters), although
c     those may also be initialized in runint or init....sjs 1/6/11
c-----------------------------------------------------------------------
c
      subroutine preint(nctype, ityc1p, usec1p, valc1p,
     .     ncpair, ityc2p, usec2p, valc2p, 
     .     nctrip, ityc3p, usec3p, valc3p, lckfil)

      include 'common_files.inc'
      include 'common_pots.inc'

c     nctype      = number of atom types with custom values
c     ityc1p(i)   = ith singleton type
c     usec1p(i,k) = whether to use a custom val. of par. k for singleton type i
c     valc1p(i,k) = custom value of parameter k for singleton type i
c     ncpair      = number of pair types with custom values
c     ityc2p(i,j) = atom type for jth half of pair type i
c     usec2p(i,k) = whether to use a custom value of par. k for pair type i
c     valc2p(i,k) = custom value of parameter k for pair type i
c     nctrip      = number of triple types with custom parameter values
c     ityc3p(i,j) = atom type for jth component of triple type i
c     usec3p(i,k) = whether to use a custom value of par. k for triple type i
c     valc3p(i,k) = custom value of parameter k for triple type i
c     lckfil      = whether to check files for custom parameter values

      integer nctype
      integer ityc1p(maxc1t)
      logical usec1p(maxc1t,maxc1p)
      real*8 valc1p(maxc1t,maxc1p)
      integer ncpair
      integer ityc2p(maxc2t,2)
      logical usec2p(maxc2t,maxc2p)
      real*8  valc2p(maxc2t,maxc2p)
c     integer nctrip
      integer ityc3p(maxc3t,3)
      logical usec3p(maxc3t,maxc3p)
      real*8  valc3p(maxc3t,maxc3p)
      logical lckfil

c     initialize common block data

      call blkdtw

c     initialize the REBO covalent bonding parameters

      call setcov(ncpair, ityc2p, usec2p, valc2p, 
     .     nctrip, ityc3p, usec3p, valc3p, lckfil)

c     initialize some LJ parameters

      call ljpars(nctype, ityc1p, usec1p, valc1p, lckfil)

      end
c
c-----------------------------------------------------------------------
c     runint initializes some stuff on the potential side, after the run
c     conditions are known but before the system state is initialized.
c-----------------------------------------------------------------------
c
      subroutine runint(igfunc)

      include 'common_files.inc'
      include 'common_pots.inc'

c     igfunc = style of angular potential to use

      integer igfunc

c     initialize angular potential stuff

      call setang(igfunc)

      end
c
c-----------------------------------------------------------------------
c     init initializes some stuff on the potential side, after the
c     system state is initialized. This includes everything that MUST be
c     initialized after the system state is known, such as atom counts.
c     It can include anything that MAY be initialized after the system
c     state is known (like many potential parameters), although those
c     may also be initialized in preint....sjs 1/6/11
c-----------------------------------------------------------------------
c
      subroutine init(lbc, ipot, llj, rbuffr, les, lvor, ljdirw, lpdirw,
     .     lljdrw, iverb)

      include 'common_files.inc'
      include 'common_pots.inc'

c     lbc         = whether we are using bond charges
c     ipot        = which potential to use
c     llj         = whether or not to use Lennard-Jones
c     rbuffr      = buffer distance
c     les         = whether or not to use elecrostatics
c     lvor        = whether or not to calculate Voronoi volumes
c     ljdirw      = whether to calculate J(r) directly (vs lookup table)
c     lpdirw      = whether to calculate pairwise terms directly
c     lljdrw     = whether to calculate lj pair terms directly
c     iverb       = verbosity level (0 = silent ... 9 = debug)

      logical lbc
      integer ipot
      logical llj
      real*8  rbuffr
      logical les
      logical lvor
      logical ljdirw
      logical lpdirw
      logical lljdrw
      integer iverb

c     initialize some REBO potential parameters

      call setin(rbuffr)

c     initialize some electrostatics parameters

      if (les) then
         call setes(lbc, ljdirw)
      endif

      if (.not. lvor) then

c        initialize some potential parameters

         call setpp(rbuffr, ipot, llj, 
     .       lpdirw, lljdrw, iverb)

      endif

      end

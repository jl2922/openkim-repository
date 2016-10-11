c-----------------------------------------------------------------------
c     setin sets some REBO potential parameters....sjs
c-----------------------------------------------------------------------
c     MD-specific initialization have been moved to initmd.
c     ...sjs 4/29/99
c-----------------------------------------------------------------------
c
      subroutine setin(rbuffr)

      include 'common_files.inc'
      include 'common_pots.inc'

c     rbuffr = buffer distance for pair list (A)

      real*8 rbuffr

c     local variables

      character name*20, unit*20
      integer itype, jtype
      logical symmet
      real*8  deldij

c     when new atom types are added here, also be sure to:
c     - fix kt and kt2 in blkdat
c     - adjust LJ parameters in ljparam
c     - specify mass in initmd
c     - increase ntypes in parameters_both.inc

      do 90 i = 1, ntypes
         noa(i) = 0
 90   continue

      pi=acos(-1.d0)

c     pre-zero some REBO parameters

      do 110 i = 1, ntypes
         do 100 j = 1, ntypes
            repmin(i,j) = 0.d0
            repmax(i,j) = 0.d0
            dijmin(i,j) = 0.d0
            pibyd(i,j)  = 1.d0
 100     continue
 110  continue

c     switching function that blends ZBL and REBO repulsion

c     for C-C, H-H, and C-H, these have been roughly fit by TE,
c     although there are some problems (nonmonotiic for C-C)

c     for Si-Si, Si-C, and Si-H, the switching acts over the range
c     where the ZBL potential goes from 70 eV to 30 eV

c     inside edge of switching function

      repmin(icarb,icarb) = 0.3d0
      repmin(ihyd,ihyd)   = 0.43d0
      repmin(icarb,ihyd)  = 0.80d0
      repmin(ihyd,icarb)  = repmin(icarb,ihyd)
      repmin(isili,isili) = 0.92d0
      repmin(isili,icarb) = 0.76d0
      repmin(icarb,isili) = repmin(isili,icarb)
      repmin(isili,ihyd) = 0.44d0
      repmin(ihyd,isili) = repmin(isili,ihyd)

c     outside edge of switching function

      repmax(icarb,icarb) = 1.17d0
      repmax(ihyd,ihyd)   = 0.67d0
      repmax(icarb,ihyd)  = 0.99d0
      repmax(ihyd,icarb)  = repmax(icarb,ihyd)
      repmax(isili,isili) = 0.8d0
      repmax(isili,icarb) = 0.8d0
      repmax(icarb,isili) = repmax(isili,icarb)
      repmax(isili,ihyd) = 0.8d0
      repmax(ihyd,isili) = repmax(isili,ihyd)

c     (ZBL is currently turned off for C-C because of problems with
c     switching range)

      repmin(icarb,icarb) = 0.d0
      repmax(icarb,icarb) = 0.d0
      
c     (ZBL is currently turned off for all interactions, because it
c     is not yet fitted well)

      repmin(ihyd,ihyd)   = 0.d0
      repmax(ihyd,ihyd)   = 0.d0
      repmin(icarb,ihyd)  = 0.d0
      repmin(ihyd,icarb)  = repmin(icarb,ihyd)
      repmax(icarb,ihyd)  = 0.d0
      repmax(ihyd,icarb)  = repmax(icarb,ihyd)

c     The values above can be overridden by custom parameter files

      symmet = .true.

      name = 'rZBLmin'
      unit = 'A'
      call cuspa2(repmin, name, unit, symmet)

      name = 'rZBLmax'
      unit = 'A'
      call cuspa2(repmax, name, unit, symmet)

c     dijmin = inside edge of switching function that turns off REBO potential
c     dijmax = outside edge of switching function that turns off REBO potential

c     there is no good default, but try anyway

      deldij = 0.25d0
      do 130 itype = 1, ntypes
         if (lrebot(itype)) then
            dijmax(itype,itype) = 2.00d0
            dijmin(itype,itype) = dijmax(itype,itype) - deldij
            do 120 jtype = itype+1, ntypes
               if (lrebot(jtype)) then
                  dijmax(itype,jtype) = 2.00d0
                  dijmax(jtype,itype) = dijmax(itype,jtype)
                  dijmin(itype,jtype) = dijmax(itype,jtype) - deldij
                  dijmin(jtype,itype) = dijmin(itype,jtype)
               endif
 120        continue
         endif
 130  continue

      dijmin(icarb,icarb) = 1.7d0
      dijmax(icarb,icarb) = 2.0d0

      dijmin(ihyd,ihyd)   = 1.1d0
      dijmax(ihyd,ihyd)   = 1.7d0

      dijmin(icarb,ihyd)  = 1.3d0
      dijmax(icarb,ihyd)  = 1.8d0

      dijmin(ihyd,icarb)  = dijmin(icarb,ihyd)
      dijmax(ihyd,icarb)  = dijmax(icarb,ihyd)

c     Liegi's values for F

      dijmin(ifluor,ifluor) = 1.58d0
      dijmax(ifluor,ifluor) = 1.65d0

c     Liegi or Brad says:
c        dijmax(icarb,ifluor) = 2.25 leads to no -F- bond
c        dijmax(icarb,ifluor) = 2.55 leads to serious -F- bonds
c        dijmax(icarb,ifluor) = 2.40 leads to serious -F- bonds
c        dijmax(icarb,ifluor) = 2.30 leads to -F- bonds

      dijmin(icarb,ifluor) = 1.55d0
      dijmax(icarb,ifluor) = 2.05d0

      dijmin(ifluor,icarb) = dijmin(icarb,ifluor)
      dijmax(ifluor,icarb) = dijmax(icarb,ifluor)

c     HF: covalent bond is ~0.92 A
c     geminal pairs at sp3 carbon are at ~2.00 A

      dijmin(ihyd,ifluor) = 1.30d0
      dijmax(ihyd,ifluor) = 1.90d0

      dijmin(ifluor,ihyd) = dijmin(ihyd,ifluor)
      dijmax(ifluor,ihyd) = dijmax(ihyd,ifluor)

c     OH: covalent bond is ~0.97 A; 2nd neighbor in H2O2 is ~1.83 A

      dijmax(ihyd,ioxy) = 1.75d0
      dijmax(ioxy,ihyd) = dijmax(ihyd,ioxy)

      dijmin(ihyd,ioxy) = dijmax(ihyd,ioxy) - deldij
      dijmin(ioxy,ihyd) = dijmin(ihyd,ioxy)

c     CO: covalent single bond is ~1.42 A
c     2nd neighbor in acetic acid is ~2.37 A
      
      dijmax(icarb,ioxy) = 2.00d0
      dijmax(ioxy,icarb) = dijmax(icarb,ioxy)
      
      dijmin(icarb,ioxy) = dijmax(icarb,ioxy) - deldij
      dijmin(ioxy,icarb) = dijmin(icarb,ioxy)

c     OF: covalent single bond is ~1.42 A typically, ~1.58 A in O2F2
c     2nd neighbor in O2F2 is ~2.29 A

      dijmax(ioxy,ifluor) = 2.00d0
      dijmax(ifluor,ioxy) = 2.00d0
      
      dijmin(ioxy,ifluor) = dijmax(ioxy,ifluor) - deldij
      dijmin(ifluor,ioxy) = dijmin(ioxy,ifluor)

c     OO: covalent single bond is ~1.48 A
c     2nd neighbor is ~2.17 A in ozone

      dijmax(ioxy,ioxy) = 2.d0
      dijmin(ioxy,ioxy) = dijmax(ioxy,ioxy) - deldij

c     Si-Si: covalent single bond is ~2.35 A
c     2nd neighbor is ~3.84 in Si(diamond)

      dijmax(isili,isili) = 3.25d0
      dijmin(isili,isili) = dijmax(isili,isili) - deldij

c     Si-H: covalent bond is 1.48 A in silane
c     2nd neighbor could be as short as 2.18 A for a tetrahedral Si-O-H

      dijmax(isili,ihyd) = 2.d0
      dijmax(ihyd,isili) = dijmax(isili,ihyd)

      dijmin(isili,ihyd) = dijmax(isili,ihyd) - deldij
      dijmin(ihyd,isili) = dijmin(isili,ihyd)

c     Si-C: covalent bond can be as short as 1.72 A for a Si=C double bond, 
c     1.84 to 1.93 A for a Si-C single bond.
c     2nd neighbor could be as short as 2.38 A for a tetrahedral C-O-Si

      dijmax(isili,icarb) = 2.25d0
      dijmax(icarb,isili) = dijmax(isili,icarb)
      
      dijmin(isili,icarb) = dijmax(isili,icarb) - deldij
      dijmin(icarb,isili) = dijmin(isili,icarb)

c     Si-O: covalent bond is 1.50 A (radical) up to ~1.66 A in silanols
c     2nd neighbor should be about 2.43 A in tetrahedral Si-O-O

      dijmax(isili,ioxy) = 2.25d0
      dijmax(ioxy,isili) = dijmax(isili,ioxy)

      dijmin(isili,ioxy) = dijmax(isili,ioxy) - deldij
      dijmin(ioxy,isili) = dijmin(isili,ioxy)

c     Si-F: covalent bond is 1.56 A in SiF4
c     2nd neighbor in should be about 3.22 A in Si2F6

      dijmax(isili,ifluor) = 2.5d0
      dijmax(ifluor,isili) = dijmax(isili,ifluor)
      
      dijmin(isili,ifluor) = dijmax(isili,ifluor) - deldij
      dijmin(ifluor,isili) = dijmin(isili,ifluor)

      dijmin(igerm,igerm) = 2.7d0
      dijmin(isili,igerm) = sqrt(dijmin(isili,isili) * 
     .     dijmin(igerm,igerm))
      dijmin(igerm,isili) = dijmin(isili,igerm)

      dijmax(igerm,igerm) = 3.0d0
      dijmax(isili,igerm) = sqrt(dijmax(isili,isili) *
     .     dijmax(igerm,igerm))
      dijmax(igerm,isili) = dijmax(isili,igerm)

c     pi / width of REBO switching function

      do 150 itype = 1, ntypes
         if (lrebot(itype)) then
            pibyd(itype,itype) = pi /
     .           (dijmax(itype,itype) - dijmin(itype,itype))
            do 140 jtype = 1, ntypes
               if (lrebot(jtype)) then
                  pibyd(itype,jtype) = pi /
     .                 (dijmax(itype,jtype) - dijmin(itype,jtype))
                  pibyd(jtype,itype) = pibyd(itype,jtype)
               endif
 140        continue
         endif
 150  continue

c     keep the squares of the outer REBO cutoffs.
c     and precalculate the square of the outer cutoff for constructing
c     the pair list

      do 210 i = 1, ntypes
         do 200 j = 1, ntypes
            dijmx2(i,j) = dijmax(i,j) ** 2
            r2rbmx(i,j) = (dijmax(i,j) + rbuffr) ** 2
 200     continue
 210  continue

      return 
      end 
c     
c-----------------------------------------------------------------------
c     setang sets some variables relating to the angular potential
c     ...sjs 2/9/11
c-----------------------------------------------------------------------
c
      subroutine setang(igfunc)
      
      include 'common_files.inc'
      include 'common_pots.inc'

c     igfunc = style of angular potential to use

      integer igfunc

c     local variables

      integer itype, jtype, ktype

c     angular parameters are set in blkdtw, but this is where we learn
c     which parameters are known

      do 120 itype = 1, ntypes
         do 110 jtype = 1, ntypes
            do 100 ktype = 1, ntypes
               lhvang(itype,jtype,ktype) = .false.
 100        continue
 110     continue
 120  continue

      if (igfunc .eq. isplin) then

         do 210 itype = 1, ntypes
            do 200 jtype = 1, ntypes
               lhvang(itype,icarb,jtype) = .true.
               lhvang(itype,ihyd,jtype) = .true.
               lhvang(itype,ifluor,jtype) = .true.
 200        continue
 210     continue

      else

      endif

      return
      end
c
c-----------------------------------------------------------------------
c     setpp calls the parameter setting routines for the various 
c     potentials.
c-----------------------------------------------------------------------
c
      subroutine setpp(rbuffr, ipot, llj, lpdirw, lljdrw, iverb)

      include 'common_files.inc'
      include 'common_pots.inc'

c     rbuffr = buffer distance for pairlist
c     ipot = potential (1 = REBO)
c     llj = whether LJ interactions are being used
c     lpdirw = whether pairwise terms are evaluated directly
c     lljdrw = whether lj pair terms are evaluated directly
c     iverb = verbosity level

      real*8 rbuffr
      integer*4 ipot
      logical llj
      logical lpdirw
      logical lljdrw
      integer iverb
      
      lljdir = lljdrw
      lpdir = lpdirw

      if (ipot .eq. irebo) then
         call param(rbuffr, iverb)
         if (.not. lpdir) then
            call mtable
         endif
      endif

      if (llj) then
         call ljcuts(rbuffr)
      endif
      return
      end 
c     
c-----------------------------------------------------------------------
c     setes sets some variables that are needed for electrostatics.
c     ...sjs 7/25/07
c-----------------------------------------------------------------------
c
      subroutine setes(lbc, ljdirw)

      include 'common_files.inc'
      include 'common_pots.inc'

c     lbc     = whether we are using bond charges
c     ljdirw  = whether to do direct J(r) calculations

      logical lbc
      logical ljdirw

c     local variables

      character name*20, unit*20
      integer iatom, idim
      logical symmet
      real*8  rmax, rmin, rprec, rtest

c     cos(0) = 1 and sin(0) = 0.

      do 110 idim = 1, ndim
         do 100 iatom = 1, npmax
            coskx(iatom,0,idim) = 1.d0
            sinkx(iatom,0,idim) = 0.d0
 100     continue
 110  continue

c     initialize parameters (default values appropriate for
c     non-fluctuating charges)

      do 200 itype = 1, ntypes
         chi0(itype) = 0.d0
         zeta(itype) = 0.d0
         iseety(itype) = .false.
 200  continue

c     chi0 is an electronegativity in units of eV/e
c     zeta is an orbital exponent in the Slater orbitals used to calculate
c     Coulomb interactions and atomic hardnesses, in units of 1/A

c     different parameters are used for the atomic charge and bond charge
c     models

      if (lbc) then

c     single-atom parameters:

c     hydrogen

c     chi(H) comes from Rappe & Goddard JPC 95 3358 '91. The zero of
c     electronegativity is arbitrary

         chi0(ihyd) = 4.5280d0
         zeta(ihyd) = 1.34476d0
         iseety(ihyd) = .true.

c     carbon

c     I think chi(C) - chi(H) was fit to reproduce the octupole moment
c     of CH4. (sjs) Indeed it is. (bjw)
c     zeta(C) was fit in conjunction with j_CC^sigma and j_CC^pi
c     to reproduce longitudinal bond polarizabilities for single, double,
c     and triple CC bonds out of Adv in Phys Org Chem Vol. 3 Chapter 1 (bjw)

         chi0(icarb) = 4.8769d0
         zeta(icarb) = 1.76115d0
         iseety(icarb) = .true.

c     oxygen
c     draft parameters taken from Rappe & Goddard

         chi0(ioxy) = 8.741d0
         zeta(ioxy) = 0.9745d0 / bohr
         iseety(ioxy) = .true.

c     fluorine

         chi0(ifluor) = 6.1703d0
         zeta(ifluor) = 2.11411d0
         iseety(ifluor) = .true.

c     silicon
c     draft parameters taken from Rappe & Goddard

         chi0(isili) = 4.168d0
         zeta(isili) = 0.7737d0 / bohr
         iseety(isili) = .true.

c     pairwise (bond) parameters:

c     default bond hardness (reverts to atomic charge model)

         do 310 itype = 1, ntypes
            do 300 jtype = 1, ntypes
               jbh(itype,jtype,isigma) = 0.d0
               jbh(itype,jtype,ipi) = 0.d0
 300        continue
 310     continue

c     The unit conversions are to convert bond hardnesses in 1/A to V/e

         jbh(ihyd,ihyd,isigma) = 0.39166d0 * fcm2A * fe2esu
     .        * fe2esu * ferg2J * fC2e
         jbh(ihyd,ifluor,isigma) = 0.66589d0 * fcm2A * fe2esu
     .        * fe2esu * ferg2J * fC2e
         jbh(ifluor,ihyd,isigma) = jbh(ihyd,ifluor,isigma)
c     j_CH was 0.26766 for methane fit
c     j_CH was 0.07741 for dodecane fit
c     j_CH was 0.12949 for hexane fit
c     Current j_CH fit is octane
         jbh(icarb,ihyd,isigma) = 0.10811d0 * fcm2A * fe2esu
     .        * fe2esu * ferg2J * fC2e
         jbh(ihyd,icarb,isigma) = jbh(icarb,ihyd,isigma)
         jbh(icarb,icarb,isigma) = 2.14009d0 * fcm2A * fe2esu
     .        * fe2esu * ferg2J * fC2e
         jbh(icarb,icarb,ipi) = 0.53052d0 * fcm2A * fe2esu
     .        * fe2esu * ferg2J * fC2e
         jbh(icarb,ifluor,isigma) = 0.48279d0 * fcm2A * fe2esu
     .        * fe2esu * ferg2J * fC2e
         jbh(ifluor,icarb,isigma) = jbh(icarb,ifluor,isigma)
         jbh(ifluor,ifluor,isigma) = 0.70606d0 * fcm2A * fe2esu
     .        * fe2esu * ferg2J * fC2e

      else

c     hydrogen

         chi0(ihyd) = 4.12415d0
         zeta(ihyd) = 1.09864 * fA2au
         iseety(ihyd) = .true.

c     carbon

c     to first approximation, CH bonds are nonpolar
c     (i.e. the CH qAIREBO model has not been fit yet)
c     zeta(C) fit to reproduce bond (total) polarizability of 0.99 A^3 for C-C
c     at r=1.54 (Adv Phys Org Chem 3 1 '65)

         chi0(icarb) = 11.5898d0
         zeta(icarb) = 2.6964d0 * fA2au
         iseety(icarb) = .true.

c     oxygen
c     draft parameters taken from Rappe & Goddard

         chi0(ioxy) = 8.741d0
         zeta(ioxy) = 0.9745d0 / bohr
         iseety(ioxy) = .true.

c     fluorine

c     chi(F) - chi(C) fit by Brennan to give q_F = -0.14757 in CF4, obtained
c     from octupole moment of 11.0e-50 C m^3 (Mol Phys 32 161 '76)
c     zeta(F) fit to reproduce longitudinal electronic polarizability of F2
c     of 12.4443 au (CPL 442 265 '07)

         chi0(ifluor) = 20.1878d0
         zeta(ifluor) = 1.78032d0 * fA2au
         iseety(ifluor) = .true.

c     silicon
c     draft parameters taken from Rappe & Goddard

         chi0(isili) = 4.168d0
         zeta(isili) = 0.7737d0 / bohr
         iseety(isili) = .true.

      endif

c     The values above can be overridden by custom parameter files

      name = 'chi0'
      unit = 'V'
      call cuspa1(chi0, name, unit)

      name = 'zeta'
      unit = '1/Angstrom'
      call cuspa1(zeta, name, unit)

      name = 'jsigma'
      unit = '1/Angstrom'
      symmet = .true.
      call cuspa2(jbh(1,1,isigma), name, unit, symmet)

      name = 'jpi'
      unit = '1/Angstrom'
      call cuspa2(jbh(1,1,ipi), name, unit, symmet)

c     valence shells, for use in J(r) calculations

      do 400 itype = 1, ntypes
         ivalsh(itype) = 0
 400  continue

c     hydrogen has a 1s valence shell

      ivalsh(ihyd) = 1

c     carbon has a 2s valence shell

      ivalsh(icarb) = 2

c     oxygen has a 2s valence shell

      ivalsh(ioxy) = 2

c     fluorine has a 2s valence shell

      ivalsh(ifluor) = 2

c     silicon has a 3s valence shell

      ivalsh(isili) = 3

c     determine the upper bound for the J(r) lookup tables,
c     as a little longer than the distance needed to obtain an
c     error of eserrt eV/e^2 or less (round to nearest rprec).

      rprec = 0.1d0
      do 510 itype = 1, ntypes
         rjmax(itype,itype) = 0.d0
         do 500 jtype = itype, ntypes
            if (.not. iseety(itype) .or. .not. iseety(jtype)) then
               rjmax(itype,jtype) = 0.d0
            else
               rmin = 0.d0
               rtest = 1.d0
 480           err = dabs(epsinv / rtest - quadj(itype,jtype,rtest))
               if (err .ge. eserrt) then
                  rmin = rtest
                  rtest = rtest * 2.d0
                  go to 480
               endif
               rmax = rtest
 499           rtest = 0.5d0 * (rmax + rmin)
               err = dabs(epsinv / rtest - quadj(itype,jtype,rtest))
               if (err .ge. eserrt) then
                  rmin = rtest
               else
                  rmax = rtest
               endif
               if (rmax - rmin .gt. rprec) then
                  go to 499
               endif
               rjmax(itype,jtype) = (int(rmax / rprec) + 1) * rprec
            endif
            rjmax(jtype,itype) = rjmax(itype,jtype)
 500     continue
 510  continue

c     fill the J(r) lookup tables, unless we are doing direct J(r)
c     calculations

      ljdir = ljdirw
      if (.not. ljdirw) then
         call Jinit
      endif

c     set the J0 hardness values from J(0).
c     make sure to use quadj() in order to get a real J(r) calculation,
c     rather than J(), which will switch off the interaction at short 
c     distances.

      J0min = bigpos
      do 600 itype = 1, ntypes
         if (iseety(itype)) then
c$$$            J0(itype) = J(itype,itype,0.d0)
            J0(itype) = quadj(itype,itype,0.d0)
            if (J0(itype) .lt. J0min .and. J0(itype) .gt. 0.d0) then
               J0min = J0(itype)
            endif
         else
            J0(itype) = 0.d0
         endif
 600  continue

      end
c
c-----------------------------------------------------------------------
c     stvolp sets everything on the potential side that depends on the
c     box size....sjs 7/26/07
c-----------------------------------------------------------------------
c
      subroutine stvolp(cube, vol, lewald, ewlkpl, ewkmxw)

      include 'common_files.inc'
      include 'common_pots.inc'

c     cube   = box size
c     vol    = box volume
c     lewald = whether to perform Ewald sums
c     ewlkpl = Ewald kappa * (min) box length
c     ewkmxw = Ewald max |k|

      real*8  cube(ndim)
      real*8  vol
      logical lewald
      real*8  ewlkpl
      real*8  ewkmxw

c     local variables

      integer idim
      real*8 boxmin
      
      if (lewald) then

c     set some Ewald variables

         r2pb3v = 2.d0 * pi / 3.d0 / vol * epsinv
         r4pb3v = 2.d0 * r2pb3v

         boxmin = bigpos
         do 100 idim = 1, ndim
            if (cube(idim) .lt. boxmin) then
               boxmin = cube(idim)
            endif
 100     continue
         ewlkap = ewlkpl / boxmin
         ewlkp2 = ewlkap * ewlkap
         rkbrp = ewlkap / sqrt(pi) * epsinv
         r2kbrp = 2.d0 * rkbrp
         ewlkmx = ewkmxw

c     update the lookup table of k vectors

         call setkvc(cube, vol)
      endif

      return
      end
c
c-----------------------------------------------------------------------
c     setcov initializes the values of the REBO covalent bonding
c     parameters, overriding the standard values with custom values
c     as needed....sjs 6/16/10
c-----------------------------------------------------------------------
c
      subroutine setcov(ncpair, ityc2p, usec2p, valc2p, 
     .     nctrip, ityc3p, usec3p, valc3p, lckfil)

      include 'common_files.inc'
      include 'common_pots.inc'

c     ncpair      = number of pair types with custom values
c     ityc2p(i,j) = atom type for jth half of pair type i
c     usec2p(i,k) = whether to use a custom value of par. k for pair type i
c     valc2p(i,k) = custom value of parameter k for pair type i
c     nctrip      = number of triple types with custom parameter values
c     ityc3p(i,j) = atom type for jth type of triple type i
c     usec3p(i,k) = whether to use a custom value of par. k for triple type i
c     valc3p(i,k) = custom value of parameter k for triple type i
c     lckfil      = whether to check files for custom parameter values

      integer ncpair
      integer ityc2p(maxc2t,2)
      logical usec2p(maxc2t,maxc2p)
      real*8  valc2p(maxc2t,maxc2p)
      integer nctrip
      integer ityc3p(maxc3t,3)
      logical usec3p(maxc3t,maxc3p)
      real*8  valc3p(maxc3t,maxc3p)
      logical lckfil

c     local variables

      character name*20, unit*20
      integer ipair, iparam, itype, jtype, ktype
      logical symmet
      real*8  betaRE, DeREBO, ReREBO, SREBO

c     make sure all parameters are initialized to some value

      do 110 itype = 1, ntypes
         do 100 jtype = 1, ntypes
            b1(itype,jtype) = 0.d0
            beta1(itype,jtype) = 0.d0
            b2(itype,jtype) = 0.d0
            beta2(itype,jtype) = 0.d0
            b3(itype,jtype) = 0.d0
            beta3(itype,jtype) = 0.d0
            capa(itype,jtype) = 0.d0
            alfa(itype,jtype) = 0.d0
            capq(itype,jtype) = 0.d0
c           flag them as missing
            lhvcov(itype,jtype) = .false.
 100     continue
 110  continue

c     Carbon-carbon

      b1(icarb,icarb)    = 12388.79197798375d0
      beta1(icarb,icarb) =     4.720452312717397d0
      b2(icarb,icarb)    =    17.56740646508968d0
      beta2(icarb,icarb) =     1.433213249951261d0
      b3(icarb,icarb)    =    30.71493208065162d0
      beta3(icarb,icarb) =     1.382691250599169d0
      capa(icarb,icarb)  = 10953.54416216992d0
      alfa(icarb,icarb)  =     4.746539060659529d0
      capq(icarb,icarb)  =     0.3134602960832605d0
      lhvcov(icarb,icarb) = .true.

c     hydrogen-hydrogen

c     New values (6/17/00), parametrized to give for H2:
c     re = 0.7461 A
c     E(re) = -4.5059 eV  --> BDE = 102.4 kcal/mol at 298 K
c     k = (d2E/dr2)|re = 34 eV/A^2

      b1(ihyd,ihyd)    = 28.2297d0
      beta1(ihyd,ihyd) =  1.708d0
      b2(ihyd,ihyd)    =  0.d0
      beta2(ihyd,ihyd) =  1.d0
      b3(ihyd,ihyd)    =  0.d0
      beta3(ihyd,ihyd) =  1.d0
      capa(ihyd,ihyd)  = 31.6731d0
      alfa(ihyd,ihyd)  =  3.536d0
      capq(ihyd,ihyd)  =  0.370d0
      lhvcov(ihyd,ihyd) = .true.

c     Original REBO values:

c      b1(ihyd,ihyd)    = 29.6325931d0
c      beta1(ihyd,ihyd) =  1.715892169856421d0
c      b2(ihyd,ihyd)    =  0.0d0
c      beta2(ihyd,ihyd) =  1.0d0
c      b3(ihyd,ihyd)    =  0.0d0
c      beta3(ihyd,ihyd) =  1.0d0
c      capa(ihyd,ihyd)  = 32.81735574722296d0
c      alfa(ihyd,ihyd)  =  3.536298648376465d0
c      capq(ihyd,ihyd)  =  0.3704714870452888d0

c     fluorine-fluorine
c     Liegi Hu's parameters

      b1(ifluor,ifluor)    = 1886.761284d0
      beta1(ifluor,ifluor) =  3.6d0
      b2(ifluor,ifluor)    =  0.d0
      beta2(ifluor,ifluor) =  1.d0
      b3(ifluor,ifluor)    =  0.d0
      beta3(ifluor,ifluor) =  1.d0
      capa(ifluor,ifluor)  = 2031.761284d0
      alfa(ifluor,ifluor)  =  4.d0
      capq(ifluor,ifluor)  =  0.5630122773d0
      lhvcov(ifluor,ifluor) = .true.

c     carbon-hydrogen

c     Original REBO values:

      b1(ihyd,icarb)    = 32.35518665873256d0
      beta1(ihyd,icarb) = 1.434458059249837d0
      b2(ihyd,icarb)    = 0.0d0
      beta2(ihyd,icarb) = 1.0d0
      b3(ihyd,icarb)    = 0.0d0
      beta3(ihyd,icarb) = 1.0d0
      capa(ihyd,icarb)  = 149.9409872288120d0
      alfa(ihyd,icarb)  = 4.102549828548784d0
      capq(ihyd,icarb)  = 0.3407757282257080d0
      lhvcov(ihyd,icarb) = .true.

      b1(icarb,ihyd)    = b1(ihyd,icarb)
      beta1(icarb,ihyd) = beta1(ihyd,icarb)
      b2(icarb,ihyd)    = b2(ihyd,icarb)
      beta2(icarb,ihyd) = beta2(ihyd,icarb)
      b3(icarb,ihyd)    = b3(ihyd,icarb)
      beta3(icarb,ihyd) = beta3(ihyd,icarb)
      capa(icarb,ihyd)  = capa(ihyd,icarb)
      alfa(icarb,ihyd)  = alfa(ihyd,icarb)
      capq(icarb,ihyd)  = capq(ihyd,icarb)
      lhvcov(icarb,ihyd) = .true.

c     Carbon-Fluorine
c     Liegi Hu's parameters

      b1(ifluor,icarb)    = 77.74643079d0
      beta1(ifluor,icarb) = 1.7d0
      b2(ifluor,icarb)    = 0.0d0
      beta2(ifluor,icarb) = 1.0d0
      b3(ifluor,icarb)    = 0.0d0
      beta3(ifluor,icarb) = 1.0d0
      capa(ifluor,icarb)  = 408.4269065d0
      alfa(ifluor,icarb)  = 4.d0
      capq(ifluor,icarb)  = 0.7953042889d0
      lhvcov(ifluor,icarb) = .true.

      b1(icarb,ifluor)    = b1(ifluor,icarb)
      beta1(icarb,ifluor) = beta1(ifluor,icarb)
      b2(icarb,ifluor)    = b2(ifluor,icarb)
      beta2(icarb,ifluor) = beta2(ifluor,icarb)
      b3(icarb,ifluor)    = b3(ifluor,icarb)
      beta3(icarb,ifluor) = beta3(ifluor,icarb)
      capa(icarb,ifluor)  = capa(ifluor,icarb)
      alfa(icarb,ifluor)  = alfa(ifluor,icarb)
      capq(icarb,ifluor)  = capq(ifluor,icarb)
      lhvcov(icarb,ifluor) = .true.

c     Si-Si
c     rough draft of parameters
c     Taken from JD Schall et al, PRB 77, 115209 (2008)

      b1(isili,isili)     = 92.74551d0
      beta1(isili,isili)  = 1.72687d0
      b2(isili,isili)     = 255.329d0
      beta2(isili,isili)  = 1.64617d0
      b3(isili,isili)     = -3.4026d0
      beta3(isili,isili)  = 132.454d0
      capa(isili,isili)   = 90.1964d0
      alfa(isili,isili)   = 2.13083d0
      capq(isili,isili)   = 15.6614d0
      lhvcov(isili,isili) = .true.

c     Si-H
c     rough draft of parameters
c     taken from Dyson & Smith, Surf Sci 355, 140 (1996), shoehorned
c     into AIREBO form

      DeREBO = 3.140d0
      SREBO = 1.8177
      betaRE = 1.6897d0
      ReREBO = 1.441d0

      beta1(isili,ihyd) = sqrt(2. / SREBO) * betaRE
      b1(isili,ihyd) = DeREBO * SREBO / (SREBO - 1.d0) *
     .     exp(beta1(isili,ihyd) * ReREBO)
      alfa(isili,ihyd) = sqrt(2. * SREBO) * betaRE
      capa(isili,ihyd) = DeREBO / (SREBO - 1.d0) *
     .     exp(alfa(isili,ihyd) * ReREBO)
      lhvcov(isili,ihyd) = .true.

      beta1(ihyd,isili) = beta1(isili,ihyd)
      b1(ihyd,isili)    = b1(isili,ihyd)
      alfa(ihyd,isili)  = alfa(isili,ihyd)
      capa(ihyd,isili)  = capa(isili,ihyd)
      lhvcov(ihyd,isili) = .true.

c     Si-C
c     rough draft of parameters
c     taken from Dyson & Smith, Surf Sci 355, 140 (1996), shoehorned
c     into AIREBO form

      DeREBO = 4.510d0
      SREBO = 1.492d0
      betaRE = 1.698d0
      ReREBO = 1.7631d0
      
      beta1(isili,icarb) = sqrt(2. / SREBO) * betaRE
      b1(isili,icarb)    = DeREBO * SREBO / (SREBO - 1.d0) *
     .     exp(beta1(isili,icarb) * ReREBO)
      alfa(isili,icarb)  = sqrt(2. * SREBO) * betaRE
      capa(isili,icarb)  = DeREBO / (SREBO - 1.d0) *
     .     exp(alfa(isili,icarb) * ReREBO)
      lhvcov(isili,icarb) = .true.

      beta1(icarb,isili) = beta1(isili,icarb)
      b1(icarb,isili)    = b1(isili,icarb)
      alfa(icarb,isili)  = alfa(isili,icarb)
      capa(icarb,isili)  = capa(isili,icarb)
      lhvcov(icarb,isili) = .true.

c     BOP-style angular potential parameters

      do 150 itype = 1, ntypes
         do 140 jtype = 1, ntypes
            do 130 ktype = 1, ntypes
               if (lrebot(jtype)) then
                  glamda(itype,jtype,ktype) = 1.d0
               else
                  glamda(itype,jtype,ktype) = 0.d0
               endif
 130        continue
 140     continue
 150  continue

      if (lckfil) then

c     The values above can be overridden by custom parameter files

         symmet = .true.

         name = 'B1'
         unit = 'eV'
         call cuspa2(b1, name, unit, symmet)

         name = 'beta1'
         unit = '1/Angstrom'
         call cuspa2(beta1, name, unit, symmet)

         name = 'B2'
         unit = 'eV'
         call cuspa2(b2, name, unit, symmet)

         name = 'beta2'
         unit = '1/Angstrom'
         call cuspa2(beta2, name, unit, symmet)

         name = 'B3'
         unit = 'eV'
         call cuspa2(b3, name, unit, symmet)

         name = 'beta3'
         unit = '1/Angstrom'
         call cuspa2(beta3, name, unit, symmet)

         name = 'A'
         unit = 'eV'
         call cuspa2(capa, name, unit, symmet)

         name = 'alpha'
         unit = '1/Angstrom'
         call cuspa2(alfa, name, unit, symmet)

         name = 'Q'
         unit = 'Angstrom'
         call cuspa2(capq, name, unit, symmet)

         symmet = .true.

         name = 'lambda'
         unit = ''
         call cuspa3(glamda, name, unit, symmet)

      endif

c     if we have been passed any custom parameter values, those
c     take priority

      do 200 ipair = 1, ncpair
         itype = ityc2p(ipair,1)
         jtype = ityc2p(ipair,2)
         if (usec2p(ipair,1)) then
            b1(itype,jtype) = valc2p(ipair,1)
            b1(jtype,itype) = b1(itype,jtype)
         endif
         if (usec2p(ipair,2)) then
            beta1(itype,jtype) = valc2p(ipair,2)
            beta1(jtype,itype) = beta1(itype,jtype)
         endif
         if (usec2p(ipair,3)) then
            b2(itype,jtype) = valc2p(ipair,3)
            b2(jtype,itype) = b2(itype,jtype)
         endif
         if (usec2p(ipair,4)) then
            beta2(itype,jtype) = valc2p(ipair,4)
            beta2(jtype,itype) = beta2(itype,jtype)
         endif
         if (usec2p(ipair,5)) then
            b3(itype,jtype) = valc2p(ipair,5)
            b3(jtype,itype) = b3(itype,jtype)
         endif
         if (usec2p(ipair,6)) then
            beta3(itype,jtype) = valc2p(ipair,6)
            beta3(jtype,itype) = beta3(itype,jtype)
         endif
         if (usec2p(ipair,7)) then
            capa(itype,jtype) = valc2p(ipair,7)
            capa(jtype,itype) = capa(itype,jtype)
         endif
         if (usec2p(ipair,8)) then
            alfa(itype,jtype) = valc2p(ipair,8)
            alfa(jtype,itype) = alfa(itype,jtype)
         endif
         if (usec2p(ipair,9)) then
            capq(itype,jtype) = valc2p(ipair,9)
            capq(jtype,itype) = capq(itype,jtype)
         endif
         if (usec2p(ipair,10)) then
            repmin(itype,jtype) = valc2p(ipair,10)
            repmin(jtype,itype) = repmin(itype,jtype)
         endif
         if (usec2p(ipair,11)) then
            repmax(itype,jtype) = valc2p(ipair,11)
            repmax(jtype,itype) = repmax(itype,jtype)
         endif
 200  continue

      do 300 ictrip = 1, nctrip
         itype = ityc3p(ictrip,1)
         jtype = ityc3p(ictrip,2)
         ktype = ityc3p(ictrip,3)
         if (usec3p(ictrip,1)) then
            glamda(itype,jtype,ktype) = valc3p(ictrip,1)
            glamda(ktype,jtype,itype) = glamda(itype,jtype,ktype)
         endif
 300  continue

      return
      end
c

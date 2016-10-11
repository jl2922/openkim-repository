c-----------------------------------------------------------------------
c     Jinit initializes some lookup tables for the J(r) function.
c     Ported from sjs' Jinit. Following that routine, it uses Steve
c     Rick's 1-dimensional integral representation of J(r), since it is
c     available for 3rd shell atoms, rather than Joel Bader's analytical
c     expressions, which would be more accurate but have not been 
c     derived (yet) for 3rd row atoms....sjs 8/2/07
c-----------------------------------------------------------------------
c     
      subroutine Jinit()

      include 'common_files.inc'
      include 'common_pots.inc'

c     local variables

      integer itype, itab, jtype
      real*8  dtspln, dvv, rc, tee, tsplin, twidth, vv, zblau, zi, zj

      write(istout, *) 'initializing J(r) data...'

c     generate lookup tables for J(r) Coulomb interaction

      do 310 itype = 1, ntypes
         do 300 jtype = itype, ntypes
            drj(itype,jtype) = rjmax(itype,jtype) / dble(njtab - 2)
            drj(jtype,itype) = drj(itype,jtype)

c           for atom types that have no zeta parameter, just
c           initialize the J(r) interaction to zero

            if (.not. iseety(itype) .or. .not. iseety(jtype)) then
               do 100 itab = 0, njtab
                  tabj(itype,jtype,itab) = 0.d0
                  tabj(jtype,itype,itab) = tabj(itype,jtype,itab)
 100           continue
            else

c     atomic number

               zi = kt2(itype)
               zj = kt2(jtype)

c     ZBL universal screening length

               zblau = zblcof * bohr / (zi**0.23d0 + zj**0.23d0)

               do 200 itab = 0, njtab
                  rc = dble(itab) * drj(itype,jtype)

c     switching function for ZBL <-> Coulomb interaction
                  
                  if (rc .gt. repmax(itype,jtype)) then
                     tsplin = 0.d0
                     dtspln = 0.d0
                  else if (rc .gt. repmin(itype,jtype)) then
                     twidth = repmax(itype,jtype) - repmin(itype,jtype)
                     tee = (rc - repmin(itype,jtype)) / twidth
                     tsplin = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
                     dtspln = 6.d0 * tee * (1.d0 - tee) / twidth
                  else
                     tsplin = 1.d0
                     dtspln = 0.d0
                  endif

                  vv = quadj(itype,jtype,rc)
                  dvv = quaddj(itype,jtype,rc)
                  
c     tabj holds J(r).
c     tabdj holds dJ/dr.
c     this includes the effect of the ZBL switching function

c     V = (1-S) * J

                  tabj(itype,jtype,itab) = (1.d0 - tsplin) * vv
                  tabj(jtype,itype,itab) = tabj(itype,jtype,itab)
                  tabdj(itype,jtype,itab) = (1.d0 - tsplin) * dvv
     .                 - dtspln * vv
                  tabdj(jtype,itype,itab) = tabdj(itype,jtype,itab)
 200           continue
            endif
 300     continue
 310  continue

      write(isterr, *) '                         ...done'

      return
      end
c
c-----------------------------------------------------------------------
c     quadj calculates the J(r) function between two atoms at a
c     specified distance, using quadrature. Integral due to Joel Bader.
c     Ported from sjs' quadj()....sjs 8/1/07
c-----------------------------------------------------------------------
c
      real*8 function quadj(itype, jtype, r)

      include 'common_files.inc'
      include 'common_pots.inc'

c     itype = atom type for first atom of pair
c     jtype = atom type for second atom of pair
c     r     = distance between atoms

      integer itype
      integer jtype
      real*8  r

c     local variables

      integer ik, ishell, jshell, kshell, lshell
      real*8  dk, rk, zetak, zetal

      ishell = ivalsh(itype)
      jshell = ivalsh(jtype)

      if (ishell .eq. 0) then
         write(isterr, *)
     .        'quadj: ERROR! valence shell unknown for atom type ', 
     .        itype
         call ioclos
         stop
      endif

      if (jshell .eq. 0) then
         write(isterr, *)
     .        'quadj: ERROR! valence shell unknown for atom type ',
     .        jtype
         call ioclos
         stop
      endif

c     make sure the smaller one comes first

      if (jshell .lt. ishell) then
         kshell = jshell
         lshell = ishell
         zetak = zeta(jtype)
         zetal = zeta(itype)
      else
         kshell = ishell
         lshell = jshell
         zetak = zeta(itype)
         zetal = zeta(jtype)
      endif

c     perform a semi-infinite 1-d integral to get J(r)
c     result is in units of eV/e^2.

      dk = rjkmax / njk
      rk = 0.d0
      quadj = 0.5d0 * rjkern(kshell, lshell, zetak, zetal, r, rk)
      do 100 ik = 1, njk - 1
         rk = dk * ik
         quadj = quadj + rjkern(kshell, lshell, zetak, zetal, r, rk)
 100  continue
      rk = dk * njk
      quadj = quadj
     .     + 0.5d0 * rjkern(kshell, lshell, zetak, zetal, r, rk)
      quadj = quadj * 2.d0 / pi * dk * epsinv

      return
      end
c
c-----------------------------------------------------------------------
c     rjkern calculates the kernel of the integral needed to evaluate
c     J(r) by quadrature.  The functions used depend on the principal
c     quantum number of the Slater orbitals being used. Ported from
c     prior sjs version....sjs 8/1/07
c-----------------------------------------------------------------------
c     
      real*8 function rjkern(ishell, jshell, zetai, zetaj, r, rk)

      include 'common_files.inc'

c     ishell = principal quantum number of first orbital
c     jshell = principal quantum number of second orbital
c     zetai  = orbital exponent of first orbital
c     zetaj  = orbital exponent of second orbital
c     r      = interatomic distance
c     rk     = value of integration variable

      integer ishell
      integer jshell
      real*8  zetai
      real*8  zetaj
      real*8  r
      real*8  rk
      
c     local variables

      real*8 rf, rfac, temp

      temp = r * rk
      if (temp .lt. 2.d-8) then
         rjkern = 1.d0
      else
         rjkern = sin(temp) / temp
      endif
      temp = rk / zetai
      rf = 1.d0 / (1.d0 + (0.5d0 * temp) ** 2) ** 2
      rfac = temp ** 2 * sqrt(rf)
      if (ishell .eq. 1) then
         rjkern = rjkern * rf
      else if (ishell .eq. 2) then
         rjkern = rjkern * rf
     .        * (1.d0 + rfac * (-3.d0 / 4.d0 + 1.d0 / 8.d0 * rfac))
      else if (ishell .eq. 3) then
         rjkern = rjkern * rf
     .        * (1.d0
     .           + rfac * (-11.d0 / 6.d0
     .                     + rfac * (17.d0 / 16.d0
     .                               + rfac * (-1.d0 / 4.d0
     .                                         + 1.d0 / 48.d0 * rfac))))
      else
         write(isterr, *)
     .        'rjkern: ERROR! cannot handle Slater shell of ', ishell
      endif

      temp = rk / zetaj
      rf = 1.d0 / (1.d0 + (0.5d0 * temp) ** 2) ** 2
      rfac = temp ** 2 * sqrt(rf)
      if (jshell .eq. 1) then
         rjkern = rjkern * rf
      else if (jshell .eq. 2) then
         rjkern = rjkern * rf
     .        * (1.d0 + rfac * (-3.d0 / 4.d0 + 1.d0 / 8.d0 * rfac))
      else if (jshell .eq. 3) then
         rjkern = rjkern * rf
     .        * (1.d0
     .           + rfac * (-11.d0 / 6.d0
     .                     + rfac * (17.d0 / 16.d0
     .                               + rfac * (-1.d0 / 4.d0
     .                                         + 1.d0 / 48.d0 * rfac))))
      else
         write(isterr, *)
     .        'rjkern: ERROR! cannot handle Slater shell of ', jshell
         call ioclos
         stop
      endif

      return
      end
c-----------------------------------------------------------------------
c     quaddj calculates dJ/dr by finite difference, with J(r) calculated
c     from quadj(), i.e. by quadrature. Ported from sjs' quaddj.
c     .....sjs 8/12/07
c-----------------------------------------------------------------------
c     This could just as easily be done by quadrature, or analytically,
c     but it's not.
c-----------------------------------------------------------------------
c
      real*8 function quaddj(itype, jtype, r)

      include 'common_files.inc'

c     itype = atom type of first atom
c     jtype = atom type of second atom
c     r     = distance between atoms

      integer itype
      integer jtype
      real*8  r

c     local variables

      dr = 1.d-5

      if (r .eq. 0.d0) then
         quaddj = 0.d0
      else
         quaddj = (quadj(itype,jtype,r+dr) - quadj(itype,jtype,r-dr))
     .        * 0.5d0 / dr
      endif

      return
      end
c
c-----------------------------------------------------------------------
c     J calculates the position-dependent part of the Coulomb
c     interaction for two partial charges.  Ported from sjs' J().
c     J(r) has units of eV/e^2....sjs 8/13/07
c-----------------------------------------------------------------------
c     
      real*8 function J(itype, jtype, r)

      include 'common_files.inc'
      include 'common_pots.inc'

c     itype = atom type of first atom
c     jtype = atom type of second atom
c     r     = distance between atoms

      integer itype
      integer jtype
      real*8  r

c     local variables

      integer itab
      real*8  dr, tee, tsplin, twidth, zblau, zi, zj
      
      zi = kt2(itype)
      zj = kt2(jtype)

c     ZBL universal screening length

      zblau = zblcof * bohr / (zi**0.23d0 + zj**0.23d0)

c     switching function for ZBL <-> Coulomb interaction

      if (r .gt. repmax(itype,jtype)) then
         tsplin = 0.d0
      else if (r .gt. repmin(itype,jtype)) then
         twidth = repmax(itype,jtype) - repmin(itype,jtype)
         tee = (r - repmin(itype,jtype)) / twidth
         tsplin = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
      else
         tsplin = 1.d0
      endif

      if (r .lt. rjmax(itype,jtype)) then

         if (ljdir) then

c           calculate J(r) directly, and apply switching function

            J = quadj(itype, jtype, r) * (1.d0 - tsplin)

         else

c           use linear interpolation in the J(r) lookup table, which
c           already has the switching function applied

            itab = int(r / drj(itype,jtype))
            dr = r - drj(itype,jtype) * itab
            J = tabj(itype,jtype,itab)
     .           + dr / drj(itype,jtype)
     .             * (tabj(itype,jtype,itab+1)-tabj(itype,jtype,itab))

         endif

      else

c        1/r is good enough. also apply switching function
         
         J = epsinv / r * tsplin

      endif

      return
      end
c
c-----------------------------------------------------------------------
c     dJdr calculates the distance derivative of the Coulomb interaction
c     for two atoms at a specified distance, as calculated in J().
c     It is evaluated from a lookup table, or as -1/r^2 beyond a cutoff.
c     Ported from SJS' dJdr()....sjs 8/15/07
c-----------------------------------------------------------------------
c     
      real*8 function dJdr(itype, jtype, r)

      include 'common_files.inc'
      include 'common_pots.inc'

c     itype = atom type for first atom
c     jtype = atom type for second atom
c     r     = distance between atoms

      integer itype
      integer jtype
      real*8  r

c     local variables

      integer itab
      real*8  dr, dtspln, tee, tsplin, twidth, zblau, zi, zj

      if (r .lt. rjmax(itype,jtype)) then

         zi = kt2(itype)
         zj = kt2(itype)

c     ZBL universal screening length

         zblau = zblcof * bohr / (zi**0.23d0 + zj**0.23d0)

c     switching function for ZBL <-> Coulomb interaction

         if (r .gt. repmax(itype,jtype)) then
            tsplin = 0.d0
            dtspln = 0.d0
         else if (r .gt. repmin(itype,jtype)) then
            twidth = repmax(itype,jtype) - repmin(itype,jtype)
            tee = (r - repmin(itype,jtype)) / twidth
            tsplin = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
            dtspln = 6.d0 * tee * (1.d0 - tee) / twidth
         else
            tsplin = 1.d0
            dtspln = 0.d0
         endif

         if (ljdir) then

c           calculate dJ/dr directly, and apply switching function

            dJdr = quaddj(itype,jtype,r) * (1.d0 - tsplin)

         else

c           use linear interpolation in the lookup table, which already has
c           the switching function applied

            itab = int(r / drj(itype,jtype))
            dr = r - drj(itype,jtype) * itab
            dJdr = tabdj(itype,jtype,itab)
     .           + dr / drj(itype,jtype)
     .             * (tabdj(itype,jtype,itab+1)-tabdj(itype,jtype,itab))

         endif
         
      else

c        -1/r^2 is good enough. also apply switching function
         
         dJdr = -epsinv / (r * r) * tsplin + epsinv / r * dtspln

      endif

      return
      end
c

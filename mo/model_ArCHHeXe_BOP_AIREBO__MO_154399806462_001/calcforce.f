c-----------------------------------------------------------------------
c     ezpotl is a wrapper for the calcforce subroutine that sets many of
c     the flags for the default cases of REBO and AIREBO potentials.
c     ...sjs 6/30/08
c-----------------------------------------------------------------------
c
      subroutine ezpotl(np, r0, lpbc, cube, mypot, fint, uu, tote)

      include 'common_files.inc'
      include 'common_pots.inc'


c     np          = number of particles
c     r0(a,d)     = position of atom a in dimension d
c     lpbc(d)     = whether to use periodic boundary conditions in dimension d
c     cube(d)     = box length in dimension d
c     mypot       = 1 for REBO, 2 for AIREBO
c     fint(a,d)   = internal force on atom a in dimension d (output)
c     uu          = internal virial (output)
c     tote        = potential energy (partial) (output)

      integer np
      real*8  r0(npmax,ndim)
      logical lpbc(ndim)
      real*8  cube(ndim)
      integer mypot
      real*8  fint(npmax,3)
      real*8  uu(ndim,ndim)
      real*8  tote

c     local variables

      integer atomof(nbcmax,2),
     .     bctype(nbcmax),
     .     iatom, igfunc, ipot, nbc
      logical lbc, les, lewald, lewsrf, llj, lnones, lqsolv,
     .     ltors
      real*8 bcchi(nbcmax), bq0(nbcmax),
     .     chi(npmax), q0(npmax),
     .     Efield(ndim)

      if (mypot .eq. 1) then

c     set varables for REBO

         llj = .false.
         ltors = .false.
         igfunc = isplin

      else if (mypot .eq. 2) then

c     set variables for AIREBO '00
         
         llj = .true.
         ltors = .true.
         igfunc = isplin

      else if (mypot .eq. 3) then

         llj = .true.
         ltors = .true.
         igfunc = ibop

      else
         write(isterr, *) "ezpotl: Potential ", mypot, " unrecognized."
         write(isterr, *) "Please use 1 for REBO, "
         write(isterr, *) "           2 for AIREBO (2000) or"
         write(isterr, *) "           3 for polyatomic AIREBO"
         call ioclos
         stop
      endif

c     these are correct for both REBO and AIREBO

      do 100 iatom = 1, np
         q0(iatom) = 0.d0
 100  continue
      lbc = .false.
      nbc = 0
      do 200 idim = 1, ndim
         Efield(idim) = 0.d0
 200  continue
      ipot = 1
      les = .false.
      lqsolv = .false.
      lnones = .true.
      lewald = .false.
      lewsrf = .false.

      call calcforce(np, r0, q0, lbc, nbc, bq0, atomof, bctype, Efield,
     .     lpbc, cube, ipot, igfunc, 
     .     llj, ltors, les, lqsolv, lnones, lewald,
     .     lewsrf, fint, chi, bcchi, uu, tote)

      return
      end
c
c-----------------------------------------------------------------------
c     calcforce evaluates the force and energy. Adapted by yl from sjs'
c     model.
c-----------------------------------------------------------------------
c
      subroutine calcforce(np, r0, q0, lbc, nbc, bq0, atomof, bctype,
     .     Efield, lpbc, cube,
     .     ipot, igfunc, 
     .     llj, ltors, les, lqsolv, lnones, lewald, lewsrf,
     .     fint, chi, bcchi, uu, tote) 

      include 'common_files.inc'
      include 'common_pots.inc'

c     np          = number of particles (input)
c     r0(a,d)     = position of atom a in dimension d (input)
c     q0(a)       = charge on atom a (input)
c     lbc         = whether there are any bond charges (input)
c     nbc         = number of bond charges (input)
c     bq0(b)      = bond charge on bond b (input)
c     atomof(b,n) = nth (n=1,2) atom of bond charge b (input)
c     bctype(b)   = type (sigma, pi) of bond charge b (input)
c     Efield(d)   = external electric field in dimension d (V/A) (input)
c     lpbc(d)     = whether to use periodic bdy conds in dimension d (input)
c     cube(d)     = box length in dimension d (input)
c     ipot        = 1 for REBO (input)
c     igfunc      = style of angular term (input)
c     llj         = whether to evaluate van der Waals interactions (input)
c     ltors       = whether to evaluate torsion interactions (input)
c     les         = whether to evaluate electrostatics interactions (input)
c     lqsolv      = whether charges will be solved by matrix methods (input)
c     lnones      = whether to evaluate nonelectrostatics interactions (input)
c     lewald      = whether to use Ewald sums (input)
c     lewsrf      = whether to use the Ewald surface term (input)
c     fint(a,d)   = internal force on atom a in dimension d (output)
c     chi(a)      = electronegativity (dV/dq, neg chem pot) of atom a (output)
c     bcchi(b)    = electroneg (dV/dq, neg electrochem pot) of bond b (output)
c     uu          = internal virial (output)
c     tote        = potential energy (partial) (output)

      integer np
      real*8  r0(npmax,ndim)
      real*8  q0(npmax)
      logical lbc
      integer nbc
      real*8  bq0(nbcmax)
      integer atomof(nbcmax,2)
      integer bctype(nbcmax)
      real*8  Efield(ndim)
      logical lpbc(ndim)
      real*8  cube(ndim)
      integer ipot
      integer igfunc
      logical llj
      logical ltors
      logical les
      logical lqsolv
      logical lnones
      logical lewald
      logical lewsrf
      real*8  fint(npmax,3)
      real*8  chi(npmax)
      real*8  bcchi(nbcmax)
      real*8  uu(ndim,ndim)
      real*8  tote

c     local variables

      integer iatom, idim, jdim, nrebo

c     rebuild the pair list if needed

      call prpcov(np, cube, r0, lpbc)

c     zero the potential energy

      tote = 0.0d0

      if (lnones) then

c        zero internal virial

         do 90 idim = 1, 3
            do 80 jdim = 1, 3
               uu(idim,jdim) = 0.d0
 80         continue
 90      continue

c        zero forces

         do 110 idim = 1, 3
            do 100 iatom = 1, np
               fint(iatom,idim) = 0.d0
 100        continue
 110     continue

c        add in the long-range LJ correction

         tote = tote + vlrc

ceatom
c$$$c     zero the individual atom energies
c$$$
c$$$      do 200 iatom = 1, np
c$$$         eatom(iatom) = 0.d0
c$$$ 200  continue
ceatom

c$$$cnrg
c$$$      eljcc = 0.d0
c$$$      eljch = 0.d0
c$$$      eljhh = 0.d0
c$$$      elj = 0.d0
c$$$      erbcc = 0.d0
c$$$      erbch = 0.d0
c$$$      erbhh = 0.d0
c$$$      erb = 0.d0
cnrg

c        number of atoms which use the REBO potential

         nrebo = 0
         do 200 itype = 1, ntypes
            if (lrebot(itype)) then
               nrebo = nrebo + noa(itype)
            endif
 200     continue

c     covalent bonding interactions

         if((ipot.eq.1).and.(nrebo.ne.0)) then
            call caguts(np, igfunc, fint, cube, ltors, r0, uu, tote)
         endif

c        van der Waals interactions

         if (llj) then
            call ljguts(np, igfunc, fint, cube, r0, uu, tote, lpbc)
         endif

      endif

c     electrostatic interactions

      if (les) then
         call esguts(np, r0, q0, lbc, nbc, bq0, atomof, bctype, 
     .        Efield, lpbc, cube,
     .        lqsolv, lewald, lewsrf, tote, fint, chi, bcchi, uu)
      endif

      return

      end
c

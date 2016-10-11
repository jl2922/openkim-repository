c-----------------------------------------------------------------------
c     ljpars sets those Lennard-Jones parameters, including lookup
c     tables, that can be set before knowing the simulation conditions
c     (such as the buffer distance)....sjs 1/7/11
c-----------------------------------------------------------------------
c
      subroutine ljpars(nctype, ityc1p, usec1p, valc1p, lckfil)

      include 'common_files.inc'
      include 'common_pots.inc'

c     nctype      = number of atom types with custom singleton parameters
c     ityc1p(i)   = ith atom type with custom singleton parameters
c     usec1p(i,k) = whether to use a custom value of par. k for ith atom type
c     valc1p(i,k) = custom value of parameter k for ith atom type
c     lckfil      = whether to check files for custom parameter values

      integer nctype
      integer ityc1p(maxc1t)
      logical usec1p(maxc1t,maxc1p)
      real*8  valc1p(maxc1t,maxc1p)
      logical lckfil

c     local variables

      character name*20, unit*20
      integer   i, itype, j, jtype, k, kmax, nsugg
      real*8    dr, dslj, dvlj, r, r6, rmax, rmin, rsqs,
     .     sigmin, sigwid, slj, swidth, tee, vlj

      do 90 i = 1, ntypes
         do 80 j = 1, ntypes
            sig(i,j) = 0.0d0
            eps(i,j) = 0.0d0
 80      continue
 90   continue

      sig(ihyd,ihyd) = 2.65d0
      eps(ihyd,ihyd) = 17.4d0

      sig(ihel,ihel) = 2.28d0
      eps(ihel,ihel) = 10.2d0

c     these are the C values from Liegi's fit, for the FC potential.
c     but we don't want to use them, as the C values should be the 
c     same in the FC potential as in the HC potential, if we want to
c     get graphite right

c      sig(icarb,icarb) = 3.02d0
c      eps(icarb,icarb) = 24.8d0

      sig(icarb,icarb) = 3.40d0
      eps(icarb,icarb) = 33.d0

c     totally random parameters for oxygen

      sig(ioxy,ioxy) = 2.75d0
      eps(ioxy,ioxy) = 30.d0

c     from Liegis fit

      sig(ifluor,ifluor) = 2.49d0
      eps(ifluor,ifluor) = 27.2d0

c     draft parameters for silicon, taken from JCP 117 1804 '02, which 
c     apparently uses UFF
c     (Glass Phys. Chem. 33 86 '07 quotes sig = 3.7, eps = 137.1)

      sig(isili,isili) = 3.826d0
      eps(isili,isili) = 202.45d0

      sig(iargon,iargon) = 3.41d0
      eps(iargon,iargon) = 119.8d0

c     approximate xenon parameters.  reference and better values needed.

      sig(ixenon,ixenon) = 4.06d0
      eps(ixenon,ixenon) = 230.d0

      if (lckfil) then

c     The values above can be overridden by custom parameter files

         name = 'sigma'
         unit = 'Angstrom'
         call cuspa1(sig, name, unit)

         name = 'epsilon'
         unit = 'K'
         call cuspa1(eps, name, unit)

      endif

c     if we have been passed any custom parameter values, those take
c     priority

      do 95 ictype = 1, nctype
         itype = ityc1p(ictype)
         if (usec1p(ictype,1)) then
            sig(itype,itype) = valc1p(ictype,1)
         endif
         if (usec1p(ictype,2)) then
            eps(itype,itype) = valc1p(ictype,2)
         endif
 95   continue

c     convert epsilon values from Kelvin to eV, and carry as 4 x epsilon

      do 100 itype = 1, ntypes
         eps(itype,itype) = 4.0d0 * eps(itype,itype) * fK2eV
 100  continue
      
c     use either Lorentz-Berthelot or OPLS combining rules for
c     mixed LJ parameters.  also set the value of rlj0 and rlj1 (the
c     inner and outer edges of the absolute LJ switch) using sig
c     and set the squares of these switching function parameters

      if (lorber) then
         do 210 itype = 1, ntypes
            do 200 jtype = 1, ntypes
               eps(itype,jtype) = 
     .              sqrt(eps(itype,itype) * eps(jtype,jtype))
               sig(itype,jtype) = 
     .              (sig(itype,itype) + sig(jtype,jtype)) / 2.d0
               rlj0(itype,jtype) = sig(itype,jtype)
               rlj1(itype,jtype) = 2.d0 ** (1.d0/6.d0) * 
     .              sig(itype,jtype)
 200        continue
 210     continue
      else
         do 230 itype = 1, ntypes
            do 220 jtype = 1, ntypes
               eps(itype,jtype) = 
     .              sqrt(eps(itype,itype) * eps(jtype,jtype))
               sig(itype,jtype) = 
     .              sqrt(sig(itype,itype) * sig(jtype,jtype))
               rlj0(itype,jtype) = sig(itype,jtype)
               rlj1(itype,jtype) = 2.d0 ** (1.d0/6.d0) *
     .              sig(itype,jtype)
 220        continue
 230     continue
      endif

c     interactions with pure LJ atoms like Ar should NOT have a switching
c     range to turn off LJ interactions.  Setting rlj0 and rlj1 to zero
c     guarantees LJ interactions are always on.

      do 240 itype = 1, ntypes
	 rlj0(ihel,itype) = 0.d0
	 rlj1(ihel,itype) = 0.d0
	 rlj0(itype,ihel) = 0.d0
	 rlj1(itype,ihel) = 0.d0
         rlj0(iargon,itype) = 0.d0
         rlj1(iargon,itype) = 0.d0
         rlj0(itype,iargon) = 0.d0
         rlj1(itype,iargon) = 0.d0
         rlj0(itype,ixenon) = 0.d0
         rlj1(itype,ixenon) = 0.d0
         rlj0(ixenon,itype) = 0.d0
         rlj1(ixenon,itype) = 0.d0
 240  continue

c     set the boundaries (in bij-space) for the switch on the 
c     contingent LJ interaction:

c     default parameters

      do 310 itype = 1, ntypes
         if (lrebot(itype)) then
            bijmin(itype,itype) = 0.5d0
            bijmax(itype,itype) = 0.9d0
            do 300 jtype = itype+1, ntypes
               if (lrebot(jtype)) then
                  bijmin(itype,jtype) = 0.5d0
                  bijmax(itype,jtype) = 0.9d0

                  bijmin(jtype,itype) = bijmin(itype,jtype)
                  bijmax(jtype,itype) = bijmax(itype,jtype)
               endif
 300        continue
         endif
 310  continue

c     parameters that have actually been fitted

      bijmin(icarb,icarb) = 0.77d0
      bijmax(icarb,icarb) = 0.81d0

      bijmin(icarb,ihyd) = 0.75d0
      bijmax(icarb,ihyd) = 0.90d0

      bijmin(ihyd,icarb) = bijmin(icarb,ihyd)
      bijmax(ihyd,icarb) = bijmax(icarb,ihyd)

      bijmin(ihyd,ihyd) = 0.32d0
      bijmax(ihyd,ihyd) = 0.42d0

c     neighbor list boundaries:

c     hardwire switching function width as 0.84 sigma (wide enough that
c     2nd derivative is never positive)

      sigcut = 3.d0
      sigwid = 0.84
      sigin  = sigcut - sigwid

      do 330 i=1,ntypes
         do 320 j=1,ntypes

c     rmaxlj is the outer edge of the LJ switching function
c     rminlj is the inner edge of the LJ switching function
c     r2mxlj is the square of the outer edge of the LJ switching function
c     sig gets carried as sigma squared

            if(eps(i,j).ne.0.0d0) then

               rmaxlj(i,j) = sigcut * sig(i,j)
               rminlj(i,j) = sigin * sig(i,j)
               r2mxlj(i,j) = rmaxlj(i,j) * rmaxlj(i,j)
               sig(i,j) = sig(i,j)*sig(i,j)
            else
               rmaxlj(i,j) = 0.d0
               rminlj(i,j) = 0.d0
               r2mxlj(i,j) = 0.d0
            endif

 320     enddo
 330  enddo

      if (.not. lljdir) then
         dellj = 0.001d0
         do 440 i = 1, ntypes
            do 430 j = 1, ntypes
               kmax = int(sqrt(r2mxlj(i,j)) / dellj)
               if (kmax + 1 .gt. maxbin) then
                  nsugg = kmax + 1
                  write(isterr, *) 'ljparam: LJ switching region for ',
     .                 'atom types ', i, ' and ', j, 'too wide ',
     .                 'for ', dellj, ' A bin width'
                  write(isterr, *) 'increase maxbin to at least ',
     .                 nsugg, ' and recompile'
                  call ioclos
                  stop
               endif
               if (kmax .ge. 0) then
                  do 410 k=kmax+1, maxbin
                     vljtb(k,i,j) = 0.d0
                     dvljtb(k,i,j) = 0.d0
 410              continue
               endif

               do 420 k = kmax, 2, -1
                  r = (k - 1) * dellj
                  if (((i .gt. 2) .or. (j .gt. 2)) .and. (r .eq. 0.d0))
     .                 r = dellj
                  rsqs = r ** 2

                  call ljpair(i, j, rsqs, r, vljtb(k,i,j),
     .                 dvljtb(k,i,j))

 420           continue

c     the LJ potential is infinite at r=0. approximate this by a big
c     number.  the derivative is negative infinity.

               vljtb(1,i,j) = bigpos
               dvljtb(1,i,j) = bigneg

 430        continue
 440     continue
      endif

      return

      end
c
c-----------------------------------------------------------------------
c     ljcuts sets any Lennard-Jones-related variables that can't be
c     defined until the simulation conditions are known....sjs 1/7/11
c-----------------------------------------------------------------------
c     
      subroutine ljcuts(rbuffr)

      include 'common_files.inc'
      include 'common_pots.inc'

      real*8 rbuffr

c     local variables

      integer itype, jtype

c     neighbor list boundaries:

c     rslj is (the square of) the outer LJ neighbor list cutoff

      do 110 itype = 1, ntypes
         do 100 jtype = 1, ntypes

            if(eps(itype,jtype) .ne. 0.0d0) then
               rslj(itype,jtype) = (rmaxlj(itype,jtype) + rbuffr) ** 2
            else
               rslj(itype,jtype) = 0.0d0
            endif

 100     enddo
 110  enddo

      return
      end

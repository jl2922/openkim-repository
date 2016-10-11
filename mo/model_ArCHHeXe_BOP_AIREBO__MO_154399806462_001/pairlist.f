c-----------------------------------------------------------------------
c     chkprs checks to see whether the (REBO) pair list has expired.
c     it should be called whenever the positions are updated.
c     ...sjs 10/9/96
c-----------------------------------------------------------------------
c
      subroutine chkprs(np, cube, r0, lpbc)
      
      include 'common_files.inc'
      include 'common_pots.inc'

c     np      = number of particles
c     cube(d) = box length in dimension d
c     r0(a,d) = position of atom a in dimension d
c     lpbc(d) = whether or not to apply periodic boundary conditions in 
c               dimension d

      integer np
      real*8  cube(ndim)
      real*8  r0(npmax,ndim)
      logical lpbc(ndim)

c     local variables

      integer iatom

c     our knowledge about the pair list will soon be fresh again

      lposol = .false.

      do 110 iatom = 1, np
         rsq = 0.d0
         do 100 idim = 1, ndim
            r = r0(iatom,idim) - r0l(iatom,idim)
            if (lpbc(idim)) then
               r = r - cube(idim) * anint(r / cube(idim))
            endif
            rsq = rsq + r * r
 100     continue
         if (rsq .gt. prtrig) then
            go to 200
         endif
 110  continue

c     neighbor lists are still fresh

      return

c     neighbor lists need to be rebuilt

 200  continue

      lcovol = .true.
      lvdwol = .true.
      return

      end
c     
c-----------------------------------------------------------------------
c     fakepr generates information for a fake REBO i,j pair based on the
c     real LJ pair passed to it.  in particular, it sets ihalf, jhalf,
c     rij, rcor, fij, dww, and lrebop.  the fake pair is numbered one 
c     higher than the current number of real REBO pairs in the pair 
c     list, and is arbitrarily set at a fixed interatomic separation.  
c     this is used in one of the switching functions for the contingent 
c     LJ interaction....sjs 5/29/97
c-----------------------------------------------------------------------
c
      subroutine fakepr(iprlst, ijljpr, rijsep, rijvec)

      include 'common_files.inc'
      include 'common_pots.inc'

      dimension rijvec(3)

c     identify the atoms in this pair from the LJ pair list

      if (iprlst .eq. iljsaf) then
         iatom = iv(ijljpr)
         jatom = jv(ijljpr)
      else
         iatom = iljwv(ijljpr)
         jatom = jljwv(ijljpr)
      endif

      itype = iat2ty(iatom)
      jtype = iat2ty(jatom)

c     set the pair number in the REBO pair list

      ijrbpr = npairs + 1

c     make the atoms identifiable from REBO routines

      ihalf(ijrbpr) = iatom
      jhalf(ijrbpr) = jatom

c     target pair separation is just inside the REBO switch

      rfixed = dijmin(itype,jtype)

      rscale = rfixed / rijsep

c     scale rij and rcor down to have the desired separation

      do 300 idim = 1,3
         rijv(ijrbpr,idim) = rijvec(idim) * rscale
 300  continue

      rcor(ijrbpr) = rfixed

      lrebop(ijrbpr) = .true.

      fij(ijrbpr) = 1.d0
      dww(ijrbpr) = 0.d0

      return
      end
c
c-----------------------------------------------------------------------
c     fillpr calculates several variables for a given i,j pair.
c     ...sjs 5/26/97
c-----------------------------------------------------------------------
c     added construction of double-length REBO pair list....sjs 7/24/98
c-----------------------------------------------------------------------
c
      subroutine fillpr(ipair, cube, r0, lpbc)

      include 'common_files.inc'
      include 'common_pots.inc'

      dimension cube(ndim)
      dimension r0(npmax,ndim)
      logical lpbc(ndim)

c     local variables

      real*8 rvec(ndim),
     .     dvrdr, vr

c     look up which atoms are in this pair
      
      i = ihalf(ipair)
      j = jhalf(ipair)

      ki = iat2ty(i)
      kj = iat2ty(j)

c     pair separation between the two atoms rij = r_i - r_j

      rsq = 0.d0
      do 200 idim = 1, ndim
         rvec(idim) = r0(i,idim) - r0(j,idim)
         if (lpbc(idim)) then
            rvec(idim) = rvec(idim) - cube(idim) * 
     .           anint(rvec(idim)/cube(idim))
         endif
         rsq = rsq + rvec(idim) * rvec(idim)
         rijv(ipair,idim) = rvec(idim)
 200  continue

      lrebop(ipair) = .false.

c     many of these variables are only used for interacting REBO pairs,
c     so can be skipped if the pair is not close enough to interact

      if (rsq .gt. dijmx2(ki,kj)) then
         return
      endif

c     rcor is the scalar r_ij distance

      rcor(ipair) = sqrt(rsq)

c      print *, "rij for i,j is=", rcor(ipair),i,j

      if (lrebot(ki) .and. lrebot(kj)) then

         lrebop(ipair) = .true.

c        build a short pair list, containing only actual REBO neighbors
c        this step, for use in looping over hydrocarbon or fluorocarbon
c        REBO pairs

         if (i .lt. j) then
            nchch = nchch + 1
            ichch(nchch) = ipair
         endif

c        also build a double-length pair list which can be used to
c        to find all actual REBO neighbors of a given atom

         n2chch = n2chch + 1
         i2chch(n2chch) = ipair
         n2strt(i+1) = n2chch + 1

      else
         lrebop(ipair) = .false.
      endif

      if (lpdir) then

c     pairwise terms should be obtained by direct evaluation

         call rebopr(ki, kj, rcor(ipair), fij(ipair), dww(ipair),
     .        vr, dvrdr, exx1(ipair), dexx1(ipair))

c     the V^R pairlist is non-redundant, so it is only filled
c     for half of the pairs

         if (i .lt. j) then
            repel(ipair) = vr
            drepel(ipair) = dvrdr
         endif

      else

c     evaluate the pairwise terms from lookup tables, using linear
c     interpolation:
c     for x <= x + dx < x + Dx, f(x+dx) = f(x) + [f(x+Dx)-f(x)]*dx/Dx.

c     rt can be at most rrebo0 [= sqrt(r2rb0)], and ddtab is 
c     rrebo0/(ntab-2), so the variable it can be no more than ntab-1, 
c     so the bounds-checking on it has been eliminated

         rt = rcor(ipair) / ddtab(ki,kj)
         it = int(rt) + 1

c     the first few of these are stored under both _ij and _ji, so they
c     can be looked up both times through the pair list

c     fij stores the switching function f_ij.  

         fij(ipair) = tabfc(ki,kj,it)
     .        + (tabfc(ki,kj,it+1) - tabfc(ki,kj,it)) * (rt-it+1)

c     dww stores the r_ij derivative of f_ij.

         dww(ipair) = tabdfc(ki,kj,it)
     .        + (tabdfc(ki,kj,it+1) - tabdfc(ki,kj,it)) * (rt-it+1)

c     exx1 stores (half of) the attractive REBO interaction V^A

         exx1(ipair) = atable(ki,kj,it) +
     .        (atable(ki,kj,it+1) - atable(ki,kj,it)) * (rt-it+1)

c     dexx1 stores the (half 1/r of the) r_ij derivative of V^A

         dexx1(ipair) = datable(ki,kj,it) +
     .        (datable(ki,kj,it+1) - datable(ki,kj,it)) * (rt-it+1)

c     the rest of these are stored in the lookup table only once, under
c     _ij, so we have to be careful to only look for them once

         if(i .ge. j) then
            return
         endif

c     repel stores the repulsive REBO interaction V^R

         repel(ipair) = rtable(ki,kj,it) +
     .        (rtable(ki,kj,it+1) - rtable(ki,kj,it)) * (rt-it+1)

c     drepel stores the r_ij derivative of V^R

         drepel(ipair) = drtable(ki,kj,it) +
     .        (drtable(ki,kj,it+1) - drtable(ki,kj,it)) * (rt-it+1)

      endif
      return
      end
c
c-----------------------------------------------------------------------
c     gljprs generates the LJ pair list
c-----------------------------------------------------------------------
c     modified to generate two LJ pair lists: one with pairs which are
c     within buffer distance(s) of being (1,2) or (1,3) neighbors (and
c     thus could become neighbors by the next pair list update) and
c     one with pairs which will not become (1,2) or (1,3) neighbors by
c     the next update....sjs 7/9/97
c-----------------------------------------------------------------------
c     watch list extended to include current and possible future (1,4)
c     pairs....sjs 7/23/98
c-----------------------------------------------------------------------
c     changed the algorithm for constructing (ordered) watch and safe
c     lists, now much faster....sjs 7/24/98
c-----------------------------------------------------------------------
c
      subroutine gljprs(np, cube, r0, lpbc)

      include 'common_files.inc'
      include 'common_pots.inc'

      integer np
      dimension cube(3)
      dimension r0(npmax,ndim)
      logical lpbc(ndim)

      logical lpair(npmax)

c     use the nabors REBO pairlist to construct the LJ watchlist,
c     consisting of all pairs which are, or could be, (1,2), (1,3),
c     or (1,4) neighbors before the next update.

      nljwpr = 0
      do 40 iatom = 1, np - 1
         do 5 jatom = 1, np
            lpair(jatom) = .false.
 5          continue
         ijbeg = nabors(iatom)
         ijend = nabors(iatom+1) - 1
         do 30 ijpair = ijbeg, ijend
            jatom = jhalf(ijpair)
            lpair(jatom) = .true.
            jkbeg = nabors(jatom)
            jkend = nabors(jatom+1) - 1
            do 20 jkpair = jkbeg, jkend
               katom = jhalf(jkpair)
               lpair(katom) = .true.
               klbeg = nabors(katom)
               klend = nabors(katom+1) - 1
               do 10 klpair = klbeg, klend
                  latom = jhalf(klpair)
                  lpair(latom) = .true.
 10          continue
 20       continue
 30    continue

c        now that we know which atoms are neighbors of iatom, place
c        them into the LJ watchlist in ascending order, and do not
c        double-count pairs

         do 35 jatom = iatom + 1, np
            if (lpair(jatom)) then
               nljwpr = nljwpr + 1
               if (nljwpr .gt. nwlmax) then
                  go to 9100
               endif
               iljwv(nljwpr) = iatom
               jljwv(nljwpr) = jatom
            endif
 35      continue
 40   continue

c     now loop over all pairs, checking to see if they are within 
c     buffer range of the LJ potential.  if so, double-check to make
c     sure they aren't already on the watch list, and add them to the
c     LJ safe pairlist.

      nljprs = 0
      iwatch = 1
      do 70 iatom = 1, np - 1
         do 60 jatom = iatom + 1, np
            itype = iat2ty(iatom)
            jtype = iat2ty(jatom)

c           skip out if there is no LJ interaction defined for this pair

            if (rslj(itype,jtype) .eq. 0.d0) then
               go to 60
            endif

c           calculate the pair distance, and skip out if it is longer
c           than the LJ cutoff plus the buffer distance

            rsqs = 0.d0
            do 50 idim = 1, 3
               rr = r0(iatom,idim) - r0(jatom,idim)
               if (lpbc(idim)) then
                  rr = rr - cube(idim) * anint(rr/cube(idim))
               endif
               rsqs = rsqs + rr * rr
               if (rsqs .gt. rslj(itype,jtype)) then
                  go to 60
               endif
 50         continue

c           we've found a valid LJ pair.  now scan through the ordered
c           watch list to see if we've already counted this pair
            
 51         continue
            if (iwatch .gt. nljwpr
     .           .or. iljwv(iwatch) .gt. iatom 
     .           .or. (iljwv(iwatch) .eq. iatom 
     .           .and. jljwv(iwatch) .gt. jatom)) then
               go to 52
            else if (iljwv(iwatch) .eq. iatom 
     .              .and. jljwv(iwatch) .eq. jatom) then
               go to 60
            else
               iwatch = iwatch + 1
               go to 51
            endif
 52         continue
            
c           it wasn't on the watch list, so add it to the safe list

            nljprs = nljprs + 1
            if (nljprs .gt. nmabig) then
               go to 9000
            endif
            iv(nljprs) = iatom
            jv(nljprs) = jatom
 60      continue
 70   continue

      lvdwol = .false.
            
      return

 9000 nsugg = int(dble(np) / iatom * nmabig)
      write(isterr, *) 'gljprs: just found safe LJ pair #', nljprs,
     .     ' vs. max of ', nmabig
      write(isterr, *) 'solution = increase nmabig to about ', nsugg,
     .     ' and recompile'
      go to 9500

 9100 nsugg = int(dble(np) / iatom * nljwpr)
      write(isterr, *) 'gljprs: just found watch LJ pair #', nljwpr,
     .     ' vs. max of ', nwlmax
      write(isterr, *) 'solution = increase nwlmax to about ',
     .     nsugg, ' and recompile'
      go to 9500

 9500 call ioclos
      stop
      end
c
c-----------------------------------------------------------------------
c     genprs builds the REBO pair list....sjs 10/10/96
c-----------------------------------------------------------------------
c     unrolled some loops and un-array-ized some length-3 vectors.
c-----------------------------------------------------------------------
c     
      subroutine genprs(np, cube, r0, lpbc)

      include 'common_files.inc'
      include 'common_pots.inc'

      integer np
      dimension cube(3)
      dimension r0(npmax,ndim)
      logical lpbc(ndim)

c     prtemp is the max. number of REBO neighbors that an atom
c     should realistically be expected to have.

      integer prtemp
      parameter (prtemp = 25)

      dimension ipr(prtemp,npmax),
     .     npr(npmax)

      do 100 irba = 1, nrba
         npr(irba) = 0
 100  continue

      kpair = 0
      do 130 irba = 1, nrba
         iatom = irblst(irba)
         nabors(iatom) = kpair + 1

         itype = iat2ty(iatom)

         rix = r0(iatom,1)
         riy = r0(iatom,2)
         riz = r0(iatom,3)

c     don't use a full N**2 loop; just use N**2/2 and build a forward-
c     looking list of pairs for future atoms to use.  should just do 
c     away with the double pairs altogether.

         do 120 jrba = irba + 1, nrba
            jatom = irblst(jrba)
         
            jtype = iat2ty(jatom)

            dx = rix - r0(jatom,1)
            if (lpbc(1)) then
               dx = dx - cube(1) * anint(dx / cube(1))
            endif
            rsq = dx * dx
            if (rsq .gt. r2rbmx(itype,jtype)) then
               go to 120
            endif
            dx = riy - r0(jatom,2)
            if (lpbc(2)) then
               dx = dx - cube(2) * anint(dx / cube(2))
            endif
            rsq = rsq + dx * dx
            if (rsq .gt. r2rbmx(itype,jtype)) then
               go to 120
            endif
            dx = riz - r0(jatom,3)
            if (lpbc(3)) then
               dx = dx - cube(3) * anint(dx / cube(3))
            endif
            rsq = rsq + dx * dx
            if (rsq .gt. r2rbmx(itype,jtype)) then
               go to 120
            endif
            if (rsq .eq. 0.d0) then
               write(isterr, *) 
     .              'genprs: atoms ', iatom, ' and ', jatom, 
     .              ' have zero separation'
               call ioclos
               stop
            endif

            kpair = kpair + 1
            ihalf(kpair) = iatom
            jhalf(kpair) = jatom

c     remember this pair so that we don't have to recalculate it
c     the other way around.  this is somewhat kludgey and will
c     go away once the pair list is cut in half everywhere

            npr(jrba) = npr(jrba) + 1
            if (npr(jrba) .gt. prtemp) then
               write(isterr, *) 'genprs: atom ', jatom, ' has ',
     .              npr(jrba), ' REBO neighbors vs. limit of ', prtemp
               write(isterr, *) 
     .              '    you could increase prtemp and recompile, ',
     .              'but something else is probably wrong'
               call ioclos
               stop
            endif
            ipr(npr(jrba),jrba) = iatom

 120     continue

c     for the i > j pairs, don't recalculate distances but just
c     remember the pairs we found when looking at them the other
c     way around

         do 125 ipair = 1, npr(irba)
            kpair = kpair + 1
            jatom = ipr(ipair,irba)
            ihalf(kpair) = iatom
            jhalf(kpair) = jatom
 125     continue

 130  continue

c     so that the last atom knows how many neighbors he has:

      nabors(np+1) = kpair + 1
      npairs = kpair

c     this is a kludge to set the nabors pair list index for non-REBO
c     atoms to the same value as the next following REBO atom (so that
c     preceding REBO atoms can find out how many neighbors they have).
c     this can disappear once I switch all loops which refer
c     to nabors to be 1,nrba loops

      do 140 inra = 1, nnra
         iatom = inrlst(inra)
         do 135 iatom2 = iatom, np
            if (lrebot(iat2ty(iatom2))) then
               nabors(iatom) = nabors(iatom2)
               go to 138
            endif
 135     continue
         nabors(iatom) = nabors(np+1)
 138     continue
 140  continue

c     npairs must be strictly less than nlmax because sometimes there
c     is a fictitious pair in the last slot

      if (npairs .ge. nlmax) then
         write(isterr, *) 'genprs: npairs value of ', npairs,
     .        ' exceeds nlmax-1 value of ', nlmax-1
         write(isterr, *) '        increase nlmax and recompile'
         call ioclos
         stop
      endif

c     reset the last-pair-list-update positions

      do 210 idim = 1, 3
         do 200 iatom = 1, np
            r0l(iatom,idim) = r0(iatom,idim)
 200     continue
 210  continue
      lposol = .false.

      lcovol = .false.

      return
      end
c
c-----------------------------------------------------------------------
c     clustr generates information about clusters of atoms, or
c     molecules.  This routine needs to be called AFTER fillpr, so that
c     the REBO-only pair list can be used....sjs 2/21/01
c-----------------------------------------------------------------------
c     
      subroutine clustr(np, cube, r0, lpbc, lrxn)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np      = number of particles
c     cube(d) = size of periodic box in dimension d
c     r0(a,d) = position of atom a in dimension d
c     lpbc    = whether to use periodic boundary conditions
c     lrxn    = whether a reaction has occurred to change molecular membership

      integer np
      real*8 cube(ndim)
      real*8 r0(npmax,ndim)
      logical lpbc(ndim)
      logical lrxn

c     local variables

      integer omolec(npmax), 
     .     iatom, ijoin, imolec, iprind, isurv, jatom

c     remember the old molecular composition

      do 50 iatom = 1, np
         omolec(iatom) = molec(iatom)
 50   continue

c     starting setup:  each atom belongs to its own molecule

      nmolec = 0
      do 100 iatom = 1, np
         molec(iatom) = iatom
         inext(iatom) = -1
         natmol(iatom) = 1
         nmolec = nmolec + 1
 100  continue

c     make sure the covalent pair list is up to date, before using it

      call prpcov(np, cube, r0, lpbc)

c     loop over all bonds, combining molecules in the process

      do 200 iprind = 1, nchch
         ijpair = ichch(iprind)
         iatom = ihalf(ijpair)
         jatom = jhalf(ijpair)

c     these atoms are bonded, so make sure they belong to the same
c     molecule

         if (molec(iatom) .ne. molec(jatom)) then
            if (molec(iatom) .lt. molec(jatom)) then
               isurv = iatom
               ijoin = jatom
            else
               isurv = jatom
               ijoin = iatom
            endif

c     find the end of the surviving molecule

            itail = isurv
 110        continue
            if (inext(itail) .ne. -1) then
               itail = inext(itail)
               go to 110
            endif

c     tack the subsumed molecule onto the end of the surviving 
c     molecule.

            natmol(molec(isurv)) = natmol(molec(isurv)) 
     .           + natmol(molec(ijoin))
            natmol(molec(ijoin)) = 0

            inext(itail) = molec(ijoin)

 120        continue
            if (inext(itail) .ne. -1) then
               itail = inext(itail)
               molec(itail) = molec(isurv)
               go to 120
            endif

            nmolec = nmolec - 1

         endif
 200  continue

c     fill a lookup table so that we can loop over molecules easily
c     also keep track of some atom and molecule counts to make life
c     easier when we do molecule-neutral variable-charge electrostatics

c     npolat = number of atoms in polyatomic molecules
c     npolmo = number of polyatomic molecules (with more than 1 atom)
c     nmonat = number of monatomic molecules

      imolec = 0
      npolat = 0
      npolmo = 0
      do 300 iatom = 1, np
	if (molec(iatom) .eq. iatom) then
          imolec = imolec + 1
          molno(imolec) = iatom
          if (natmol(iatom) .gt. 1) then
             npolmo = npolmo + 1
             npolat = npolat + natmol(iatom)
          endif
        endif
 300  continue
      nmonat = np - npolat

      if (imolec .ne. nmolec) then
        write(isterr, *) 'clustr: miscounted molecules!'
        call ioclos()
        stop
      endif

c     check to see if the molecular composition has changed
c     (note that this will (intentionally) not detect reactions that
c     do not change molecular membership; its main use is to
c     determine whether the charge-neutral groups have changed)

      lrxn = .false.
      do 400 iatom = 1, np
         if (molec(iatom) .ne. omolec(iatom)) then
            lrxn = .true.
            go to 410
         endif
 400  continue

c     mark the status flag as fresh

      lmolol = .false.
      
 410  return
      end
c
c-----------------------------------------------------------------------
c     ggnprs builds a general distance-based pair list, using the 
c     supplied (squared) cutoff distance....sjs 12/29/03
c-----------------------------------------------------------------------
c     This could be expanded later to include type-dependent cutoffs
c     (once Voronoi volumes include spherical caps).
c-----------------------------------------------------------------------
c     
      subroutine ggnprs(np, cube, r0, lpbc, r2cut, nbrlst, ilist, jlist,
     .     nprs, niprs, r0p)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np         = number of particles
c     cube(d)    = box size in dimension d (Angstrom)
c     r0(a,d)    = position of atom a in dimension d (Angstrom)
c     lpbc(d)    = whether to use periodic boundary conditions in dimension d
c     r2cut      = square of cutoff (Angstrom^2)
c     nbrlist(a) = position on pairlist of first pair involving atom a
c     ilist(p)   = first (i) atom in pair # p
c     jlist(p)   = second (j) atom in pair # p
c     nprs       = maximum number of pairs 
c     niprs      = maximum number of neighbors per atom
c     r0p(a,d)   = position of atom a in dimension d when pairlist was built

      integer np
      dimension cube(ndim)
      dimension r0(npmax,ndim)
      logical lpbc(ndim)
      dimension nbrlst(npmax)
      dimension ilist(nprs)
      dimension jlist(nprs)
      dimension r0p(npmax,ndim)

      dimension ipr(niprs,npmax),
     .     npr(npmax)

      do 100 iatom = 1, np
         npr(iatom) = 0
 100  continue

      kpair = 0
      do 130 iatom = 1, np
         nbrlst(iatom) = kpair + 1

         rix = r0(iatom,1)
         riy = r0(iatom,2)
         riz = r0(iatom,3)

c     don't use a full N**2 loop; just use N**2/2 and build a forward-
c     looking list of pairs for future atoms to use.  should just do 
c     away with the double pairs altogether.

         do 120 jatom = iatom + 1, np
            dx = rix - r0(jatom,1)
            if (lpbc(1)) then
               dx = dx - cube(1) * anint(dx / cube(1))
            endif
            rsq = dx * dx
            if (rsq .gt. r2cut) then
               go to 120
            endif
            dx = riy - r0(jatom,2)
            if (lpbc(2)) then
               dx = dx - cube(2) * anint(dx / cube(2))
            endif
            rsq = rsq + dx * dx
            if (rsq .gt. r2cut) then
               go to 120
            endif
            dx = riz - r0(jatom,3)
            if (lpbc(3)) then
               dx = dx - cube(3) * anint(dx / cube(3))
            endif
            rsq = rsq + dx * dx
            if (rsq .gt. r2cut) then
               go to 120
            endif
            if (rsq .eq. 0.d0) then
               write(isterr, *) 
     .              'ggnprs: atoms ', iatom, ' and ', jatom, 
     .              ' have zero separation'
               call ioclos
               stop
            endif

            kpair = kpair + 1
            ilist(kpair) = iatom
            jlist(kpair) = jatom

c     remember this pair so that we don't have to recalculate it
c     the other way around.  this is somewhat kludgey and will
c     go away once the pair list is cut in half everywhere

            npr(jatom) = npr(jatom) + 1
            if (npr(jatom) .gt. niprs) then
               write(isterr, *) 'ggnprs: atom ', jatom, ' has ',
     .              npr(jatom), ' neighbors vs. limit of ', niprs
               write(isterr, *) 
     .              '    Increase niprs and recompile (unless this ',
     .              'is an unphysical result'
               call ioclos
               stop
            endif
            ipr(npr(jatom),jatom) = iatom

 120     continue

c     for the i > j pairs, don't recalculate distances but just
c     remember the pairs we found when looking at them the other
c     way around

         do 125 ipair = 1, npr(iatom)
            kpair = kpair + 1
            jatom = ipr(ipair,iatom)
            ilist(kpair) = iatom
            jlist(kpair) = jatom
 125     continue

 130  continue

c     so that the last atom knows how many neighbors he has:

      nbrlst(np+1) = kpair + 1

c     npairs must be strictly less than nlmax because sometimes there
c     is a fictitious pair in the last slot

      if (kpair .ge. nprs) then
         write(isterr, *) 'ggnprs: kpair value of ', kpair,
     .        ' exceeds nprs-1 value of ', nprs-1
         write(isterr, *) '        ensure that ggnprs is called with ',
     .        'an adequately large value of nprs'
         call ioclos
         stop
      endif

c     reset the last-pair-list-update positions

      do 210 idim = 1, 3
         do 200 iatom = 1, np
            r0p(iatom,idim) = r0(iatom,idim)
 200     continue
 210  continue

      return
      end
c

c-----------------------------------------------------------------------
c     gisold checks to see whether a pair list has expired, by 
c     determining whether any of the atoms have moved more than the
c     specified (squared) distance since the pair list was generated.
c     It is general enough for use on a variety of different pair lists.
c     ...sjs 12/29/03
c-----------------------------------------------------------------------
c
      logical function gisold(np, cube, r0, lpbc, r0p, r2trip)
      
      include 'common_files.inc'
      include 'common_pots.inc'

      integer np
      dimension cube(3)
      dimension r0(npmax,ndim)
      logical lpbc(ndim)
      dimension r0p(npmax,ndim)

      do 110 iatom = 1, np
         rsq = 0.d0
         do 100 idim = 1, ndim
            r = r0(iatom,idim) - r0p(iatom,idim)
            if (lpbc(idim)) then
               r = r - cube(idim) * anint(r / cube(idim))
            endif
            rsq = rsq + r * r
 100     continue
         if (rsq .gt. r2trip) then
            go to 200
         endif
 110  continue

c     neighbor list is still fresh

      gisold = .false.
      return

c     neighbor list needs to be redone

 200  continue

      gisold = .true.
      return

      end
c-----------------------------------------------------------------------
c     prpcov prepares the single-use pair lists and lookup tables for
c     covalent (REBO) interactions. These are consolidated here rather
c     than being called individually in calcforce and below, now that
c     qsolv also makes use of the covalent pair lists....sjs 8/1/07
c-----------------------------------------------------------------------
c     
      subroutine prpcov(np, cube, r0, lpbc)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np      = number of particles
c     cube(d) = box length in dimension d
c     r0(a,d) = position of atom a in dimension d
c     lpbc(d) = whether to use periodic boundary conditions in dimension d

      integer np
      real*8  cube(ndim)
      real*8  r0(npmax,ndim)
      logical lpbc(ndim)

c     local variables

      integer iatom, ijpair

c     check to see whether the pair list needs updating. only bother if
c     positions have moved since the last time we checked, and if we 
c     don't already know that it is expired

      if (lposol .and. .not. lcovol) then
         call chkprs(np, cube, r0, lpbc)
      endif

c     regenerate the covalent pair list if the old one has expired

      if (lcovol) then
         call genprs(np, cube, r0, lpbc)
      endif

c     fill up some variables for each pair, and also build a better
c     REBO-only, hydro/fluorocarbon-only pair list. These are only good
c     for this state, and need to be reconstructed every time the positions
c     change.

      nchch = 0
      n2chch = 0
      n2strt(1) = 1
      do 190 iatom = 2, np
         n2strt(iatom) = 0
 190  continue
     
      do 200 ijpair = 1, npairs
         call fillpr(ijpair, cube, r0, lpbc)
 200  continue

c     add this so that the last atom knows where to stop

      n2strt(np+1) = n2chch + 1
     
c     in case any atoms had no neighbors, go back and clean up
     
      do 205 iatom = 2, np
         if (n2strt(iatom) .eq. 0) then
            n2strt(iatom) = n2strt(iatom-1)
         endif
 205  continue

      return
      end
c
!-----------------------------------------------------------------------
!     bldstt sets up some state-specific variables and lookup tables.
!     it does NOT set up the pair lists, and should be called before
!     bldpls....sjs 3/16/12
!-----------------------------------------------------------------------
      subroutine bldstt(np, iatno, cube, lpbc, lewald, ewlkpl, ewkmxw)

      include 'common_files.inc'
      include 'common_pots.inc'

!     Arguments
!     np       = number of atoms
!     iatno(a) = atomic number of atom a
!     cube(d)  = box length in dimension d (A)
!     lpc(d)   = whether to use periodic boundary conditions in dimension d
!     lewald   = whether to perform Ewald sums
!     ewlkpl   = Ewald kappa * (min) box length
!     ewkmxw   = Ewald max |k|

      integer np
      integer iatno(npmax)
      real*8  cube(ndim)
      logical lpbc(ndim)
      logical lewald
      real*8  ewlkpl
      real*8  ewkmxw

!     local variables

      integer iatom, idim, itype
      logical lfinit

!     initialize some counters

      do 90 itype = 1, ntypes
         noa(itype) = 0
 90   continue

      nrba = 0
      nnra = 0

      do 100 iatom = 1, np

!     keep track of the internal atom type

         iat2ty(iatom) = kt(iatno(iatom))

!     complain if we have never heard of this atom
         
         if (iat2ty(iatom) .eq. 0) then
            write(isterr, *) 'bldstt: unknown atom type for atom ',
     .           iatom
            call ioclos
            stop
         endif

!     keep a count of how many atoms of each type we have

         noa(iat2ty(iatom)) = noa(iat2ty(iatom)) + 1

!     set up lists of REBO and non-REBO atoms

         if (lrebot(iat2ty(iatom))) then
            nrba = nrba + 1
            irblst(nrba) = iatom
         else
            nnra = nnra + 1
            inrlst(nnra) = iatom
         endif

 100  continue

!     set the box length to zero if not periodic, and calculate the
!     system volume if relevant

      lfinit = .true.
      vol = 1.d0
      do 300 idim = 1, ndim
         lfinit = lfinit .and. lpbc(idim)
         if (.not. lpbc(idim)) then
            cube(idim) = 0.d0
         endif
         vol = vol * cube(idim)
 300  continue
      if (.not. lfinit) then
         vol = bigpos
      endif

!     set everything that depends on box shape and volume

      call stvolp(cube, vol, lewald, ewlkpl, ewkmxw)

      return
      end
!
!-----------------------------------------------------------------------
!     build sets up the pair lists, state-specific variables, and 
!     state-specific lookup tables....sjs 3/16/10
!-----------------------------------------------------------------------
!
      subroutine build(np, r0, iatno, cube, lpbc, lewald, ewlkpl,
     .     ewkmxw, lvor)

      include 'common_files.inc'
      include 'common_pots.inc'

!     np       = number of atoms
!     r0(a,d)  = Cartesian coordinate of atom a in dimension d
!     iatno(a) = atomic number of atom a
!     cube(d)  = box length in dimension d (A)
!     lpbc(d)  = whether to use periodic boundary conditions in dimension d
!     lewald   = whether to perform Ewald sums
!     ewlkpl   = Ewald kappa * (min) box length
!     ewkmxw   = Ewald max |k|
!     lvor     = whether to calculate Voronoi volumes

      integer np
      real*8 r0(npmax, ndim)
      integer iatno(npmax)
      real*8 cube(ndim)
      logical lpbc(ndim)
      logical lewald
      real*8  ewlkpl
      real*8  ewkmxw
      logical lvor

!     set up the stuff that doesn't involve pair lists

      call bldstt(np, iatno, cube, lpbc, lewald, ewlkpl, ewkmxw)

!     and then set up the pair lists

      call bldpls(np, r0, cube, lpbc,
     .     lvor)

      return
      end
!
!-----------------------------------------------------------------------
!     bldpls sets up the pair lists by calculating them from the atom
!     positions....sjs 3/16/12
!-----------------------------------------------------------------------
!
      subroutine bldpls(np, r0, cube, lpbc,
     .     lvor)

      include 'common_files.inc'
      include 'common_pots.inc'

!     np       = number of atoms
!     r0(a,d)  = Cartesian coordinate of atom a in dimension d
!     cube(d)  = box length in dimension d (A)
!     lpc(d)   = whether to use periodic boundary conditions in dimension d
!     lvor     = whether to calculate Voronoi volumes

      integer np
      real*8  r0(npmax, ndim)
      real*8  cube(ndim)
      logical lpbc(ndim)
      logical lvor

c     local variables

      integer iatom, idim, itype
      logical lfinit
      real*8  vol

c     store the initial positions
      do 210 idim = 1, ndim
         do 200 iatom = 1, np
            r0l(iatom,idim) = r0(iatom,idim)
 200     continue
 210  continue
      lposol = .false.

      lcovol = .true.
      lvdwol = .true.

      if (.not. lvor) then
         call prpcov(np, cube, r0, lpbc)
         call gljprs(np, cube, r0, lpbc)
      endif

      return
      end
!
!-----------------------------------------------------------------------
!     pakpls repacks the pair list info sent from the KIM API into the
!     form needed. Must be called after bldstt....sjs 3/16/12
!-----------------------------------------------------------------------
!
      subroutine pakpls(knpair, kpairs, np, cube, r0, lpbc)

      include 'common_files.inc'
      include 'common_pots.inc'

!     Arguments
!     knpair      = number of KIM pairs
!     kpairs(p,n) = nth parter of KIM pair p
!     np          = number of atoms
!     cube(d) = box length in dimension d
!     r0(a,d) = position of atom a in dimension d
!     lpbc(d) = whether to use periodic boundary conditions in dimension d
      
      integer knpair
      integer kpairs(nlmax,2)
      integer np
      real*8  cube(ndim)
      real*8  r0(npmax,ndim)
      logical lpbc(ndim)

!     local variables

!     maximum number of REBO neighbors that an atom should realistically
!     be expected to have, using an inflated KIM pair cutoff

      integer prkim
      parameter (prkim = 100)

      integer ipr(prkim,npmax),
     .     npr(npmax), 
     .     iatom, iatom2, idim, ijbeg, ijend, ijpair,
     .     ilast, inra, ipair, irba, itype, iwatch,
     .     jkbeg, jkend, jkpair, jpair, jtype,
     .     katom, klbeg, klend, klpair, latom
      logical lpair(npmax)
      real*8  dx, rix, riy, riz, rr, rsq, rsqs

      do irba = 1, nrba
         npr(irblst(irba)) = 0
      end do

!     build a redundant pair list from the non-redundant pair list fed to us by
!     KIM

!     we assume here that the pair list is sorted by the first dimension.
!     i.e. pair (i,m) will always come before (j,n) when i < j

      npairs = 0
      ilast = -1
      do 120 ipair = 1, knpair
         iatom = kpairs(ipair,1)
         itype = iat2ty(iatom)

         if (lrebot(itype)) then

            rix = r0(iatom,1)
            riy = r0(iatom,2)
            riz = r0(iatom,3)

            if (iatom .ne. ilast) then
               if (ilast .ne. -1) then
                  ! we're on to a new iatom section of the list. add the 
                  ! redundant pairs
                  do jpair = 1, npr(iatom-1)
                     npairs = npairs + 1
                     jatom = ipr(jpair,iatom-1)
                     ihalf(npairs) = iatom-1
                     jhalf(npairs) = jatom
                  end do
               endif
               ilast = iatom
               nabors(iatom) = npairs + 1
            endif

            jatom = kpairs(ipair,2)

            jtype = iat2ty(jatom)

            if (lrebot(jtype)) then

               ! don't trust KIM's neighbor list. It may have used the LJ cutoff,
               ! and may have used a buffer (skin) width. Screen pairs to see if
               ! they are within REBO cutoffs plus buffer (which is set to zero
               ! here) before adding to the nabors-based REBO pair list

               dx = rix - r0(jatom,1)
               if (lpbc(1)) then
                  dx = rix - cube(1) * anint(dx / cube(1))
               endif
               rsq = dx * dx
               if (rsq .gt. r2rbmx(itype,jtype)) then
                  go to 120
               endif
               dx = riy - r0(jatom,2)
               if (lpbc(2)) then
                  dx = dx - cube(2) * anint(dx / cube(2))
               endif
               rsq = rsq + dx * dx
               if (rsq .gt. r2rbmx(itype,jtype)) then
                  go to 120
               endif
               dx = riz - r0(jatom,3)
               if (lpbc(3)) then
                  dx = dx - cube(3) * anint(dx / cube(3))
               endif
               rsq = rsq + dx * dx
               if (rsq .gt. r2rbmx(itype,jtype)) then
                  go to 120
               endif
               if (rsq .eq. 0.d0) then
                  write(isterr, *)
     .                 'genprs: atoms ', iatom, ' and ', jatom,
     .                 ' have zero separation'
                  call ioclos
                  stop
               endif

               ! this pair is within the REBO cutoff, add it to the REBO pairlist

               npairs = npairs + 1
               ihalf(npairs) = iatom
               jhalf(npairs) = jatom

               ! also remember this pair so that we don't have to recalculate it the
               ! other way around. this is somewhat kludgey and will go away once the
               ! pair list is cut in half everywhere

               npr(jatom) = npr(jatom) + 1
               if (npr(jatom) .gt. prkim) then
                  write(isterr, *) 'pakpls: atom ', jatom, ' has ',
     .                 npr(jatom),
     .                 'REBO neighbors on the KIM list vs. limit of ',
     .                 prkim
                  write(isterr, *)
     .                 '    you could increase prkim and recompile, ',
     .                 'but something else is probably wrong'
                  call ioclos
                  stop
               endif
               ipr(npr(jatom),jatom) = iatom
            endif
         endif
 120  continue
      ilast = iatom

      do iatom = ilast, np
         if (lrebot(iat2ty(iatom))) then
            nabors(iatom) = npairs + 1
            do jpair = 1, npr(iatom)
               npairs = npairs + 1
               jatom = ipr(jpair, iatom)
               ihalf(npairs) = iatom
               jhalf(npairs) = jatom
            end do
         endif
      end do

!     so that the last atom knows how many neighbors he has:

      nabors(np+1) = npairs + 1

!     this is a kludge to set the nabors pair list index for non-REBO atoms to
!     the same values as the nex following REBO atom (so that preceding REBO
!     atoms can find out how many neighbors they have). this can disappear once
!     all loops which refer to nabors are converted to 1,nrba loops

      do inra = 1, nnra
         iatom = inrlst(inra)
         do iatom2 = iatom, np
            if (lrebot(iat2ty(iatom2))) then
               nabors(iatom) = nabors(iatom2)
               go to 138
            endif
         end do
         nabors(iatom) = nabors(np + 1)
 138     continue
      end do

!     npairs must be strictly less than nlmax because sometimes there is a 
!     fictitious pair in the last slot

      if (npairs .ge. nlmax) then
         write(isterr, *) 'pakpls: npairs value of ', npairs,
     .        ' exceeds nlmax-1 value of ', nlmax-1
         write(isterr, *) '        increase nlmax and recompile'
         stop
      endif

!     fill up some variables for each pair, and also build a better, REBO-only,
!     hydro/fluorocarbon-only pair list. These are only good for this state,
!     and need to be reconstructed every time the positions change.

      nchch = 0
      n2chch = 0
      n2strt(1) = 1
      do iatom = 2, np
         n2strt(iatom) = 0
      end do

      do ijpair = 1, npairs
         call fillpr(ijpair, cube, r0, lpbc)
      end do

!     add this so that the last atom knows where to stop

      n2strt(np + 1) = n2chch + 1

!     in case any atoms had no neighbors, go back and clean up

      do iatom = 2, np
         if (n2strt(iatom) .eq. 0) then
            n2strt(iatom) = n2strt(iatom - 1)
         endif
      end do

!     use the nabors REBO pairlist to construct the LJ watchlist, consisting of
!     all pairs which are (1,2), (1,3), or (1,4) neighbors. 

!     should actually use the n2strt pairlist instead, since the nabors pairlist
!     built from KIM contains too many pairs

      nljwpr = 0
      do iatom = 1, np - 1
         do jatom = 1, np
            lpair(jatom) = .false.
         end do
         ijbeg = nabors(iatom)
         ijend = nabors(iatom+1) - 1
         do ijpair = ijbeg, ijend
            jatom = jhalf(ijpair)
            lpair(jatom) = .true.
            jkbeg = nabors(jatom)
            jkend = nabors(jatom+1) - 1
            do jkpair = jkbeg, jkend
               katom = jhalf(jkpair)
               lpair(katom) = .true.
               klbeg = nabors(katom)
               klend = nabors(katom+1) - 1
               do klpair = klbeg, klend
                  latom = jhalf(klpair)
                  lpair(latom) = .true.
               end do
            end do
         end do

         ! now that we know which atoms are neighbors of iatom, place them into
         ! the LJ watchlist in ascending order, and do not double-count pairs

         do jatom = iatom + 1, np
            if (lpair(jatom)) then
               nljwpr = nljwpr + 1
               if (nljwpr .gt. nwlmax) then
                  go to 9100
               endif
               iljwv(nljwpr) = iatom
               jljwv(nljwpr) = jatom
            endif
         end do

      end do

      ! now loop over all pairs, checking to see if they are within buffer range
      ! of the LJ potential. if so, double-check to make sure they aren't
      ! already on the watch list, and add them to the LJ safe pairlist.

      nljprs = 0
      iwatch = 1
      do iatom = 1, np - 1
         do 60 jatom = iatom + 1, np
            itype = iat2ty(iatom)
            jtype = iat2ty(jatom)

            ! skip out if there is no LJ interaction defined for this pair

            if (rslj(itype,jtype) .eq. 0.d0) then
               go to 60
            endif

            ! calculate the pair distance, and skip out if it is longer than the
            ! LJ cutoff plus the buffer distance (which is zero when using KIM)

            rsqs = 0.d0
            do idim = 1, ndim
               rr = r0(iatom,idim) - r0(jatom,idim)
               if (lpbc(idim)) then
                  rr = rr - cube(idim) * anint(rr/cube(idim))
               endif
               rsqs = rsqs + rr * rr
               if (rsqs .gt. rslj(itype,jtype)) then
                  go to 60
               endif
            end do

            ! we've found a valid LJ pair. now scan through the ordered watch
            ! watch list to see if we've alrady counted this pair

 51         continue
            if (iwatch .gt. nljwpr
     .           .or. iljwv(iwatch) .gt. iatom
     .           .or. (iljwv(iwatch) .eq. iatom
     .                 .and. jljwv(iwatch) .gt. jatom)) then
               go to 52
            else if (iljwv(iwatch) .eq. iatom
     .               .and. jljwv(iwatch) .eq. jatom) then
               go to 60
            else
               iwatch = iwatch + 1
               go to 51
            endif
 52         continue

            ! it wasn't on the watch list, so add it to the safe list

            nljprs = nljprs + 1
            if (nljprs .gt. nmabig) then
               go to 9000
            endif
            iv(nljprs) = iatom
            jv(nljprs) = jatom
 60      continue
      end do

      return

 9000 nsugg = int(dble(np) / iatom * nmabig)
      write(isterr, *) 'pakpls: just found safe LJ pair #', nljprs,
     .     ' vs. max of ', nmabig
      write(isterr, *) 'solution = increase nmabig to about ', nsugg,
     .     ' and recompile'
      go to 9500

 9100 nsugg = int(dble(np) / iatom * nljwpr)
      write(isterr, *) 'pakpls: just found watch LJ pair #', nljwpr,
     .     ' vs. max of ', nwlmax
      write(isterr, *) 'solution = increase nwlmax to about ',
     .     nsugg, ' and recompile'
      go to 9500

 9500 call ioclos
      stop
      
      end
!

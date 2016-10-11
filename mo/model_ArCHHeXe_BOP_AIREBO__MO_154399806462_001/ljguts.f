c-----------------------------------------------------------------------
c     ljguts
c-----------------------------------------------------------------------
c
      subroutine ljguts(np, igfunc, fint, cube, r0, uu, tote, lpbc)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np        = number of particles (input)
c     igfunc    = style of angular term (input)
c     fint(a,d) = internal force on atom a in dimension d (input/output)
c     cube(d)   = box length in dimension d (input)
c     r0(a,d)   = position of atom a in dimension d (input)
c     uu(d,d')  = d,d' element of internal virial tensor (input/output)
c     tote      = potential energy (input/output)
c     lpbc(d)   = whether to use periodic boundary conds in dimension d (input)

      integer np
      integer igfunc
      real*8  fint(npmax,ndim)
      real*8  cube(ndim)
      real*8  r0(npmax,ndim)
      real*8  uu(ndim,ndim)
      real*8  tote
      logical lpbc(ndim)

c     local variables

      logical lcntgt

      real*8 rrs(3)

c     reconstruct the LJ neighbor list if it is outdated

      if (lvdwol) then
         call gljprs(np, cube, r0, lpbc)
      endif

      pvdw=0.d0

c     first loop over the "safe" pairs, which are guaranteed not to
c     be first or second or third neighbors

      iprlst = iljsaf
      wwval = 0.d0

      do 900 ijpair = 1, nljprs
         iatom = iv(ijpair)
         jatom = jv(ijpair)
         itype = iat2ty(iatom)
         jtype = iat2ty(jatom)

c     calculate the pair distance and skip out if it is bigger than
c     the LJ cutoff

         do 100 idim = 1, ndim
            rrs(idim) = r0(iatom,idim) - r0(jatom,idim)
            if (lpbc(idim)) then
               rrs(idim) = rrs(idim) - cube(idim) 
     .              * anint(rrs(idim) / cube(idim))
            endif
 100     continue

         rsqs = rrs(1)**2+rrs(2)**2+rrs(3)**2
         if (rsqs .gt. r2mxlj(itype,jtype)) then
            go to 900
         endif

         rijsep = sqrt(rsqs)

         if (lljdir) then

c     calculate LJ interactions directly

            call ljpair(itype,jtype,rsqs,rijsep,vdw,dvdw)

         else

c     calculate LJ interactions from a table lookup with linear interpolation

            rt = rijsep / dellj
            ii = int(rt)+1

c        vdw is the LJ potential for this pair.  dvdw is (-1/r d/dr) of
c        the LJ potential.

            vdw = vljtb(ii,itype,jtype) +
     .           (vljtb(ii+1,itype,jtype) - vljtb(ii,itype,jtype)) *
     .           (rt - ii +1)
            dvdw = dvljtb(ii,itype,jtype) +
     .           (dvljtb(ii+1,itype,jtype) - dvljtb(ii,itype,jtype)) *
     .           (rt - ii + 1)
         endif

c     fsplin is the connectivity switching function.  C in my notes.
c     since these atoms are not first or second or third neighbors,
c     they are not connected and fsplin is 1 (leaving the LJ
c     interactions on)

         fsplin = 1.d0

c     rsplin is the cubic spline switching function that switches from
c     the contingent LJ interaction at low r to the absolute LJ
c     interaction at high r.  S(t_r(r)).
c     drspln is -1/r d/dr of rsplin

         lcntgt = .false.
         if (rijsep .gt. rlj1(itype,jtype)) then
            rsplin = 0.d0
            drspln = 0.d0
         else if (rijsep .gt. rlj0(itype,jtype)) then
            dr = rijsep - rlj0(itype,jtype)
            swidth = rlj1(itype,jtype) - rlj0(itype,jtype)
            tee = dr / swidth
            rsplin = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
            drspln = 6.d0 * tee * (1.d0 - tee) / rijsep / swidth
            lcntgt = .true.
         else
            rsplin = 1.d0
            drspln = 0.d0
            lcntgt = .true.
         endif

c     bsplin is the cubic spline switching function that switches from
c     full contingent LJ contributions at low bij to no contingent LJ
c     contributions at high bij.  S(t'(b)) in my notes.
c     all of the db/dr forces are taken care of inside spoof.

c     initialize bsplin only because NaN would be bad

         bsplin = 0.d0

         if (lcntgt) then

            bij = getbij(iprlst, ijpair, rijsep, rrs, wwval, rsplin,
     .           fsplin, igfunc, vdw, fint, cube, r0, uu, lpbc)

            if (bij .lt. bijmin(itype,jtype)) then
               bsplin = 1.d0

            else if (bij .lt. bijmax(itype,jtype)) then
               dbij = bij - bijmin(itype,jtype)
               swidth = bijmax(itype,jtype) - bijmin(itype,jtype)
               tee = dbij / swidth
               bsplin = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
            else
               bsplin = 0.d0
            endif
         endif

c     the LJ contribution to the energy is:
c       [1 - S(t(r))] C V^LJ + S(t(r)) S(t'(b)) C V^LJ
c          = [ 1 - S(t(r)) + S(t(r)) S'(t(b)) ] C V^LJ
c          = { S(t(r)) [ S(t'(b)) - 1 ] + 1 } C V^LJ

c     fsplin is guaranteed to be 1 here, so it's left out.

         pvdw = pvdw + (rsplin * (bsplin - 1.d0) + 1.d0) * vdw

c     From the LJ energy:
c       { S(t(r))   [      S(t'(b))       - 1 ] + 1 }   C     V^LJ
c     the d/dr derivative of the energy is:
c       dS/dt dt/dr [      S(t'(b))       - 1 ]         C     V^LJ   (1)
c     +   S(t(r))     dS/dt dt'/db db/dr                C     V^LJ   (2)
c     + { S(t(r))   [      S(t'(b))       - 1 ] + 1 } dC/dr   V^LJ   (3)
c     + { S(t(r))   [      S(t'(b))       - 1 ] + 1 }   C   dV^LJ/dr (4)

c     Line (2) is taken care of in spoof().  Line (3) is zero here
c     since dC/dr = 0.  And C = 1, which simplifies (1) and (4) a bit.

         dfac = drspln * (bsplin - 1.d0) * vdw
     .        + (rsplin * (bsplin - 1.d0) + 1.d0) * dvdw

         do 800 idim = 1, 3
            dr = dfac * rrs(idim)
            fint(iatom,idim) = fint(iatom,idim) + dr
            fint(jatom,idim) = fint(jatom,idim) - dr
            do 799 jdim = 1, 3
               uu(jdim,idim) = uu(jdim,idim) + dr * rrs(jdim)
 799        continue
 800     continue
 900  continue

c     next loop over the "watch list" -- the LJ pairs that are, or could
c     be, 1st, 2nd, or 3rd neighbors.  these don't get full LJ
c     treatment, and might have torsional interactions.

      iprlst = iljwat

      do 1900 ijpair = 1, nljwpr
         iatom = iljwv(ijpair)
         jatom = jljwv(ijpair)
         itype = iat2ty(iatom)
         jtype = iat2ty(jatom)

c     calculate the pair distance.  since this pair is in the watch
c     list, and thus inside twice the REBO buffered cutoff, it is
c     also inside the LJ cutoff, so we don't have to check.

c     if this pair is on the REBO pair list, then rcor and rij
c     are quicker ways to get to rijsep and rrs.  The lookup can be
c     made to be worth the effort, since a pass through the pair list
c     has to be made later anyway.

         do 1100 idim = 1, ndim
            rrs(idim) = r0(iatom,idim) - r0(jatom,idim)
            if (lpbc(idim)) then
               rrs(idim) = rrs(idim) - cube(idim)
     .              * anint(rrs(idim) / cube(idim))
            endif
 1100    continue

         rsqs = rrs(1)**2+rrs(2)**2+rrs(3)**2

         rijsep = sqrt(rsqs)

         if (lljdir) then

c     calculate LJ interactions directly

            call ljpair(itype,jtype,rsqs,rijsep,vdw,dvdw)

         else

c     calculate LJ interaction from a table lookup with linear interpolation

            rt = rijsep / dellj
            ii = int(rt)+1

            vdw = vljtb(ii,itype,jtype) +
     .           (vljtb(ii+1,itype,jtype) - vljtb(ii,itype,jtype)) *
     .           (rt - ii +1)
            dvdw = dvljtb(ii,itype,jtype) +
     .           (dvljtb(ii+1,itype,jtype) - dvljtb(ii,itype,jtype)) *
     .           (rt - ii + 1)
         endif

c     rsplin is the cubic spline switching function that switches from
c     the contingent LJ interaction at low r to the absolute LJ
c     interaction at high r.  S(t_r(r)).
c     drspln is -1/r d/dr of rsplin

         lcntgt = .true.
         if (rijsep .lt. rlj0(itype,jtype)) then
            rsplin = 1.d0
            drspln = 0.d0
         else if (rijsep .lt. rlj1(itype,jtype)) then
            dr = rijsep - rlj0(itype,jtype)
            swidth = rlj1(itype,jtype) - rlj0(itype,jtype)
            tee = dr / swidth
            rsplin = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
            drspln = 6.d0 * tee * (1.d0 - tee) / rijsep / swidth
         else
            rsplin = 0.d0
            drspln = 0.d0
            lcntgt = .false.
         endif

c     fsplin is the cubic spline switching function that switches from
c     full contingent LJ contribution at low f^c values (poor
c     connectivity) to no contingent LJ contribution at high f^c values
c     (full connectivity).  S(t_b(b*))

c     this must be calculated for all pairs now that it is applied to
c     non-contingent LJ interactions as well.

         fsplin = 1.d0

         ibond1 = 0
         ibond2 = 0
         ibond3 = 0
         ichain = 99
         wwval = 0.d0

c        loop over all REBO neighbors of the first atom in this pair

         ikbeg = n2strt(iatom)
         ikend = n2strt(iatom+1) - 1
         do 1410 ikprin = ikbeg, ikend
            ikpair = i2chch(ikprin)
            katom = jhalf(ikpair)

c           if the neighbor is the other half of the current LJ
c           pair, then assign it a spline.  also set the bond weight
c           factor, regardless of whether this pair is the strongest
c           connection

            if (katom .eq. jatom) then
               tfcspl = 1.d0 - fij(ikpair)
               wwval = fij(ikpair)

c           only keep this value of the spline if this
c           connection is stronger than any previous ones

c           direct connections win all ties

               if (tfcspl .lt. fsplin .or.
     .              (tfcspl .eq. fsplin .and. ichain .gt. 2)) then
                  fsplin = tfcspl
                  ibond1 = ikpair
                  ibond2 = 0
                  ibond3 = 0

c                 remember that our best connection so far was
c                 through a direct pair

                  ichain = 2
               endif

c     if we found jatom we know there's no better path that starts
c     through jatom, so don't bother doing the 2nd neighbor search

            else

c             also loop over all REBO neighbors of this
c             REBO neighbor

               klbeg = n2strt(katom)
               klend = n2strt(katom+1) - 1
               do 1400 klprin = klbeg, klend
                  klpair = i2chch(klprin)
                  latom = jhalf(klpair)

c                if we have looped back, give up on this path

                  if (latom .eq. iatom) then
                     go to 1400
                  endif

c                if the neighbor is the other half of the
c                current LJ pair, then assign it a spline

                  if (latom .eq. jatom) then
                     tfcspl = 1.d0 - fij(ikpair) * fij(klpair)

c                   only keep this value of the spline if
c                   this connection is stronger than any
c                   previous ones

c     in case of a tie, geminal beats vicinal

                     if (tfcspl .lt. fsplin .or.
     .                    (tfcspl .eq. fsplin .and. ichain .gt. 3)) then
                        fsplin = tfcspl
                        ibond1 = ikpair
                        ibond2 = klpair
                        ibond3 = 0

c                      remember that the best connection so far
c                      was through a single intermediary

                        ichain = 3
                     endif

c     if we found jatom we know it won't be a neighbor of itself so
c     don't bother doing the 3rd neighbor search

                  else

c     also loop over all 3rd neighbors

                     lmbeg = n2strt(latom)
                     lmend = n2strt(latom+1) - 1
                     do 1390 lmprin = lmbeg, lmend
                        lmpair = i2chch(lmprin)
                        matom = jhalf(lmpair)

c     if the neighbor is the other half of
c     the current LJ pair, then assign it a
c     spline

                        if (matom .eq. jatom) then
                           tfcspl = 1.d0 - fij(ikpair)
     .                          * fij(klpair) * fij(lmpair)

c     only keep this value of the spline
c     if this connection is stronger than
c     any previous ones

c     in case of a tie, this vicinal connection will not win

                           if (tfcspl .lt. fsplin) then
                              fsplin = tfcspl
                              ibond1 = ikpair
                              ibond2 = klpair
                              ibond3 = lmpair

c     remember that the best connection
c     so far was through two
c     intermediaries

                              ichain = 4
                           endif
                        endif
 1390                continue
                  endif
 1400          continue
            endif
 1410    continue

c     we still haven't figured out dC/dr for the connecting bonds, but
c     we'll save that until we know S(t'(b)).  we had to calculate
c     fsplin before calling getbij(), however.

c     if we're close enough to have to worry about contingent LJ
c     interactions, and any connecting bonds are partially dissociated,
c     we need S(t'(b)) = bsplin

         if (lcntgt .and. fsplin .ne. 0.d0) then

            bij = getbij(iprlst, ijpair, rijsep, rrs, wwval, rsplin,
     .           fsplin, igfunc, vdw, fint, cube, r0, uu, lpbc)

c     bsplin is the cubic spline switching function that switches from
c     full contingent LJ contributions at low bij to no contingent LJ
c     contributions at high bij.

            if (bij .lt. bijmin(itype,jtype)) then
               bsplin = 1.d0
            else if (bij .lt. bijmax(itype,jtype)) then
               dbij = bij - bijmin(itype,jtype)
               swidth = bijmax(itype,jtype) - bijmin(itype,jtype)
               tee = dbij / swidth
               bsplin = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
            else
               bsplin = 0.d0
            endif

         else


c           if the pair does not have a short enough separation to
c           need contingent LJ interactions, 
c           or is connected by 1, 2, or 3 solid bonds, 
c           then bsplin is irrelevant.
c           But we might as well make sure it's not INF or NaNQ,
c           because that would mess things up.

            bsplin = 0.d0

         endif

c        now that we know bsplin we can put in the dC/dr forces.
c        in particular, since the LJ energy looks like:
c          { S(t(r)) [ S(t'(b)) - 1 ] + 1 }   C   V^LJ
c        one of the derivative terms looks like:
c            S(t(r)) [ S(t'(b)) - 1 ] + 1 } dC/dr V^LJ

         if (fsplin .eq. 1.d0 .or. fsplin .eq. 0.d0) then

c           the pair is not connected. dC/dr = 0

         else

            if (ichain .eq. 2) then

c              the pair is directly connected.  C = 1-f(r_ij), so
c              dC/dr_ij = -df/dr|r_ij
c              ddrfac is -1/r_ij (S (S' - 1) + 1) dC/dr_ij V^LJ

               prefac = (rsplin * (bsplin - 1.d0) + 1.d0) * vdw
               ddrfac = prefac * dww(ibond1) / rijsep
               do 1770 idim = 1, 3
                  dr = ddrfac * rrs(idim)
                  fint(iatom,idim) = fint(iatom,idim) + dr
                  fint(jatom,idim) = fint(jatom,idim) - dr
                  do 1769 jdim = 1, 3
                     uu(jdim,idim) = uu(jdim,idim) + dr * rrs(jdim)
 1769             continue
 1770          continue
            else if (ichain .eq. 3) then

c              the pair is connected through an intermediate atom.
c              C = 1 - f(r_ik) f(r_kj), so we need
c              -df/dr|r_ik f(r_kj) and -f(r_ik) df/dr|r_kj.

               prefac = (rsplin * (bsplin - 1.d0) + 1.d0) * vdw
               katom = jhalf(ibond1)

c              first do the d/dr_ik derivatives.
c              ddrfac is -1/r_ik (S (S' - 1) + 1) dC/dr_ik V^LJ

               ddrfac = prefac * dww(ibond1) / rcor(ibond1)
     .              * fij(ibond2)
               do 1780 idim = 1, 3
                  dr = ddrfac * rijv(ibond1,idim)
                  fint(iatom,idim) = fint(iatom,idim) + dr
                  fint(katom,idim) = fint(katom,idim) - dr
                  do 1779 jdim = 1, 3
                     uu(jdim,idim) = uu(jdim,idim)
     .                    + dr * rijv(ibond1,jdim)
 1779             continue
 1780          continue

c              next do the d/dr_kj derivatives.
c              ddrfac is -1/r_kj (S (S' - 1) + 1) dC/dr_kj V^LJ

               ddrfac = prefac * fij(ibond1) * dww(ibond2)
     .              / rcor(ibond2)
               do 1790 idim = 1, 3
                  dr = ddrfac * rijv(ibond2,idim)
                  fint(katom,idim) = fint(katom,idim) + dr
                  fint(jatom,idim) = fint(jatom,idim) - dr
                  do 1789 jdim = 1, 3
                     uu(jdim,idim) = uu(jdim,idim)
     .                    + dr * rijv(ibond2,jdim)
 1789             continue
 1790          continue
            else

c              the pair is connected through two intermediate atoms.
c              C = 1 - f(r_ik) f(r_kl) f(r_lj), so we need
c              -df/dr|r_ik f(r_kl) f(r_lj) and
c              -f(r_ik) df/dr|r_kl f(r_lj) and
c              -f(r_ik) f(r_kl) df/dr|r_lj.

               prefac = (rsplin * (bsplin - 1.d0) + 1.d0) * vdw
               katom = jhalf(ibond1)
               latom = jhalf(ibond2)

c              first do the d/dr_ik derivatives
c              ddrfac is -1/r_ik (S (S' - 1) + 1) dC/dr_ik V^LJ

               ddrfac = prefac * dww(ibond1) / rcor(ibond1)
     .              * fij(ibond2) * fij(ibond3)
               do 1792 idim = 1, 3
                  dr = ddrfac * rijv(ibond1,idim)
                  fint(iatom,idim) = fint(iatom,idim) + dr
                  fint(katom,idim) = fint(katom,idim) - dr
                  do 1791 jdim = 1,3
                     uu(jdim,idim) = uu(jdim,idim)
     .                    + dr * rijv(ibond1,jdim)
 1791             continue
 1792          continue

c              next do the d/dr_kl derivatives.
c              ddrfac is -1/r_kl (S (S' - 1) + 1) dC/dr_kl V^LJ

               ddrfac = prefac * fij(ibond1) * dww(ibond2)
     .              / rcor(ibond2) * fij(ibond3)
               do 1794 idim = 1, 3
                  dr = ddrfac * rijv(ibond2,idim)
                  fint(katom,idim) = fint(katom,idim) + dr
                  fint(latom,idim) = fint(latom,idim) - dr
                  do 1793 jdim = 1, 3
                     uu(jdim,idim) = uu(jdim,idim)
     .                    + dr * rijv(ibond2,jdim)
 1793             continue
 1794          continue

c              next do the d/dr_lj derivatives.
c              ddrfac is -1/r_lj (S (S' - 1) + 1) dC/dr_lj V^LJ

               ddrfac = prefac * fij(ibond1) * fij(ibond2)
     .              * dww(ibond3) / rcor(ibond3)
               do 1796 idim = 1, 3
                  dr = ddrfac * rijv(ibond3,idim)
                  fint(latom,idim) = fint(latom,idim) + dr
                  fint(jatom,idim) = fint(jatom,idim) - dr
                  do 1795 jdim = 1, 3
                     uu(jdim,idim) = uu(jdim,idim)
     .                    + dr * rijv(ibond3,jdim)
 1795             continue
 1796          continue
            endif
         endif

c     the LJ contribution to the energy is:
c       [1 - S(t(r))] C V^LJ + S(t(r)) S(t'(b)) C V^LJ
c          = [ 1 - S(t(r)) + S(t(r)) S(t'(b)) ] C V^LJ
c          = { S(t(r)) [ S(t'(b)) - 1 ] + 1 } C V^LJ

         pvdw = pvdw + (rsplin * (bsplin - 1.d0) + 1.d0) * fsplin * vdw

c     From the LJ energy:
c       { S(t(r))   [      S(t'(b))      - 1 ] + 1 }   C     V^LJ
c     the d/dr derivative of the energy is:
c       dS/dt dt/dr [      S(t'(b))      - 1 ]     }   C     V^LJ   (1)
c     +   S(t(r))     dS/dt dt'/db db/dr               C     V^LJ   (2)
c     + { S(t(r))   [      S(t'(b))      - 1 ] + 1 } dC/dr   V^LJ   (3)
c     + { S(t(r))   [      S(t'(b))      - 1 ] + 1 }   C   dV^LJ/dr (4)

c     Line (2) is taken care of in spoof().  Line (3) is taken care of
c     above.

         dfac = drspln * (bsplin - 1.d0) * fsplin * vdw +
     .        (rsplin * (bsplin - 1.d0) + 1.d0) * fsplin * dvdw

         do 1800 idim = 1, 3
            dr = dfac * rrs(idim)
            fint(iatom,idim) = fint(iatom,idim) + dr
            fint(jatom,idim) = fint(jatom,idim) - dr
            do 1799 jdim = 1, 3
               uu(jdim,idim) = uu(jdim,idim) + dr * rrs(jdim)
 1799       continue
 1800    continue
 1900 continue

      tote = tote + pvdw

      return
      end

c-----------------------------------------------------------------------
c     ljpair calculates the pairwise interaction of LJ part of the
c     airebo potential, if system is small enough  - MF 01/04/2011
c-----------------------------------------------------------------------

      subroutine ljpair(itype, jtype, rsqs, rijsep, vdw, dvdw)
     
      include 'common_files.inc'
      include 'common_pots.inc'

c     itype  = atom type of i (input)
c     jtype  = atom type of j (input)
c     rsqs   = squared distance between atom i, atom j (input)
c     rijsep = distance between atoms i and j (input)
c     vdw    = switched LJ interaction (output)
c     dvdw   = -1/r d/dr of the switched LJ interaction (output)

      integer itype
      integer jtype
      real*8 rsqs
      real*8 vdw
      real*8 dvdw

c     local variables

c     r6 = r^6 term in LJ calc
c     sigwid = hardwired
c     slj = switch part of interaction
c     vlj = vlj interaction, no switch

      real*8 drij, dslj, dvlj, r6,
     .     slj, swidth, tee, vlj
  
c     vlj is the LJ potential
c     dvlj is (-1/r d/dr) of the LJ potential

      r6 = (sig(itype,jtype) / rsqs) ** 3
      vlj = eps(itype,jtype) * r6 * (r6 - 1.d0) 

      dvlj = eps(itype,jtype) / rsqs * r6 * (12.d0 * r6 - 6.d0)

c     slj is the LJ switching function at long distances
c     dslj is (-1/r d/dr) of the LJ switching function

      if (rijsep .gt. rmaxlj(itype,jtype)) then
         slj = 0.d0
         dslj = 0.d0
      else if (rijsep .gt. rminlj(itype,jtype)) then
         drij = rijsep - rminlj(itype,jtype)
         swidth = rmaxlj(itype,jtype) - rminlj(itype,jtype)
         tee = drij / swidth
         slj = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
         dslj = 6.d0 * tee * (1.d0 - tee) / rijsep / swidth
      else
         slj = 1.d0
         dslj = 0.d0
      endif

c     vdw is the LJ pairwise interaction between the atoms (before
c     adaptive switches)
c     dvdw is (-1/r) d/dr of the LJ interaction

      vdw = vlj * slj
      dvdw = dvlj * slj + vlj * dslj

      return
      end

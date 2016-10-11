c-----------------------------------------------------------------------
c     dihedrals calculates the energy and forces associated with the
c     torsional potential for dihedral angles.  (There are additional
c     terms for dihedral rotations around double bonds.  These are
c     accounted for in the bond order, in the pibond routine.)
c     ...sjs 7/29/98
c-----------------------------------------------------------------------
c     It would possibly be better to include these interactions as a
c     penalty to the bond order, rather than as a purely energetic
c     potential component.
c-----------------------------------------------------------------------
c     Much of this work is duplicated in pibond.  These calculations
c     should be moved into pibond if only to cut computation time.
c-----------------------------------------------------------------------
c
      subroutine dihed(ijpair, fint, uu, tote)

      include 'common_files.inc'
      include 'common_pots.inc'

      dimension fint(npmax, 3)
      dimension uu(ndim, ndim)
      real*8 tote

      parameter (ncc = 20)

      dimension iklist(ncc), jllist(ncc)

      dimension cosjik(ncc), cosijl(ncc), sinjik(ncc), sinijl(ncc),
     .     tspjik(ncc), dtsjik(ncc), tspijl(ncc), dtsijl(ncc),
     .     cjik(3), cijl(3),
     .     dndij(3), dndik(3), dndjl(3),
     .     rijvec(3), rikvec(3), rilvec(3), rjkvec(3), rjlvec(3),
     .     tmpvec(3)

      iatom = ihalf(ijpair)
      jatom = jhalf(ijpair)
      rijvec(1) = rijv(ijpair,1)
      rijvec(2) = rijv(ijpair,2)
      rijvec(3) = rijv(ijpair,3)

c     k            we need to loop over all REBO neighbors k of i
c      \           and all REBO neighbors l of j to find the
c       i--j       4-atom groups which define the dihedral angles
c           \
c            l

c     r_ij bond

      rij = rcor(ijpair)
      rij2 = rij ** 2

c     do a preliminary loop over k and l separately to build some
c     lists

      nk = 0
      ikbeg = n2strt(iatom)
      ikend = n2strt(iatom+1) - 1
      do 100 ik = ikbeg, ikend
         ikpair = i2chch(ik)
         katom = jhalf(ikpair)
         if (katom .eq. jatom) then
            go to 100
         endif

c        r_ik

         rik = rcor(ikpair)
         rik2 = rik ** 2

c        r_jk

         rjk2 = (rijv(ikpair,1) - rijv(ijpair,1)) ** 2
     .        + (rijv(ikpair,2) - rijv(ijpair,2)) ** 2
     .        + (rijv(ikpair,3) - rijv(ijpair,3)) ** 2

c        we know 3 sides, calculate angle

         costmp = 0.5d0 * (rij2 + rik2 - rjk2) / rij / rik

c        seemingly impossible, but happens due to roundoff:

         if (costmp .lt. -1.d0) then
            costmp = -1.d0
         else if (costmp .gt. 1.d0) then
            costmp = 1.d0
         endif
         sintmp = sqrt(1.d0 - costmp ** 2)

c     tsplin is the cubic spline switching function that turns off the
c     torsion interactions at nearly linear bond angles.
c     S(t_theta(cos(theta))).  dtspln is -dS/d(cos theta).

         if (costmp .gt. tthmax) then
            tsplin = 0.d0
            dtspln = 0.d0
         else if (costmp .gt. tthmin) then
            dtthet = costmp - tthmin
            twidth = tthmax - tthmin
            tee = dtthet / twidth
            tsplin = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
            dtspln = 6.d0 * tee * (1.d0 - tee) / twidth
         else
            tsplin = 1.d0
            dtspln = 0.d0
         endif
         nk = nk + 1
         iklist(nk) = ikpair
         cosjik(nk) = costmp
         sinjik(nk) = sintmp
         tspjik(nk) = tsplin
         dtsjik(nk) = dtspln
 100  continue

      nl = 0
      jlbeg = n2strt(jatom)
      jlend = n2strt(jatom+1) - 1
      do 200 jl = jlbeg, jlend
         jlpair = i2chch(jl)
         latom = jhalf(jlpair)
         if (latom .eq. iatom) then
            go to 200
         endif

c        r_jl

         rjl = rcor(jlpair)
         rjl2 = rjl ** 2

c        r_il

         ril2 = (rijv(jlpair,1) + rijv(ijpair,1)) ** 2
     .        + (rijv(jlpair,2) + rijv(ijpair,2)) ** 2
     .        + (rijv(jlpair,3) + rijv(ijpair,3)) ** 2

c        we know 3 sides, calculate angle

         costmp = 0.5d0 * (rij2 + rjl2 - ril2) / rij / rjl

c        seemingly impossible, but happens due to roundoff error:

         if (costmp .lt. -1.0d0) then
            costmp = -1.0d0
         else if (costmp .gt. 1.d0) then
            costmp = 1.d0
         endif
         sintmp = sqrt(1.d0 - costmp ** 2)

c     tsplin is the cubic spline switching function that turns off the
c     torsion interactions at nearly linear bond angles.
c     S(t_theta(cos(theta))).  dtspln is -dS/d(cos theta).

         if (costmp .gt. tthmax) then
            tsplin = 0.d0
            dtspln = 0.d0
         else if (costmp .gt. tthmin) then
            dtthet = costmp - tthmin
            twidth = tthmax - tthmin
            tee = dtthet / twidth
            tsplin = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
            dtspln = 6.d0 * tee * (1.d0 - tee) / twidth
         else
            tsplin = 1.d0
            dtspln = 0.d0
         endif
         nl = nl + 1
         jllist(nl) = jlpair
         cosijl(nl) = costmp
         sinijl(nl) = sintmp
         tspijl(nl) = tsplin
         dtsijl(nl) = dtspln
 200  continue

      do 1010 ik = 1, nk
         ikpair = iklist(ik)
         katom = jhalf(ikpair)
         ktype = iat2ty(katom)
         rikvec(1) = rijv(ikpair,1)
         rikvec(2) = rijv(ikpair,2)
         rikvec(3) = rijv(ikpair,3)
         rik = rcor(ikpair)
         rik2 = rik ** 2
         rjkvec(1) = rikvec(1) - rijvec(1)
         rjkvec(2) = rikvec(2) - rijvec(2)
         rjkvec(3) = rikvec(3) - rijvec(3)
         rjk2 = rjkvec(1) ** 2 + rjkvec(2) ** 2 + rjkvec(3) ** 2
         rjk = sqrt(rjk2)
         do 1000 jl = 1, nl

c     if splines are 1 and spline derivatives are 0
c     (eg. for 2 collinear bonds), then all forces and 
c     energies will be zero and we can skip out

c     sometimes dS/dt > 0 because of roundoff error, so be tolerant

            if ((tspjik(ik) .eq. 1.d0
     .           .or. (tspijl(jl) .eq. 1.d0
     .                 .and. dtsijl(jl) .le. dtttol))
     .          .and. (tspijl(jl) .eq. 1.d0
     .                  .or. (tspjik(ik) .eq. 1.d0 
     .                        .and. dtsjik(ik) .le. dtttol))) then
               go to 1000
            endif

            jlpair = jllist(jl)
            latom = jhalf(jlpair)
            ltype = iat2ty(latom)

c     ignore cyclic trimers:

            if (latom .eq. katom) then
               go to 1000
            endif

            rjlvec(1) = rijv(jlpair,1)
            rjlvec(2) = rijv(jlpair,2)
            rjlvec(3) = rijv(jlpair,3)
            rjl = rcor(jlpair)
            rjl2 = rjl ** 2
            rilvec(1) = rjlvec(1) + rijvec(1)
            rilvec(2) = rjlvec(2) + rijvec(2)
            rilvec(3) = rjlvec(3) + rijvec(3)
            ril2 = rilvec(1) ** 2 + rilvec(2) ** 2 + rilvec(3) ** 2
            ril = sqrt(ril2)

c           cross products:
c           cjik = r_ji x r_ik = r_ik x r_ij
c           cijl = r_ij x r_jl

            call cross(cjik, rikvec, rijvec)
            call cross(cijl, rijvec, rjlvec)

c           magnitudes of the cross products:
c           cjikm = |r_ji x r_ik| = r_ij r_ik sin(theta_jik)
c           cijlm = |r_ij x r_jl| = r_ij r_jl sin(theta_ijl)

            cjikm = rcor(ijpair) * rcor(ikpair) * sinjik(ik)
            cijlm = rcor(ijpair) * rcor(jlpair) * sinijl(jl)

c           cwnum = (r_ji x r_ik) . (r_ij x r_jl)
c           cwdnom = |r_ji x r_ik| |r_ij x r_jl|

            cwnum = dot(cjik, cijl, ndim)
            cwdnom = cjikm * cijlm

c           the torsion angle is undefined for 2 collinear bonds.  If the
c           switching function ranges are set appropriately (i.e. if 
c           tthmin = -1 or larger), then the collinear bonds will have
c           skipped over this section of code.  But check just in case.

            if (cwdnom .eq. 0.d0) then
               write(isterr, *) 'dihed: dihedral angle is undefined ',
     .              'for collinear bonds.  Make sure tthmin >= -1.0'
               call ioclos
               stop
            endif

c           omega = angle between e_jik and e_ijl (NOT dihedral angle)
c           phi = dihedral angle
c           omega = pi + phi (or pi - phi, depending on a handedness choice)
c           cos(omega) = cos(pi +/- phi) = -cos(phi)
c                      = num / denom = e_jik . e_ijl
c           e_jik = (r_ji x r_ik) / |r_ji x r_ik|
c           e_ijl = (r_ij x r_jl) / |r_ij x r_jl|

            cw = cwnum / cwdnom

c           The torsional potential for dihedral rotations have the
c           form
c              V(phi) = tor0 + tor10 * cos^10(phi/2)
c           It's easier to calculate omega = phi + pi than to
c           calculate phi, so
c              V(phi) = tor0 + tor10 * (cos^2(phi/2))^5
c                     = tor0 + tor10 * (1/2 (1 + cos(phi)))^5
c                     = tor0 + tor10 * (1/2 (1 - cos(omega)))^5
c           When added to the potential, the torsion term is
c           scaled by the 3 component bond weights
c              V = V(phi) f_ik f_ij f_jl

            phipc = tor0(ktype,ltype) + tor10(ktype,ltype)
     .           * (0.5d0 * (1.d0 - cw)) ** 5
            tote = tote + phipc
     .           * fij(ikpair) * fij(ijpair) * fij(jlpair)
     .           * (1.d0 - tspjik(ik)) * (1.d0 - tspijl(jl))

c           Take care of some easy derivative terms:
c              V(phi) d(f_ik)/d(r_ik) f_ij f_jl (1-S_jik) (1-S_ijl)
c              V(phi) f_ik d(f_ij)/d(r_ij) f_jl (1-S_jik) (1-S_ijl)
c              V(phi) f_ik f_ij d(f_jl)/d(r_jl) (1-S_jik) (1-S_ijl)

            fcikpc = -phipc * dww(ikpair) * fij(ijpair) * fij(jlpair)
     .           / rcor(ikpair) * (1.d0 - tspjik(ik))
     .           * (1.d0 - tspijl(jl))
            do 300 idim = 1, ndim
               fcpc = fcikpc * rijv(ikpair,idim)
               fint(iatom,idim) = fint(iatom,idim) + fcpc
               fint(katom,idim) = fint(katom,idim) - fcpc
               do 295 jdim = 1, ndim
                  uu(jdim,idim) = uu(jdim,idim) + fcpc * rikvec(jdim)
 295           continue
 300        continue

            fcijpc = -phipc * fij(ikpair) * dww(ijpair) * fij(jlpair)
     .           / rcor(ijpair) * (1.d0 - tspjik(ik))
     .           * (1.d0 - tspijl(jl))
            do 310 idim = 1, ndim
               fcpc = fcijpc * rijv(ijpair,idim)
               fint(iatom,idim) = fint(iatom,idim) + fcpc
               fint(jatom,idim) = fint(jatom,idim) - fcpc
               do 305 jdim = 1, ndim
                  uu(jdim,idim) = uu(jdim,idim) + fcpc * rijvec(jdim)
 305           continue
 310        continue

            fcjlpc = -phipc * fij(ikpair) * fij(ijpair) * dww(jlpair)
     .           / rcor(jlpair) * (1.d0 - tspjik(ik))
     .           * (1.d0 - tspijl(jl))
            do 320 idim = 1, ndim
               fcpc = fcjlpc * rijv(jlpair,idim)
               fint(jatom,idim) = fint(jatom,idim) + fcpc
               fint(latom,idim) = fint(latom,idim) - fcpc
               do 315 jdim = 1, ndim
                  uu(jdim,idim) = uu(jdim,idim) + fcpc * rjlvec(jdim)
 315           continue
 320        continue

c           The dV(phi)/dr derivatives are tougher.
c           d(cos omega)/d(r_ij)
c              = d(num)/d(r_ij) 1/denom
c                 - num / denom^2 d(denom)/d(r_ij)

c           cos omega = num / denom

c           num = (rji x rik) . (rij x rjl)
c           denom = |rji x rik| . |rij x rjl|
c                 = rij rik sin(theta_jik) rij rjl sin(theta_ijl)

c           d(num)/d(x_ij) = - (rik x (rij x rjl))_x
c                               - ((rji x rik) x rjl)_x
c                  = ((rij x rjl) x rik)_x + (rjl x (rji x rik))_x
c           d(num)/d(x_ik) = (rij x (rij x rjl))_x
c           d(num)/d(x_jl) = ((rji x rik) x rij)_x

            call cross(dndij, cijl, rikvec)
            call cross(tmpvec, rjlvec, cjik)
            dndij(1) = dndij(1) + tmpvec(1)
            dndij(2) = dndij(2) + tmpvec(2)
            dndij(3) = dndij(3) + tmpvec(3)
            call cross(dndik, rijvec, cijl)
            call cross(dndjl, cjik, rijvec)

c           Using law of cosines,
c           d(cos(theta_jik))/d(rij) =
c              (rij^2 - rik^2 + rjk^2) / (2 rij^2 rik)
c           d(cos(theta_jik))/d(rik) =
c              (-rij^2 + rik^2 + rjk^2) / (2 rij rik^2)
c           d(cos(theta_jik))/d(rjk) = - rjk / (rij rik)

c           d(cos(theta_ijl))/d(rij) =
c              (rij^2 - rjl^2 + ril^2) / (2 rij^2 rjl)
c           d(cos(theta_ijl))/d(rjl) =
c              (-rij^2 + rjl^2 + ril^2) / (2 rij rjl^2)
c           d(cos(theta_ijl))/d(ril) = - ril / (rij rjl)

            dcidij = (rij2 - rik2 + rjk2) / (2.d0 * rij2 * rik)
            dcidik = (rik2 - rij2 + rjk2) / (2.d0 * rij * rik2)
            dcidjk = -rjk / (rij * rik)

            dcjdji = (rij2 - rjl2 + ril2) / (2.d0 * rij2 * rjl)
            dcjdjl = (rjl2 - rij2 + ril2) / (2.d0 * rij * rjl2)
            dcjdil = -ril / (rij * rjl)

c           d(sin(theta))/d(r) = - cos / sin d(cos)/d(r)

            dsidij = -cosjik(ik) / sinjik(ik) * dcidij
            dsidik = -cosjik(ik) / sinjik(ik) * dcidik
            dsidjk = -cosjik(ik) / sinjik(ik) * dcidjk

            dsjdji = -cosijl(jl) / sinijl(jl) * dcjdji
            dsjdjl = -cosijl(jl) / sinijl(jl) * dcjdjl
            dsjdil = -cosijl(jl) / sinijl(jl) * dcjdil

c           d|rji x rik|/d(rij) = rik sin(theta_jik)
c                                    + rij rik d(sin)/d(rij)
c           d|rji x rik|/d(rik) = rij sin(theta_jik)
c                                    + rij rik d(sin)/d(rik)
c           d|rji x rik|/d(rjk) = rij rik d(sin)/d(rjk)

c           likewise for the |rij x rjl| derivatives

            dxidij = rik * sinjik(ik) + rij * rik * dsidij
            dxidik = rij * sinjik(ik) + rij * rik * dsidik
            dxidjk =                    rij * rik * dsidjk

            dxjdji = rjl * sinijl(jl) + rij * rjl * dsjdji
            dxjdjl = rij * sinijl(jl) + rij * rjl * dsjdjl
            dxjdil =                    rij * rjl * dsjdil

c           d(denom)/d(r) = d|rji x rik|/d(r) |rij x rjl|
c                             + |rji x rik| d|rij x rjl|/d(r)

            ddndij = dxidij * cijlm + cjikm * dxjdji
            ddndik = dxidik * cijlm
            ddndjk = dxidjk * cijlm
            ddndjl =                  cjikm * dxjdjl
            ddndil =                  cjikm * dxjdil

c           cos omega = num / denom
c           d(cos omega)/d(denom) = - num / denom^2

            dcwddn = -cwnum / cwdnom ** 2

c           d(cos omega)/d(num) = 1/denom

            dcwdn = 1.d0 / cwdnom

c           V(phi) = tor0 + tor10 * (1/2 (1 - cos(omega)))^5
c                  = tor0 + 1/32 tor10 (1 - cos(omega))^5
c           dV(phi)/d(cos omega) = -5/32 tor10 (1 - cos(omega))^4

            dvpdcw = -5.d0 / 32.d0 * tor10(ktype,ltype)
     .           * (1.d0 - cw) ** 4

c           V = V(phi) fij fik fjl (1-S_jik) (1-S_ijl)
c           dV/d(V(phi)) = fij fik fjl (1-S_jik) (1-S_ijl)

            dvdv = fij(ijpair) * fij(ikpair) * fij(jlpair)
     .           * (1.d0 - tspjik(ik)) * (1.d0 - tspijl(jl))

c           forces are a little tricky because we have d(num)/d(xij)
c           and d(denom)/d(rij).

c           d(V(phi))/d(xij) = dV/d(cos omega) d(cos omega)/d(xij)
c              = dV(phi)/d(cos) ( d(cos)/d(num) d(num)/d(xij)
c                   + d(cos)/d(denom) d(denom)/d(rij) d(rij)/d(xij)

c           Fij ~= -dV/d(V(phi)) d(V(phi))/d(xij)
c           (mass factor is taken care of later)

            fcpc = -dvdv * dvpdcw
            do 400 idim = 1, ndim

               fcijpc = fcpc * (dcwdn * dndij(idim)
     .              + dcwddn * ddndij * rijvec(idim) / rij)
               fint(iatom,idim) = fint(iatom,idim) + fcijpc
               fint(jatom,idim) = fint(jatom,idim) - fcijpc

               do 395 jdim = 1, ndim
                  uu(jdim,idim) = uu(jdim,idim) + fcijpc * rijvec(jdim)
 395           continue
               fcikpc = fcpc * (dcwdn * dndik(idim)
     .              + dcwddn * ddndik * rikvec(idim) / rik)
               fint(iatom,idim) = fint(iatom,idim) + fcikpc
               fint(katom,idim) = fint(katom,idim) - fcikpc

               do 396 jdim = 1, ndim
                  uu(jdim,idim) = uu(jdim,idim) + fcikpc * rikvec(jdim)
 396           continue
               fcjkpc = fcpc * dcwddn * ddndjk * rjkvec(idim) / rjk
               fint(jatom,idim) = fint(jatom,idim) + fcjkpc
               fint(katom,idim) = fint(katom,idim) - fcjkpc
 
               do 397 jdim = 1, ndim
                  uu(jdim,idim) = uu(jdim,idim) + fcjkpc * rjkvec(jdim)
 397           continue
               fcjlpc = fcpc * (dcwdn * dndjl(idim)
     .              + dcwddn * ddndjl * rjlvec(idim) / rjl)
               fint(jatom,idim) = fint(jatom,idim) + fcjlpc
               fint(latom,idim) = fint(latom,idim) - fcjlpc

               do 398 jdim = 1, ndim
                  uu(jdim,idim) = uu(jdim,idim) + fcjlpc * rjlvec(jdim)
 398           continue
               fcilpc = fcpc * dcwddn * ddndil * rilvec(idim) / ril
               fint(iatom,idim) = fint(iatom,idim) + fcilpc
               fint(latom,idim) = fint(latom,idim) - fcilpc

               do 399 jdim = 1, ndim
                  uu(jdim,idim) = uu(jdim,idim) + fcilpc * rilvec(jdim)
 399           continue
 400        continue

c     the only thing left is the d(S(t(cos(theta))))/dx derivatives.

c     V^phi f_ik f_ij f_jl (-dS/dt) dt/d(cos theta_jik)
c            x d(cos theta_jik)/dr (1-S_ijl)
c     V^phi f_ik f_ij f_jl (1-S_jik) (-dS/dt) dt/d(cos theta_ijl)
c            x d(cos theta_ijl)/dr

            fcpc = -phipc * fij(ikpair) * fij(ijpair) * fij(jlpair)
     .           * dtsjik(ik) * (1.d0 - tspijl(jl))
            do 540 idim = 1, 3
               fcijpc = fcpc * dcidij * rijvec(idim) / rij
               fint(iatom,idim) = fint(iatom,idim) + fcijpc
               fint(jatom,idim) = fint(jatom,idim) - fcijpc
               do 500 jdim = 1, 3
                  uu(jdim,idim) = uu(jdim,idim) + fcijpc * rijvec(jdim)
 500           continue
               fcikpc = fcpc * dcidik * rikvec(idim) / rik
               fint(iatom,idim) = fint(iatom,idim) + fcikpc
               fint(katom,idim) = fint(katom,idim) - fcikpc
               do 510 jdim = 1, 3
                  uu(jdim,idim) = uu(jdim,idim) + fcikpc * rikvec(jdim)
 510           continue
               fcjkpc = fcpc * dcidjk * rjkvec(idim) / rjk
               fint(jatom,idim) = fint(jatom,idim) + fcjkpc
               fint(katom,idim) = fint(katom,idim) - fcjkpc
               do 520 jdim = 1, 3
                  uu(jdim,idim) = uu(jdim,idim) + fcjkpc * rjkvec(jdim)
 520           continue
 540        continue

            fcpc = -phipc * fij(ikpair) * fij(ijpair) * fij(jlpair)
     .           * (1.d0 - tspjik(ik)) * dtsijl(jl)
            do 580 idim = 1, 3
               fcijpc = fcpc * dcjdji * rijvec(idim) / rij
               fint(iatom,idim) = fint(iatom,idim) + fcijpc
               fint(jatom,idim) = fint(jatom,idim) - fcijpc
               do 550 jdim = 1, 3
                  uu(jdim,idim) = uu(jdim,idim) + fcijpc * rijvec(jdim)
 550           continue
               fcjlpc = fcpc * dcjdjl * rjlvec(idim) / rjl
               fint(jatom,idim) = fint(jatom,idim) + fcjlpc
               fint(latom,idim) = fint(latom,idim) - fcjlpc
               do 560 jdim = 1, 3
                  uu(jdim,idim) = uu(jdim,idim) + fcjlpc * rjlvec(jdim)
 560           continue
               fcilpc = fcpc * dcjdil * rilvec(idim) / ril
               fint(iatom,idim) = fint(iatom,idim) + fcilpc
               fint(latom,idim) = fint(latom,idim) - fcilpc
               do 570 jdim = 1, 3
                  uu(jdim,idim) = uu(jdim,idim) + fcilpc * rilvec(jdim)
 570           continue
 580        continue
 1000    continue
 1010 continue

      return
      end
c

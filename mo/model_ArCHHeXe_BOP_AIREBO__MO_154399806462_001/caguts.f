c-----------------------------------------------------------------------
c     caguts is the guts of the REBO potential.  (ca for carbon)
c     its main purpose is to generate a pair list for the REBO atom
c     pairs and get the various pieces of the REBO potential from lookup
c     tables or other routines....sjs 5/22/97
c-----------------------------------------------------------------------
c
      subroutine caguts(np, igfunc, fint, cube, ltors, r0, uu, tote)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np        = number of atoms (input)
c     igfunc    = style of angular term (input)
c     fint(a,d) = force on atom a in dimensiond (input/output)
c     cube(d)   = box length in dimension d (input)
c     ltors     = whether to include torsions (input)
c     r0(a,d)   = position of atom a in dimension d (input)
c     uu(d,d')  = d,d' element of internal virial tensor (input/output)
c     tote      = potential energy (output)

      integer np
      integer igfunc
      real*8 fint(npmax,3)
      real*8  cube(3)
      logical ltors
      real*8  r0(npmax,ndim)
      real*8  uu(ndim,ndim)
      real*8  tote

c     repulsive REBO interaction

      do 220 iprind = 1, nchch
         ijpair = ichch(iprind)

         iatom = ihalf(ijpair)
         jatom = jhalf(ijpair)

         tote = tote + repel(ijpair)

ceatom
c$$$c        split the pair energy between the two atoms
c$$$
c$$$         eatom(iatom) = eatom(iatom) + 0.5d0 * repel(ijpair)
c$$$         eatom(jatom) = eatom(jatom) + 0.5d0 * repel(ijpair)
ceatom

c        rpp holds -r_ij/|r_ij| * d(V^R_ij)/d|r_ij| = F^R_ij

         do 210 idim = 1, 3
            rpp = drepel(ijpair) * rijv(ijpair,idim)

            fint(iatom,idim) = fint(iatom,idim) + rpp
            fint(jatom,idim) = fint(jatom,idim) - rpp
            do 209 jdim = 1, 3
               uu(jdim,idim) = uu(jdim,idim) + rpp * rijv(ijpair,jdim)
 209        continue
 210     continue

 220  continue
     
c     Find atom type-specific coordination numbers of each atom,
c     allowing partial neighbors when they're in the switching region.

      do 510 iatom = 1, np
         ibegin = nabors(iatom)
         iend = nabors(iatom+1) - 1

         do 495 itype = 1, ntypes
            rnpls1(iatom,itype) = 1.d0
 495     continue

         do 500 ijpair = ibegin, iend
            jatom = jhalf(ijpair)
            if (lrebop(ijpair)) then
               rnpls1(iatom,iat2ty(jatom))
     .              = rnpls1(iatom,iat2ty(jatom)) + fij(ijpair)
            endif
 500     continue
 510  continue

c     sum over REBO pairs, and for each one
c     (1) call pibond to get the -bij V^A interaction, and
c     (2) call dihed to get the torsional interactions

      do 600 iprind = 1, nchch
         ijpair = ichch(iprind)

         call pibond(ijpair, igfunc, btot, fint, uu, tote)

         if (ltors .and.
     .        ijty(iat2ty(ihalf(ijpair)),iat2ty(jhalf(ijpair)))
     .        .eq. ijcc) then
            call dihed(ijpair, fint, uu, tote)
         endif
 
 600  continue

      return
c
      end
c
c-----------------------------------------------------------------------
c     getbij is a routine which takes a possibly non-legitimate REBO
c     pair (i.e. one that is out of the REBO cutoff) and temporarily
c     makes it look like it is within the cutoff at a specified radius.
c     the bij term in the REBO potential is then calculated for use in
c     the contingent LJ interaction....sjs 5/29/97
c-----------------------------------------------------------------------
c
      function getbij(iprlst, ijljpr, rijsep, rijvec, wwval, rsplin,
     .     fsplin, igfunc, vdw, fint, cube, r0, uu, lpbc)

      include 'common_files.inc'
      include 'common_pots.inc'

c     iprlst    = which pairlist to spoof: iljsaf or iljwat (input)
c     ijljpr    = which pair to spoof (input)
c     rijsep    = scalar distance between non-spoofed atoms (input)
c     rijvec(d) = distance between non-spoofed atoms in dimension d (input)
c     wwval     = bond connectivity value for the non-spoofed pair (input)
c     rsplin    = S(r) switching function (input)
c     fsplin    = connectivity switching function (input)
c     igfunc    = style of angular term to use (input)
c     vdw       = non-screened vdW interaction (input)
c     fint(a,d) = internal force on atom a in dimension d (input/output)
c     cube(d)   = box length in dimension d (input)
c     r0(a,d)   = position of atom a in dimension d (input)
c     uu(d,d')  = d,d' element of internal virial tensor (input/output)
c     lpbc(d)   = whether to use periodic boundary conds in dimension d (input)

      integer iprlst
      integer ijljpr
      real*8  rijsep
      real*8  rijvec(ndim)
      real*8  wwval
      real*8  rsplin
      real*8  fsplin
      integer igfunc
      real*8  vdw
      real*8  fint(npmax,ndim)
      real*8  cube(ndim)
      real*8  r0(npmax,ndim)
      real*8  uu(ndim,ndim)
      logical lpbc(ndim)

      getbij = 0.d0

c     set up the pair information needed for a spoofed call to pibond

      call fakepr(iprlst, ijljpr, rijsep, rijvec)

      ijrbpr = npairs + 1

      if (iprlst .eq. iljsaf) then
         iatom = iv(ijljpr)
         jatom = jv(ijljpr)
      else
         iatom = iljwv(ijljpr)
         jatom = jljwv(ijljpr)
      endif

      itype = iat2ty(iatom)
      jtype = iat2ty(jatom)

c     adjust the coordination number to take into account the new
c     fictitious neighbor that we've just created.

c     this does NOT count any i-k or j-l REBO interactions that would
c     be created in this fictitious pair.

      rnpls1(iatom,jtype) = rnpls1(iatom,jtype) + 1.d0 - wwval
      rnpls1(jatom,itype) = rnpls1(jatom,itype) + 1.d0 - wwval

c     call pibond to get bij

      call spoof(ijrbpr, igfunc, getbij, rsplin, fsplin, vdw,
     .           fint, cube, r0, uu, lpbc)

c     bij comes back too large by a factor of two

      getbij = 0.5d0 * getbij

c     repair rnpls1

      rnpls1(iatom,jtype) = rnpls1(iatom,jtype) - (1.d0 - wwval)
      rnpls1(jatom,itype) = rnpls1(jatom,itype) - (1.d0 - wwval)

      return

      end
c
c-----------------------------------------------------------------------
c     h2pot is a custom routine to evaluate the potential energy of
c     a H2 molecule at a specified bond length. Used mainly as a wrapper
c     because pair potential lookup tables are not accessible from
c     wrapper routines....sjs 5/23/07
c-----------------------------------------------------------------------
c
      real*8 function h2pot(r)

      include 'common_files.inc'
      include 'common_pots.inc'

      real*8 r

c     local variables

      integer it
      real*8  rt

      rt = r / ddtab(ihyd,ihyd)
      it = int(rt) + 1

      h2pot = 0.d0
      h2pot = h2pot + rtable(ihyd,ihyd,it)
     .     + (rtable(ihyd,ihyd,it+1)-rtable(ihyd,ihyd,it))
     .     * (rt - it + 1)
      h2pot = h2pot - 2.d0 * (atable(ihyd,ihyd,it)
     .     + (atable(ihyd,ihyd,it+1) - atable(ihyd,ihyd,it))
     .     * (rt - it + 1))

      return
      end
c



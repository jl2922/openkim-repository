c-----------------------------------------------------------------------
c     subroutine getmu calculates the system dipole moment.
c     ...sjs 5/29/08
c-----------------------------------------------------------------------
c
      subroutine getmu(np, r0, q0, smu)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np          = number of particles
c     r0(i,d)     = Cartesian coordinates of atom i in dimension d
c     q0(a)       = charge of atom a
c     smu(d)      = system dipole in dimension d (e A)

      integer np
      real*8  r0(npmax,ndim)
      real*8  q0(npmax)
      real*8  smu(ndim)

c     local variables

      integer iatom, idim

c     calculate the system dipole moment
c     M = sum_i qi ri

      do 110 idim = 1, ndim
         smu(idim) = 0.d0
         do 100 iatom = 1, np
            smu(idim) = smu(idim) + q0(iatom) * r0(iatom,idim)
 100     continue
 110  continue

      return
      end
c
c-----------------------------------------------------------------------
c     setkvc sets up the array of k vectors, to be used in Ewald
c     summation. Ported from earlier sjs version....sjs 7/25/07
c-----------------------------------------------------------------------
c
      subroutine setkvc(cube, vol)

      include 'common_files.inc'
      include 'common_pots.inc'

c     cube = box length
c     vol  = box volume

      real*8 cube(ndim)
      real*8 vol

c     local variables

      integer idim, jdim, nbig
      real*8 rk2prt(0:ndim),
     .     ewkmx2, rkpfac, twopi, volfac

c     initialize some stuff

      ewkmx2 = ewlkmx * ewlkmx
      rkpfac = 0.25d0 / ewlkp2
      twopi = 2.d0 * pi
      do 100 idim = 1, ndim
         rntok(idim) = twopi / cube(idim)
 100  continue
      volfac = 4.d0 * pi / vol

c     find the largest possible lattice vectors in each direction

      nbig = 0
      do 200 idim = 1, ndim
         nlamax(idim) = int(ewlkmx / rntok(idim))
         if (nlamax(idim) .gt. nbig) then
            nbig = nlamax(idim)
         endif
 200  continue

      if (nbig .gt. maxk1d) then
         write(isterr, *) 'ERROR! largest lattice vector component is ',
     .        nbig, ', which exceeds the maximum of ', maxk1d
         write(isterr, *) '       decrease box size, or increase ',
     .        'maxk1d and recompile.'
         call ioclos
         stop
      endif

c     counting only one of any inversion-symmetric pair, find all
c     vectors smaller than ewlmx, except the zero vector

      if (ndim .ne. 3) then
         write(isterr, *) 'ERROR! ndim = ', ndim
         write(isterr, *) '       Ewald k vector calculations in ',
     .        'setkvc only work in 3 dimensions'
         call ioclos
         stop
      endif

c     loop over all lattice vectors in a cubic lattice

c     initialize

      numkvc = 1
      rk2prt(0) = 0.d0
      do 300 idim = 1, ndim-1
         nvec(numkvc,idim) = 0
         rkvec(numkvc,idim) = rntok(idim) * nvec(numkvc,idim)
         rk2prt(idim) = rk2prt(idim-1)
     .        + rkvec(numkvc,idim) * rkvec(numkvc,idim)
 300  continue

c     start with (0,0,1)

      nvec(numkvc,ndim) = 1
      rkvec(numkvc,ndim) = rntok(ndim) * nvec(numkvc,ndim)
      rk2prt(ndim) = rk2prt(ndim-1)
     .     + rkvec(numkvc,ndim) * rkvec(numkvc,ndim)

c     keep this k vector only if it is inside the cutoff

 310  if (rk2prt(ndim) .le. ewkmx2) then
         rk2fac(numkvc) = rk2prt(ndim) * rkpfac
c        2 * 1/2 * 1/pi/V * 4 pi^2 / k^2 * exp(-k^2 / 4 K^2)
         rk2exp(numkvc) = volfac / rk2prt(ndim) * exp(-rk2fac(numkvc))
         numkvc = numkvc + 1
         if (numkvc .gt. maxkvc) then
            write(isterr, *) 'setkvc: ERROR! more than ', maxkvc,
     .           ' k vectors found.'
            write(isterr, *) '        increase maxkvc or decrease max ',
     .        'k vector in .in file and recompile'
            call ioclos
            stop
         endif
         do 350 idim = 1, ndim
            nvec(numkvc,idim) = nvec(numkvc-1,idim)
            rkvec(numkvc,idim) = rkvec(numkvc-1,idim)
 350     continue
      endif
      
c     go to the next lattice vector

c     increment last element of n vector
      nvec(numkvc,ndim) = nvec(numkvc,ndim) + 1
c     roll over, if needed
      if (nvec(numkvc,ndim) .gt. nlamax(ndim)) then
         do 400 idim = ndim - 1, 1, -1
            nvec(numkvc,idim) = nvec(numkvc,idim) + 1
            if (nvec(numkvc,idim) .le. nlamax(idim)) then
               go to 410
            endif
 400     continue
         go to 500
 410     rkvec(numkvc,idim) = rntok(idim) * nvec(numkvc,idim)
         rk2prt(idim) = rk2prt(idim-1)
     .        + rkvec(numkvc,idim) * rkvec(numkvc,idim)
         do 420 jdim = idim+1, ndim
            nvec(numkvc,jdim) = -nlamax(jdim)
            rkvec(numkvc,jdim) = rntok(jdim) * nvec(numkvc,jdim)
            rk2prt(jdim) = rk2prt(jdim-1)
     .           + rkvec(numkvc,jdim) * rkvec(numkvc,jdim)
 420     continue
      else
         rkvec(numkvc,ndim) = rntok(ndim) * nvec(numkvc,ndim)
         rk2prt(ndim) = rk2prt(ndim-1) 
     .        + rkvec(numkvc,ndim) * rkvec(numkvc,ndim)
      endif

      go to 310

c     the last one didn't work

 500  numkvc = numkvc - 1

      end
c      
c-----------------------------------------------------------------------
c     esguts evaluates the electrostatic part of the potential
c-----------------------------------------------------------------------
c     
      subroutine esguts(np, r0, q0, lbc, nbc, bq0, atomof, bctype, 
     .     Efield, lpbc,
     .     cube, lqsolv, lewald, lewsrf, tote, fint, chi, bcchi, uu)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np          = number of particles
c     r0(i,d)     = Cartesian coordinates of atom i in dimension d
c     q0(a)       = charge of atom a
c     lbc         = whether there are any bond charges
c     nbc         = number of bond charges
c     bq0(b)      = bond charge on bond b
c     atomof(b,n) = nth (n=1,2) atom of bond charge b
c     bctype(b)   = type (sigma, pi) of bond charge b
c     Efield(d)   = external electric field in dimension d (V/A)
c     lpbc(d)     = whether to apply periodic boundary conds in dimension d
c     cube(d)     = length of periodic box in dimension d
c     lqsolv      = whether matrix methods will be used to solve for charges
c     lewald      = whether to use Ewald summation
c     lewsrf      = whether to use the Ewald surface term
c     tote        = potential energy of the system
c     fint(i,d)   = internal force on atom i in dimension d
c     chi(i)      = electronegativity (dV/dq, neg electrochem pot) of atom i
c     bcchi(b)    = electronegativity (dV/dq, neg electrochem pot) of bond b
c     uu(d,d')    = internal virial tensor for dimensions d,d'

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
      logical lqsolv
      logical lewald
      logical lewsrf
      real*8  tote
      real*8  fint(npmax,ndim)
      real*8  chi(npmax)
      real*8  bcchi(nbcmax)
      real*8  uu(ndim,ndim)

c     external function

      real*8 bh, dbh, J

c     local variables

      integer iatom, ibc, 
c$$$     .     ibcind, 
     .     idim, ihead, itail, itype,
     .     jatom, jbc, jdim, jhead, jtail, jtype
c$$$      logical isbc(npmax)
      real*8  dr(ndim),
     .     dljdr, fcpc, littlj, onebyr, pepc, qiqj, qpetrm, rdotE, 
     .     rij, rsq

c     zero the electronegativities (neg. charge forces)

      do 10 iatom = 1, np
         chi(iatom) = 0.d0
 10   continue

c     zero the arrays used to store the Coulomb interaction, if we
c     need to fill it in order to evaluate charges with matrix methods

      if (lqsolv) then
         do 30 iatom = 1, np
            qvec(iatom) = 0.d0
            do 20 jatom = 1, np
               qmat(jatom,iatom) = 0.d0
 20         continue
 30      continue
      endif

      if (lbc) then

c     zero the bond charge electronegativities (neg. charge forces)

         do 32 ibc = 1, nbc
            bcchi(ibc) = 0.d0
 32      continue

c     zero the arrays used to store the Coulomb interaction between bond
c     charges, if we need to fill it in order to evaluate bond charges with
c     matrix methods

         if (lqsolv) then
            do 34 ibc = 1, nbc
               bcqvec(ibc) = 0.d0
               do 33 jbc = 1, nbc
                  bcqmat(jbc,ibc) = 0.d0
 33            continue
 34         continue
         endif

      endif

c     single charge contributions to the electrostatic energy
c     and charge forces. There is no position dependence, so no
c     spatial forces.

c     qi are atomic charges

c     Vii = chi^0 qi + 1/2 J^0 qi^2
c     chi = dVii/dqi = chi^0 + J^0 qi
c     J = d2V/dqi2 = J^0
c     b = linear term in q = chi^0

c     calculate these even if using bond charges, as they will be used
c     later to calculate
c     chi_a = dVii/dqa = dVii/dqi dqi/dqa
c           = cia ( chi^0 + J^0 qi )

      do 50 iatom = 1, np
         itype = iat2ty(iatom)
         tote = tote + chi0(itype) * q0(iatom)
     .        + 0.5d0 * J0(itype) * q0(iatom) * q0(iatom)
         chi(iatom) = chi(iatom) + chi0(itype)
     .        + J0(itype) * q0(iatom)
         qvec(iatom) = qvec(iatom) + chi0(itype)
 50   continue

      if (lqsolv) then
         do 60 iatom = 1, np
            itype = iat2ty(iatom)
            qmat(iatom,iatom) = qmat(iatom,iatom) + J0(itype)
 60      continue
      endif

c     interaction between the charges and the electric field.

c     qi are atomic charges

c     Vi = -q ri.E
c     chi = dVi/dqi = -ri.E
c     F = -dV/dxi = qi E
c     b = linear term in qi = -ri.E

c     calculate these even if using bond charges, as they will be used
c     later to calculate
c     chi_a = dVi/dqa = dVi/dqi dqi/dqa
c           = -cia ri.E

      do 80 iatom = 1, np
         rdotE = 0.d0
         do 70 idim = 1, ndim
            rdotE = rdotE + r0(iatom,idim) * Efield(idim)
            fint(iatom,idim) = fint(iatom,idim)
     .           + q0(iatom) * Efield(idim)
 70      continue
         tote = tote - q0(iatom) * rdotE
         chi(iatom) = chi(iatom) - rdotE
         qvec(iatom) = qvec(iatom) - rdotE
 80   continue

c     bond hardness interaction

      if (lbc) then

c     q are bond charges

c     V = 1/2 q^2 j(|r-r'|)
c     chi = dV/dq = q j(|r-r'|)
c     F = -dV/dx = 1/2 q^2 dj/dr (x-x') / r
c     F' = -dV/dx' = 1/2 q^2 dj/dr (x'-x) / r
c     b = linear term in q = 0
c     J = d2V/dqi^2 = j(|r-r'|)

         do 88 ibc = 1, nbc
            iatom = atomof(ibc,1)
            jatom = atomof(ibc,2)
            itype = iat2ty(iatom)
            jtype = iat2ty(jatom)
            rsq = 0.d0
            do 82 idim = 1, ndim
               dr(idim) = r0(iatom,idim) - r0(jatom,idim)
               if (lpbc(idim)) then
                  dr(idim) = dr(idim)
     .                 - cube(idim) * anint(dr(idim) / cube(idim))
               endif
               rsq = rsq + dr(idim) * dr(idim)
 82         continue
            rij = sqrt(rsq)
            write(isterr, *) 'store rij(bond chg) instead of recalcing'
            littlj = bh(itype, jtype, bctype(ibc), rij)
            tote = tote + 0.5d0 * bq0(ibc) * bq0(ibc) * littlj
            bcchi(ibc) = bcchi(ibc) + bq0(ibc) * littlj
            dljdr = dbh(itype, jtype, bctype(ibc), rij)
            fcrpc = 0.5d0 * bq0(ibc) * bq0(ibc) * dljdr / rij
            do 85 idim = 1, ndim
               fcpc = fcrpc * dr(idim)
               fint(iatom,idim) = fint(iatom,idim) + fcpc
               fint(jatom,idim) = fint(jatom,idim) - fcpc
 85         continue
            if (lqsolv) then
               bcqmat(ibc,ibc) = bcqmat(ibc,ibc) + littlj
            endif
 88      continue

      endif
      	
c     evaluate the Coulombic interaction between all atom pairs, storing
c     it in an array in case it is needed later to solve for the charges

      write(isterr, *) 'eventually, double loop needs to be ',
     .     'optimized to avoid recalculating distances'

      do 310 iatom = 1, np
         itype = iat2ty(iatom)

c     look ahead to see which j atoms correspond to a bond charge
c     hopefully this extra work makes up for not having to recalculate
c     bond lengths of bond charges.  For big enough systems, it won't.

c     This is work in progress...

c$$$         do 89 jatom = 1, np
c$$$            isbc(jatom) = .false.
c$$$ 89      continue
c$$$         do 90 ibcind = 1, nhbcof(iatom)
c$$$            ibc = hbcof(iatom,ibcind)
c$$$            isbc(atomof(ibc,2)) = .true.
c$$$ 90      continue
c$$$         do 92 ibcind = 1, ntbcof(iatom)
c$$$            ibc = tbcof(iatom,ibcind)
c$$$            isbc(atomof(ibc,2)) = .true.
c$$$ 92      continue

         do 300 jatom = iatom + 1, np
            jtype = iat2ty(jatom)
            rsq = 0.d0
            do 100 idim = 1, ndim
               dr(idim) = r0(iatom,idim) - r0(jatom,idim)
               if (lpbc(idim)) then
                  dr(idim) = dr(idim)
     .                 - cube(idim) * anint(dr(idim) / cube(idim))
               endif
               rsq = rsq + dr(idim) * dr(idim)
 100        continue
            rij = sqrt(rsq)

            qiqj = q0(iatom) * q0(jatom)

            if (.not. lewald) then
               
               pepc = J(itype, jtype, rij)

               qpetrm = qiqj * pepc
               
c     V = qi qj J(rij)
c     or in terms of bond charges,
c      = sum(cia qa) sum(cjb qb) J(rij),
c     where cia = 1 if atom i is first atom of bond charge a
c           cia = -1 if atom i is second (primed) atom of bond charge a
c           cia = 0 if atom i does not participate in bond charge a.
c     The first form is easier even if using bond charges

               tote = tote + qpetrm

c     qi are atomic charges, qa are bond charges

c     chi_i = dVij/dqi = qj J(rij)
c     chi_j = dVij/dqj = qi J(rij)
c     calculate these even if using bond charges, as they will be used
c     later to calculate
c     chi_a = dVij/dqa = dVij/dqi dqi/dqa + dVij/dqj dqj/dqa
c           = cia qj J(rij) + cja qi J(rij)
c           = cia chi_i + cja chi_j

               if (iseety(itype)) then
                  chi(iatom) = chi(iatom) + q0(jatom) * pepc
               endif
               if (iseety(jtype)) then
                  chi(jatom) = chi(jatom) + q0(iatom) * pepc
               endif
                  
c     J_ij = d2Vij/(dq_i dq_j)
c     calculate this even if using bond charges, as it will be
c     used later to calculate
c     J_ab = d2Vij/(dq_a dq_b) = d2Vij/(dq_i dq_j) dq_i/dq_a dq_j dq_b
c          = c_ia c_jb J(rij)

               if (lqsolv) then
                  qmat(iatom,jatom) = qmat(iatom,jatom) + pepc
                  qmat(jatom,iatom) = qmat(jatom,iatom) + pepc
               endif

c              forces
c              Fi = -dV/dxi = -qi qj dJ(rij)/d(rij) drij/dxi
c                           = -qi qj dJ/dr (xi - xj) / rij
c              Fj = -dV/dxj = -qi qj dJ/dr (xj - xi) / rij

               fcrpc = -qiqj * dJdr(itype, jtype, rij) / rij
               do 120 idim = 1, ndim
                  fcpc = fcrpc * dr(idim)
                  fint(iatom,idim) = fint(iatom,idim) + fcpc
                  fint(jatom,idim) = fint(jatom,idim) - fcpc
                  do 110 jdim = 1, ndim
                     uu(jdim,idim) = uu(jdim,idim) + fcpc * dr(jdim)
 110              continue
 120           continue

            else

c              Vij = qi qj J(rij)
c                  = qi qj { [ J(rij) - 1/rij ] + 1/rij }
c              The qi qj [ J(rij)  - 1/rij ] is evaluated by direct
c              summation, and the qi qj / rij contribution is evaluated
c              by Ewald summation.

c              For the first term,
c              Vij = qi qj [ J(rij) - 1/rij ]

               onebyr = epsinv / rij

               if (rij .le. rjmax(itype,jtype)) then

                  pepc = J(itype, jtype, rij) - onebyr

                  qiqj = q0(iatom) * q0(jatom)
                  qpetrm = qiqj * pepc

                  tote = tote + qpetrm

c     chi_i = dVij/dqi = qj [ J(rij) - 1/rij ]
c     chi_j = dVij/dqj = qi [ J(rij) - 1/rij ]

c     chi_a = dVij/dqi dqi/dqa dVij/dqj dqj/dqa
c           = (cia qj + cja qi) [ J(rij) - 1/rij ]
c           = cia chi_i + cja chi_j

                  if (iseety(itype)) then
                     chi(iatom) = chi(iatom) + q0(jatom) * pepc
                  endif
                  if (iseety(jtype)) then
                     chi(jatom) = chi(jatom) + q0(iatom) * pepc
                  endif

c     Jij = d2Vij/(dqi dqj) = J(rij) - 1/rij

                  if (lqsolv) then
                     qmat(iatom,jatom) = qmat(iatom,jatom) + pepc
                     qmat(jatom,iatom) = qmat(jatom,iatom) + pepc
                  endif

c                 forces
c                 Fi = -dV/dxi = -qi qj d[ J(rij) - 1/rij ]/drij drij/dxi
c                              = -qi qj [ dJ/dr + 1/r^2 ] (xi - xj) / rij
c                 Fj = -dV/dxj = -qi qj [ dJ/dr + 1/r^2 ] (xj - xi) / rij

                  fcrpc = -qiqj * (dJdr(itype, jtype, rij)
     .                 + onebyr / rij) / rij
                  do 210 idim = 1, ndim
                     fcpc = fcrpc * dr(idim)
                     fint(iatom,idim) = fint(iatom,idim) + fcpc
                     fint(jatom,idim) = fint(jatom,idim) - fcpc
                     do 200 jdim = 1, ndim
                        uu(jdim,idim) = uu(jdim,idim) + fcpc * dr(jdim)
 200                 continue
 210              continue

               endif

c              For the second term,
c              Vij = qi qj / rij
c                  = qi qj { erfc(K rij) / rij
c                          +sum_k[ 4pi^2/k^2 exp(-k^2/4K^2) cos(k.rij)]/(pi V)}
c              via the Ewald sum

c              For the real space term,
c              Vij = qi qj erfc(K rij) / rij

               pepc = herfc(ewlkap * rij) * onebyr

               qpetrm = qiqj * pepc

               tote = tote + qpetrm

c              chi_i = dVij/dqi = qj erfc(K rij) / rij
c              chi_j = dVij/dqj = qi erfc(K rij) / rij

               if (iseety(itype)) then
                  chi(iatom) = chi(iatom) + q0(jatom) * pepc
               endif
               if (iseety(jtype)) then
                  chi(jatom) = chi(jatom) + q0(iatom) * pepc
               endif

c     J_ij = d2Vij/(dqi dqj) = erfc(K rij) / rij

               if (lqsolv) then
                  qmat(iatom,jatom) = qmat(iatom,jatom) + pepc
                  qmat(jatom,iatom) = qmat(jatom,iatom) + pepc
               endif

c              forces
c              Fi = -dVij/dxi = -qi qj d[ erfc(K rij) / rij ]/drij drij/dxi
c                           =  qi qj [ erfc(K rij) / rij
c                                       + 2 K / sqrt(pi) exp(-K^2 rij^2) ]
c                                    (xi - xj) / rij^2
c              Fj = -dVij/dxj =  qi qj [ erfc(K rij) / rij
c                                       + 2 K / sqrt(pi) exp(-K^2 rij^2) ]
c                                    (xj - xi) / rij^2

               fcrpc = qiqj * (pepc
     .              + 2.d0 * ewlkap / sqrt(pi) * epsinv
     .                * exp(-ewlkp2 * rij * rij)) / rij / rij

               do 250 idim = 1, ndim
                  fcpc = fcrpc * dr(idim)
                  fint(iatom,idim) = fint(iatom,idim) + fcpc
                  fint(jatom,idim) = fint(jatom,idim) - fcpc
                  do 240 jdim = 1, ndim
                     uu(jdim,idim) = uu(jdim,idim) + fcpc * dr(jdim)
 240              continue
 250           continue

            endif

 300     continue
 310  continue

      if (lewald) then

c        precalculate cos(n k . r) and sin(n k . r) for many n and all r

         call cossin(np,r0)

c        reciprocal space part of Ewald sum

         call ewlrcp(np, r0, q0, lqsolv, tote, fint, chi, uu)

c        self-energy correction 

         call ewlslf(np, q0, lqsolv, tote, chi)

c        surface energy, called for vacuum boundary conditions

         if (lewsrf) then
            call ewlsrf(np, r0, q0, lqsolv, tote, fint, chi, uu)
         endif

      endif

      if (lbc) then

         do 410 ibc = 1, nbc

c     Anything accumulated in chi(i) is dV/dqi for some atom i.
c     Use these to calculate
c     dV/dqa = sum_i { dV/dqi dqi/dqa }
c            = sum_i { cia dV/dqi }
c            = dV/dqa - dV/dqa'
c     for bond charges a, comprising atoms a and a'

            ihead = atomof(ibc,1)
            itail = atomof(ibc,2)
            bcchi(ibc) = bcchi(ibc) + chi(ihead) - chi(itail)

c     Anything accumulated in qvec(i) is the coefficient multiplying
c     q_i for some atom i.
c     Use these to calculate the linear coefficients for bond charges a,
c     V = sum_i { b_i q_i } + nonlinear terms
c       = sum_i { b_i sum_a { cia qa } }
c       = sum_a { sum_i { b_i cia qa } }
c       = sum_a { b_a qa }
c     where
c     b_a = sum_i { b_i cia }
c         = b_a - b_a'
c     for bond charge a comprising atoms a and a'

            bcqvec(ibc) = bcqvec(ibc) + qvec(ihead) - qvec(itail)

            if (lqsolv) then

               do 400 jbc = ibc, nbc

c     Anything accumulated in qmat(i,j) is d2V/(dqi dqj) for atoms i,j.
c     Use these to calculate
c     d2V/(dqa dqb) = sum_i { sum_j { d2V/(dqi dqj) dqi/dqa  dqj/dqb } }
c                   = sum_i { sum_j { cia cjb d2V/(dqi dqj) } }
c                   = d2V/(dqa dqb) - d2V/(dqa dqb') - d2V/(dqa' dqb)
c                       + d2V/(dqa' dqb')
c     for bond charges a and b, comprising atoms a,a' and b,b'
               
                  jhead = atomof(jbc,1)
                  jtail = atomof(jbc,2)
                  bcqmat(ibc,jbc) = bcqmat(ibc,jbc) + qmat(ihead,jhead)
     .                 - qmat(ihead,jtail) - qmat(itail,jhead)
     .                 + qmat(itail,jtail)
                  bcqmat(jbc,ibc) = bcqmat(ibc,jbc)

 400           continue

            endif

 410     continue

      endif

      end
c
c-----------------------------------------------------------------------
c     ewlrcp evaluates the reciprocal space part of the Ewald sum.
c     ported from earlier sjs version....sjs 11/21/07
c-----------------------------------------------------------------------
c
      subroutine ewlrcp(np, r0, q0, lqsolv, tote, fint, chi, uu)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np        = number of particles
c     r0(i,d)   = Cartesian coordinates of atom i in dimension d
c     q0(i)     = charge of atom i
c     lqsolv    = whether matrix methods will be used to solve for charges
c     tote      = potential energy of the system
c     fint(i,d) = internal force on atom i in dimension d
c     chi(i)    = electronegativity (dV/dq, neg electrochem pot) of atom i
c     uu(d,d')  = internal virial tensor for dimensions d,d'

      integer np
      real*8  r0(npmax,ndim)
      real*8  q0(npmax)
      logical lqsolv
      real*8  tote
      real*8  fint(npmax,ndim)
      real*8  chi(npmax)
      real*8  uu(ndim,ndim)

c     local variables

      integer iatom, idim, ikvec, itype, jatom
      real*8  coskr(npmax), sinkr(npmax),
     .     coskyz, fcpc, fckpc, fckpc2, pepc,
     .     qcossm, qmatrm, qsinsm, sinkyz

c     This only works in 3 dimensions.
      if (ndim .ne. 3) then
         write(isterr, *)
     .        'ewlrcp: Ewald sums are not implemented for ndim = ',
     .        ndim, ' dimensions'
         write(isterr, *) 'change ndim to 3 or rewrite ewlrcp'
         call ioclos
         stop
      endif

c     loop over k vectors

      do 150 ikvec = 1, numkvc

         qcossm = 0.d0
         qsinsm = 0.d0

c        This is 2 * 1/2 * 1 / (pi V) * 4 pi^2 / k^2 * exp(-k^2 / 4K^2)

         pepc = rk2exp(ikvec) * epsinv

         do 100 iatom = 1, np

c           cos(ky y + kz z) = cos(ky y) cos(kz z) - sin(ky y) sin(kz z)

            coskyz = coskx(iatom,nvec(ikvec,2),2)
     .               * coskx(iatom,nvec(ikvec,3),3)
     .               - sinkx(iatom,nvec(ikvec,2),2)
     .                 * sinkx(iatom,nvec(ikvec,3),3)

c           sin(ky y + kz z) = sin(ky y) cos(kz z) + cos(ky y) sin(kz z)

            sinkyz = sinkx(iatom,nvec(ikvec,2),2)
     .               * coskx(iatom,nvec(ikvec,3),3)
     .               + coskx(iatom,nvec(ikvec,2),2)
     .                 * sinkx(iatom,nvec(ikvec,3),3)

c           cos(k.r) = cos(kx x + ky y + kz z)
c                    = cos(kx x) cos(ky y + kz z) - sin(kx x) sin(ky y + kz z)

            coskr(iatom) = coskx(iatom,nvec(ikvec,1),1) * coskyz
     .           - sinkx(iatom,nvec(ikvec,1),1) * sinkyz

c           sin(k.r) = sin(kx x) cos(ky y + kz z) + cos(kx x) sin(ky y + kz z)

            sinkr(iatom) = sinkx(iatom,nvec(ikvec,1),1) * coskyz
     .           + coskx(iatom,nvec(ikvec,1),1) * sinkyz

c           qcossm = sum_i qi cos(k.ri)
c           qsinsm = sum_i qi sin(k.ri)

            qcossm = qcossm + q0(iatom) * coskr(iatom)
            qsinsm = qsinsm + q0(iatom) * sinkr(iatom)

 100     continue

c        V_ks = 2 sum_k>0 c(k) 1/2 sum_i sum_j qi qj cos(k.rij)
c             = sum_k>0 c(k) [ sum_i qi exp(i k.ri) ]^2
c             = sum_k>0 c(k) { [ sum_i qi cos(k.ri) ]^2
c                              + [ sum_i qi sin(k.ri) ]^2 }
c        with
c           c(k) = 1/pi/V 4 pi^2 / k^2 exp(-k^2 / 4K^2)

         tote = tote + pepc
     .        * (qcossm * qcossm + qsinsm * qsinsm)

         do 140 iatom = 1, np

            itype = iat2ty(iatom)

c           chi_i = dV/dqi
c                 = sum_k>0 c(k) 2 { cos(k.ri) [ sum_j qj cos(k.rj) ]
c                                    + sin(k.ri) [ sum_j qj sin(k.rj) ] }

            fckpc = pepc * 2.d0

            if(iseety(itype)) then
               chi(iatom) = chi(iatom) + fckpc
     .              * (coskr(iatom) * qcossm
     .                 + sinkr(iatom) * qsinsm)
            endif

c     maybe more effcient to take the if statement out of the iatom loop
c     and write a new double loop

c           d2V / dqi dqj = sum_k>0 c(k) 2 [ cos(k.ri) cos(k.rj)
c                                            + sin(k.ri) sin(k.rj) ]

            if (lqsolv) then
               qmat(iatom,iatom) = qmat(iatom,iatom) + fckpc
     .              * (coskr(iatom) * coskr(iatom)
     .                 + sinkr(iatom) * sinkr(iatom))
               do 110 jatom = iatom+1, np
                  qmatrm = fckpc
     .                 * (coskr(iatom) * coskr(jatom)
     .                    + sinkr(iatom) * sinkr(jatom))
                  qmat(iatom,jatom) = qmat(iatom,jatom) + qmatrm
                  qmat(jatom,iatom) = qmat(jatom,iatom) + qmatrm
 110           continue
            endif

c           forces
c           Fi = -dV_self/dxi
c              = sum_k>0 c(k) 2 qi kx { sin(k.ri) [ sum_j qj cos(k.rj) ]
c                                       - cos(k.ri) [ sum_j qj sin(k.rj) ] }

            fckpc2 = fckpc * q0(iatom)
     .           * (sinkr(iatom) * qcossm - coskr(iatom) * qsinsm)

            do 130 idim = 1, ndim
               fcpc = fckpc2 * rkvec(ikvec,idim)
               fint(iatom,idim) = fint(iatom,idim) + fcpc
               do 120 jdim = 1, ndim
                  uu(jdim,idim) = uu(jdim,idim) + fcpc * r0(idim,jdim)
 120           continue
 130        continue

 140     continue

 150  continue

      end
c
c-----------------------------------------------------------------------
c     ewlslf evaluates the self energy correction to the Ewald sum.
c     ported from earlier SJS version....sjs 11/26/07
c-----------------------------------------------------------------------
c
      subroutine ewlslf(np, q0, lqsolv, tote, chi)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np        = number of particles
c     q0(i)     = charge of atom i
c     lqsolv    = whether matrix methods will be used to solve for charges
c     tote      = potential energy of the system
c     chi(i)    = electronegativity (dV/dq, neg electrochem pot) of atom i

      integer np
      real*8  q0(npmax)
      logical lqsolv
      real*8  tote
      real*8  chi(npmax)

c     local variables

      integer iatom, itype
      real*8  sqs

      sqs = 0.d0
      do 100 iatom = 1, np
         sqs = sqs + q0(iatom) * q0(iatom)
 100  continue

c     V_self = -K / sqrt(pi) * sum_i q_i^2

      tote = tote - rkbrp * sqs

c     charge forces
c     chi_i = dV/dq_i = -2 K / sqrt(pi) q_i

      do 200 iatom = 1, np
         itype = iat2ty(iatom)
         if (iseety(itype)) then
            chi(iatom) = chi(iatom) - r2kbrp * q0(iatom)
         endif
 200  continue

c     d2V/dq_i/dq_j = -2 K / sqrt(pi) delta_ij

      if (lqsolv) then
         do 300 iatom = 1, np
            itype = iat2ty(iatom)
            if (iseety(itype)) then
               qmat(iatom,iatom) = qmat(iatom,iatom) - r2kbrp
            endif
 300     continue
      endif

c     forces
c     Fi = dV/dxi = 0

      end
c
c-----------------------------------------------------------------------
c     ewlsrf calculates the surface energy term in the Ewald sum, and
c     associated forces.  This term is included in the Ewald energy if
c     the surrounding medium is vacuum. For a conducting medium, this
c     term is not included. Ported from an earlier version.
c     ...sjs 11/26/07
c-----------------------------------------------------------------------
c
      subroutine ewlsrf(np, r0, q0, lqsolv, tote, fint, chi, uu)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np        = number of particles
c     r0(i,d)   = Cartesian coordinates of atom i in dimension d
c     q0(i)     = charge of atom i
c     lqsolv    = whether matrix methods will be used to solve for charges
c     tote      = potential energy of the system
c     fint(i,d) = internal force on atom i in dimension d
c     chi(i)    = electronegativity (dV/dq, neg electrochem pot) of atom i
c     uu(d,d')  = internal virial tensor for dimensions d,d'

      integer np
      real*8  r0(npmax,ndim)
      real*8  q0(npmax)
      logical lqsolv
      real*8  tote
      real*8  fint(npmax,ndim)
      real*8  chi(npmax)
      real*8  uu(ndim,ndim)

c     local variables

      integer iatom, idim, jatom, jdim
      real*8  smu(ndim),
     .     dotprd, fcpc

c     calculate the system dipole moment

      call getmu(np, r0, q0, smu)

c     V_surf = 2 pi / 3 / V * M^2
c            = 2 pi / 3 / V * (sum_i qi ri)^2

      dotprd = 0.d0
      do 220 idim = 1, ndim
         dotprd = dotprd + smu(idim) * smu(idim)
 220  continue
      tote = tote + r2pb3v * dotprd

c     charge forces
c     chi_i = dV/dqi = 2 pi / 3 / V * 2 M . ri

      do 310 iatom = 1, np
         itype = iat2ty(iatom)
         if (iseety(itype)) then
            dotprd = 0.d0
            do 300 idim = 1, ndim
               dotprd = dotprd + smu(idim) * r0(iatom,idim)
 300        continue
            chi(iatom) = chi(iatom) + r4pb3v * dotprd
         endif
 310  continue

c     d2V/dqi/dqj = 2 pi / 3 / V * 2 ri . rj

      if (lqsolv) then
         do 420 iatom = 1, np
            do 410 jatom = iatom, np
               dotprd = 0.d0
               do 400 idim = 1, ndim
                  dotprd = dotprd + r0(iatom,idim) * r0(jatom,idim)
 400           continue
               qmat(iatom,jatom) = r4pb3v * dotprd
 410        continue
 420     continue
      endif

c     forces
c     -dV/dxi = -2 pi / 3V * 2 Mx qi

      do 520 idim = 1, ndim
         fcprt = -r4pb3v * smu(idim)
         do 510 iatom = 1, np
            fcpc = fcprt * q0(iatom)
            fint(iatom,idim) = fint(iatom,idim) + fcpc
            do 500 jdim = 1, ndim
               uu(jdim,idim) = uu(jdim,idim) + fcpc * r0(iatom,jdim)
 500        continue
 510     continue
 520  continue

      end
c
c-----------------------------------------------------------------------
c     cossin evaluates and stores cos(kx) and sin(kx) in a lookup table.
c     Ported from sjs from sjs' gtcos, not lh's cossin....7/30/07
c-----------------------------------------------------------------------
c
      subroutine cossin(np,r0)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np = number of particles
c     r0 = Cartesian coordinates of atoms

      integer np
      real*8 r0(npmax,ndim)

c     local variables

      integer iatom
      real*8  k1drpc

c     first, precalculate cos(2 pi n x / L) and sin(2 pi n x / L) for
c     various atoms, integers n, and dimensions. These will be used to
c     construct cos(k.r) and sin(k.r)

c     These cos and sin terms are calculated using deMoivre's identity,
c        (cos nx + i sin nx) = (cos x + i sin x)^n

c     waste some O(N) storage to store cos(-nx) and sin(-nx), in order
c     to save some CPU time elsewhere

      do 110 idim = 1, ndim
         do 100 iatom = 1, np
            k1drpc = rntok(idim) * r0(iatom,idim)
            coskx(iatom,1,idim) = cos(k1drpc)
            sinkx(iatom,1,idim) = sin(k1drpc)
            coskx(iatom,-1,idim) = coskx(iatom,1,idim)
            sinkx(iatom,-1,idim) = -sinkx(iatom,1,idim)
 100     continue
 110  continue

      do 220 idim = 1, ndim
         do 210 n = 2, nlamax(idim)
            do 200 iatom = 1, np
               coskx(iatom,n,idim) =
     .              coskx(iatom,n-1,idim) * coskx(iatom,1,idim)
     .              - sinkx(iatom,n-1,idim) * sinkx(iatom,1,idim)
               sinkx(iatom,n,idim) = 
     .              coskx(iatom,n-1,idim) * sinkx(iatom,1,idim)
     .              + sinkx(iatom,n-1,idim) * coskx(iatom,1,idim)
               coskx(iatom,-n,idim) = coskx(iatom,n,idim)
               sinkx(iatom,-n,idim) = -sinkx(iatom,n,idim)
 200        continue
 210     continue
 220  continue

      return
      end
c
c-----------------------------------------------------------------------
c     qsolv solves algebraically for the charges that minimize the 
c     energy, given J q = chi.
c     ...sjs 7/31/07
c-----------------------------------------------------------------------
c
      subroutine qsolv(np, cube, r0, lpbc, q0)

      include 'common_files.inc'
      include 'common_pots.inc'

c     np      = number of particles
c     cube(d) = box length in dimension d (A)
c     r0(a,d) = position of atom a in dimension d (A)
c     lpbc(d) = whether to use periodic boundary conditions in dimension d
c     q0(a)   = charge on atom a (e)

      integer np
      real*8  cube(ndim)
      real*8  r0(npmax,ndim)
      logical lpbc(ndim)
      real*8  q0(npmax)

c     local variables

      integer iatom, imlind, imolec, iqind, jqind, istatus
      logical lrxn
      real*8  qmatsm(nqmat,nqmat),
     .     qvecsm(npmax)

c     this is done to make the large qmatsm array statically
c     allocated, preventing some compilers (Absoft on MacOS X)
c     from allocating it on (and overflowing) the stack as an
c     automatic array

      save qmatsm

      if (np .gt. nqmat) then
         write(isterr, *) 'qsolv: np = ', np, ' atoms exceeds npmax = ',
     .        npmax
         write(isterr, *)
     .        '       increase npmax and recompile, or optimize ',
     .        'charges iteratively'
         call ioclos
         stop
      endif

c     make sure the molecular membership info is up to date

      if (lmolol) then
         call clustr(np, cube, r0, lpbc, lrxn)
      endif

c     for now, we assume molecule neutrality. implementing
c     intermolecular charge transfer or non-neutral molecules 
c     will require some changes.

c     for now, we assume all charges are variable (if any are).
c     mixed rigid and variables charges will require some changes.

c     consequently, arrays are not NxN, but (N-n)x(N-n), where N
c     is the number of variable charges and n is the number of
c     charge-constrained entities containing variable charges

c     construct a list used to loop over only the independent charges
c     (which may have changed since the last time we were here)

      nqind = 0
      do 100 imlind = 1, nmolec
         imolec = molno(imlind)
         iatom = imolec
 110     if (inext(iatom) .ne. -1) then
            nqind = nqind + 1
            iatom = inext(iatom)
            iqi2at(nqind) = iatom
            go to 110
         endif
 100  continue

c     qvec is already a full-sized length N vector containing the linear
c     terms multiplying the atomic charges, b_i = chi^0_i - E_ext . r_i,
c     i.e. the terms in the electrochemical potential that are independent
c     of other charges.
c     Condense this to a vector where every element is
c     -(b_i - b_1)
c     where the index "1" refers to the head atom in the corresponding
c     molecule.

      do 300 iqind = 1, nqind
         qvecsm(iqind) = qvec(molec(iqi2at(iqind)))
     .        - qvec(iqi2at(iqind))
 300  continue

c     qmat is already a full-sized NxN array of Jij = d2V/dqi/dqj
c     condense this to an array where every term is
c     Jij - J1j - Ji1 + J11
c     with the index "1" meaning the head atom in the corresponding
c     molecule

c     it is wasteful to use a second array instead of repacking qmat,
c     especially since these are memory-hogging arrays.
c     it would be complicated to do this, however, since we don't
c     know much about the atom ordering

      do 410 iqind = 1, nqind
         iatom = iqi2at(iqind)
         ihead = molec(iatom)
         do 400 jqind = 1, nqind
            jatom = iqi2at(jqind)
            jhead = molec(jatom)
            qmatsm(iqind,jqind) = qmat(iatom,jatom)
     .           - qmat(ihead,jatom) - qmat(iatom,jhead)
     .           + qmat(ihead,jhead)
 400     continue
 410  continue

c     solve for the charges using Gaussian elimination

      call dgesvw(nqind, nqmat, qmatsm, qvecsm, istatus)
      if (istatus .ne. 0) then
         write(isterr, *) 'qsolv: unable to solve for charges. ',
     .        'singular matrix?'
         call ioclos
         stop
      endif

c     extract the charges, solving for the dependent charges as we go.

      do 500 imlind = 1, nmolec
         imol = molno(imlind)
         q0(imol) = 0.d0
 500  continue

      do 510 iqind = 1, nqind
         iatom = iqi2at(iqind)
         ihead = molec(iatom)
         q0(iatom) = qvecsm(iqind)
         q0(ihead) = q0(ihead) - q0(iatom)
 510  continue

      return
      end
c
c-----------------------------------------------------------------------
c     bqsolv solves algebraically for the bond charges that minimize
c     the energy, given J q = chi.
c     ...sjs 6/13/08
c-----------------------------------------------------------------------
c
      subroutine bqsolv(nbc, bq0)

      include 'common_files.inc'
      include 'common_pots.inc'

c     nbc = number of bond charges
c     bq0(b) = bond charge on bond b (e)

      integer nbc
      real*8  bq0(nbc)

c     local variables

      integer ibc, istatus

      if (nbc .gt. nbcmax) then
         write(isterr, *) 'bqsolv: nbc = ', nbc,
     .        ' bond charges exceeds nbcmax = ', nbcmax
         write(isterr, *) '        increase nbcmax and recompile, or ',
     .        'optimize bond charges iteratively'
         call ioclos
         stop
      endif

c     bcqvec contains the coefficients of the terms linear in the bond
c     charge, b_a = chi^0_a - chi^0_a' - E_ext . (r_a - r_a'),
c     i.e. the terms in the electrochemical potential that are independent
c     of charge.

c     Change the sign of these terms.
c     (This should probably be done at the time bcqvec is filled.)

      do 50 ibc = 1, nbc
         bcqvec(ibc) = -bcqvec(ibc)
 50   continue
      
c     solve for the bond charges using Gaussian elimination

      call dgesvw(nbc, nqmat, bcqmat, bcqvec, istatus)

      if (istatus .ne. 0) then
         write(isterr, *) 'bqsolv: unable to solve for bond charges. ',
     .        'singular matrix?'
         call ioclos
         stop
      endif

      do 100 ibc = 1, nbc
         bq0(ibc) = bcqvec(ibc)
 100  continue

      return
      end
c
c-----------------------------------------------------------------------
c     bh calculates the bond hardness for a particular bond charge.
c     ...sjs 6/12/08
c-----------------------------------------------------------------------
c
      real*8 function bh(itype, jtype, btype, dist)

      include 'common_files.inc'
      include 'common_pots.inc'

c     itype = atom type of first atom in bond charge
c     jtype = atom type of second atom in bond charge
c     btype = type of bond charge (sigma or pi)
c     dist  = distance between atoms in bond charge (A)

      integer itype
      integer jtype
      integer btype
      real*8  dist

c     local variables

      integer it
      real*8  fc, rt

      if (btype .eq. isigma) then

         if (dist .le. dijmin(itype,jtype)) then

            bh = jbh(itype,jtype,btype)

         else if (dist .ge. dijmax(itype,jtype)) then
            
            bh = bigpos

         else

c     as the bond dissociates, make the bond hardness diverge to infinity

            rt = dist / ddtab(itype,jtype)
            it = int(rt) + 1

            if (lpdir) then
               write(isterr, *) 'bh: update for direct fc evaluation'
               call ioclos
               stop
            endif

            fc = tabfc(itype,jtype,it)
     .           + (tabfc(itype,jtype,it+1) - tabfc(itype,jtype,it))
     .           * (rt - it + 1)
            if (fc .le. 0.d0) then
               bh = bigpos
            else
               bh = jbh(itype,jtype,btype) + 1.d0 / fc
            endif

         endif

      else if (btype .eq. ipi) then

c     no cutoffs for pi bonds

         bh = jbh(itype,jtype,btype)

      else

         write(isterr, *) 'bh: unrecognized bond charge type ', btype
         call ioclos
         stop

      endif

      return
      end
c
c-----------------------------------------------------------------------
c     dbh calculates dj/dr for the bond hardness function j(r).
c     ...sjs 6/12/08
c-----------------------------------------------------------------------
c
      real*8 function dbh(itype, jtype, btype, dist)

      include 'common_files.inc'
      include 'common_pots.inc'

c     itype = atom type of first atom in bond charge
c     jtype = atom type of second atom in bond charge
c     btype = type of bond charge (sigma or pi)
c     dist  = ditance between atoms in bond charge (A)

      integer itype
      integer jtype
      integer btype
      real*8  dist

c     local variables

      integer it
      real*8  dfc, fc, rt

c     The bh() function uses j=1/f_c for the bond hardness function.
c     So dj/dr = d(1/f_c)/dr = -1/f_c^2 df_c/dr

      rt = dist / ddtab(itype,jtype)
      it = int(rt) + 1

      if (dist .le. dijmin(itype,jtype)) then

         dbh = 0.d0

      else if (dist .ge. dijmax(itype,jtype)) then

         dbh = 0.d0

      else

         if (lpdir) then
            write(isterr, *) 'dbh: update for direct fc evaluation'
            call ioclos
            stop
         endif

         fc = tabfc(itype,jtype,it)
     .        + (tabfc(itype,jtype,it+1) - tabfc(itype,jtype,it))
     .          * (rt - it + 1)
         dfc = tabdfc(itype,jtype,it)
     .        + (tabdfc(itype,jtype,it+1) - tabdfc(itype,jtype,it))
     .          * (rt - it + 1)

         if (fc .le. 0) then
            dbh = bigpos
         else
            dbh = -1.d0 / fc**2 * dfc
         endif

      endif
      return
      end
c

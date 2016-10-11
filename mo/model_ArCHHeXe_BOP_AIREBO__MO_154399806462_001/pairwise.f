c-----------------------------------------------------------------------
c     rebopr evaluates the pairwise REBO functions V^R(r), V^A(r) and
c     f_c(r)....sjs 5/7/10
c-----------------------------------------------------------------------
c
      subroutine rebopr(itype, jtype, rij,
     .     fc, dfcdr, vr, dvrdr, va, dvadr)

      include 'common_files.inc'
      include 'common_pots.inc'

c     itype = atom type of atom i (input)
c     jtype = atom type of atom j (input)
c     rij   = r_ij = distance between atoms i and j (output)
c     fc    = f_c(r_ij) (output)
c     dfcdr = df_c/dr at r_ij (output)
c     vr    = V^R(r_ij) (output)
c     dvrdr = dV^R/dr (output)
c     va    = V^A(r_ij) (output)
c     dvadr = dV^A/dr at r_ij (output)

      integer itype
      integer jtype
      real*8  rij
      real*8  fc
      real*8  dfcdr
      real*8  vr
      real*8  dvrdr
      real*8  va
      real*8  dvadr

c     local variables

      integer iterm
      real*8 dzblph(4), zblphi(4),
     .     df1, df2,
     .     dtemp, dtspln, dtzblf,
     .     dva1, dva2, dva3, dvcoul, dvm, dvrebo, dvv, dvzbl,
     .     ff1, ff2, tee, tsplin, twidth, tzblph, x,
     .     va1, va2, va3, vcoul, vrebo, vv, vzbl, zblau, zi, zj

      if (dijmax(itype,jtype) .eq. 0.d0) then
c     this is not a REBO pair
         fc = 0.d0
         dfcdr = 0.d0
         vr = 0.d0
         dvrdr = 0.d0
         va = 0.d0
         dvadr = 0.d0
      else

c     cutoff function

         if (rij .lt. dijmax(itype,jtype)
     .        .and. rij .gt. dijmin(itype,jtype)) then
            dtemp = pibyd(itype,jtype) * (rij - dijmin(itype,jtype))
            fc = (1.d0 + cos(dtemp)) * 0.5d0
            dfcdr = -pibyd(itype,jtype) * 0.5d0 * sin(dtemp)
         else if (rij .le. dijmin(itype,jtype)) then
            fc = 1.d0
            dfcdr = 0.d0
         else if (rij .ge. dijmax(itype,jtype)) then
            fc = 0.d0
            dfcdr = 0.d0
         else
            write(isterr, *) 'ERROR! dijmin(',itype,',',jtype,') = ',
     .           dijmin(itype,jtype), ' > dijmax(', itype, ',', jtype,
     .           ') = ', dijmax(itype,jtype)
            call ioclos
            stop
         endif
         
c     atomic number, for use in ZBL potential

         zi = kt2(itype)
         zj = kt2(jtype)

c     ZBL universal screening length
      
         zblau = zblcof * bohr / (zi**0.23d0 + zj**0.23d0)

         rsq = rij * rij

c     switching function for ZBL <-> REBO repulsion

         if (rij .ge. repmax(itype,jtype)) then
            tsplin = 0.d0
            dtspln = 0.d0
         else if (rij .gt. repmin(itype,jtype)) then
            twidth = repmax(itype,jtype) - repmin(itype,jtype)
            tee = (rij - repmin(itype,jtype)) / twidth
            tsplin = 1.d0 - tee ** 2 * (3.d0 - 2.d0 * tee)
            dtspln = 6.d0 * tee * (1.d0 - tee) / twidth
         else if (rij .le. repmin(itype,jtype)) then
            tsplin = 0.d0
            dtspln = 0.d0
         else
            write(isterr, *) 'ERROR! repmin(', itype, ',', jtype,
     .           ') = ', repmin(itype,jtype), ' > repmax(', itype, ',',
     .           jtype, ') = ', repmax(itype,jtype)
            call ioclos
            stop
         endif

c     attractive pair terms

         va1 = b1(itype,jtype) * exp(-beta1(itype,jtype)*rij)
         dva1 = -beta1(itype,jtype) * va1

         va2 = b2(itype,jtype) * exp(-beta2(itype,jtype)*rij)
         dva2 = -beta2(itype,jtype) * va2

         va3 = b3(itype,jtype) * exp(-beta3(itype,jtype)*rij)
         dva3 = -beta3(itype,jtype) * va3

         vv = (va1 + va2 + va3) * 0.5d0
         dvv = (dva1 + dva2 + dva3) * 0.5d0

         vrebo = vv * fc
         dvrebo = dvv * fc + vv * dfcdr

c     va holds -1/2 V^A_ij
c     dvadr holds -1/2 1/r dV^A/dr
c     These include the effects of both the REBO cutoff switching function
c     at the outer edge and the ZBL switching function at the inner edge

c     V^A = (1-S) * V^REBO

         va = (1.d0 - tsplin) * vrebo
         dvadr = ((1.d0 - tsplin) * dvrebo - dtspln * vrebo) / rij

c     repulsive pair terms

         ff1 = capa(itype,jtype) * exp(-alfa(itype,jtype)*rij)
         df1 = -alfa(itype,jtype) * ff1
         
         ff2 = (1.d0 + capq(itype,jtype) / rij)
         df2 = -capq(itype,jtype) / rsq

         vv = ff1 * ff2
         dvm = df1 * ff2 + ff1 * df2

c     ZBL universal screening potential (J. F. Ziegler, J. P. Biersack
c     and U. Littmark, The Stopping and Range of Ions in Matter,
c     Pergamon: NY 1977)

c     reduced radial coordinate

         x = rij / zblau

c     universal screening function

         tzblph = 0.d0
         dtzblf = 0.d0
         do 280 iterm = 1, 4
            zblphi(iterm) = zbla(iterm) * exp(-zblb(iterm) * x)
            dzblph(iterm) = -zblb(iterm) / zblau
     .           * zblphi(iterm)
            tzblph = tzblph + zblphi(iterm)
            dtzblf = dtzblf + dzblph(iterm)
 280     continue

c     Coulomb potential

         vcoul = zi * zj / rij
     .        * fe2esu * fe2esu / fA2cm * ferg2J * fJ2eV
         dvcoul = -vcoul / rij

c     vr holds V^R_ij
c     dvrdr holds -1/r dV^R/dr
c     These include the effects of both the REBO cutoff switching function
c     at the outer edge, and the ZBL switching function at the inner edge

c     V^R = S * V^ZBL + (1-S) * V^REBO
c         = V^REBO + S * (V^ZBL - V^REBO)

         vzbl = tzblph * vcoul * fc
         dvzbl = dtzblf * vcoul * fc + tzblph * dvcoul * fc
     .        + tzblph * vcoul * dfcdr

         vrebo = vv * fc
         dvrebo = dvm * fc + vv * dfcdr

         vr = vrebo + tsplin * (vzbl - vrebo)
         dvrdr = -(dvrebo + dtspln * (vzbl - vrebo)
     .        + tsplin * (dvzbl - dvrebo)) / rij

      endif

      end
c

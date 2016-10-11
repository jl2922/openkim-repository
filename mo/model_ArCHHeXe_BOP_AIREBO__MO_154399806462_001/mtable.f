c-----------------------------------------------------------------------
c     mtable fills lookup tables for the pair pieces of the REBO
c     potential:
c        tabfc contains the switching function fc
c        tabdfc contains the derivative of the switching function dfc/dr
c        atable contains -1/2 times the attractive piece V^A
c        datable contains -1/2 -1/r dV^A/dr
c        rtable contains the repulsive piece V^R
c        drtable contains -1/r dV^R/dr
c-----------------------------------------------------------------------
c     modified to extend the lookup tables out to where they need to be
c     evaluated to do the LJ switch, and included a non-switched version
c     of the lookup table....sjs 2/15/97
c-----------------------------------------------------------------------
c     
      subroutine mtable

      include 'common_files.inc'
      include 'common_pots.inc'

c     local variables

      integer ki, kj

c     generate lookup tables for bond-order potentials

      do 1000 ki = 1, ntypes
         do 900 kj = ki, ntypes
            ddtab(ki,kj) = dijmax(ki,kj) / dble(ntab - 2)
            ddtab(kj,ki) = ddtab(ki,kj)

            if (dijmax(ki,kj) .eq. 0.d0) then
               do 200 itab = 2, ntab-1
                  tabfc(ki,kj,itab)=0.0d0
                  tabfc(kj,ki,itab) = tabfc(ki,kj,itab)
                  tabdfc(ki,kj,itab) = 0.0d0
                  tabdfc(kj,ki,itab) = tabdfc(ki,kj,itab)
                  atable(ki,kj,itab) = 0.0d0
                  atable(kj,ki,itab) = atable(ki,kj,itab)
                  datable(ki,kj,itab) = 0.0d0
                  datable(kj,ki,itab) = datable(ki,kj,itab)
                  rtable(ki,kj,itab) = 0.0d0
                  rtable(kj,ki,itab) = rtable(ki,kj,itab)
                  drtable(ki,kj,itab) = 0.0d0
                  drtable(kj,ki,itab) = drtable(ki,kj,itab)
 200           continue
            else

               rc = 0.d0

               do 300 itab = 2, ntab-1
               
                  rc = rc + ddtab(ki,kj)

                  call rebopr(ki, kj, rc,
     .                 tabfc(ki,kj,itab), tabdfc(ki,kj,itab),
     .                 rtable(ki,kj,itab), drtable(ki,kj,itab),
     .                 atable(ki,kj,itab), datable(ki,kj,itab))

                  tabfc(kj,ki,itab) = tabfc(ki,kj,itab)
                  tabdfc(kj,ki,itab) = tabdfc(ki,kj,itab)

                  rtable(kj,ki,itab) = rtable(ki,kj,itab)
                  drtable(kj,ki,itab) = drtable(ki,kj,itab)

                  atable(kj,ki,itab) = atable(ki,kj,itab)
                  datable(kj,ki,itab) = datable(ki,kj,itab)

 300           continue
            endif 

c     the ( , ,1) element of a lookup table corresponds to r=0.  the 
c     ( , ,2) element corresponds to r=ddtab.

            atable(ki,kj,1) = 0.5d0 *
     .           (b1(ki,kj) + b2(ki,kj) + b3(ki,kj))
            atable(kj,ki,1) = atable(ki,kj,1)

c     the r=0 element of datable should be infinite, since it includes
c     a 1/r factor.  it is instead set to the r=ddtab element

            datable(ki,kj,1) = datable(ki,kj,2)
            datable(kj,ki,1) = datable(ki,kj,1)

c     the r=0 element of rtable should be infinite.  it is set instead
c     to the r=ddtab element.  ditto for d/dr.

            rtable(ki,kj,1) = rtable(ki,kj,2)
            rtable(kj,ki,1) = rtable(ki,kj,1)
            drtable(ki,kj,1) = drtable(ki,kj,2)
            drtable(kj,ki,1) = drtable(ki,kj,1)

            tabfc(ki,kj,1) = 1.0d0
            tabfc(kj,ki,1) = 1.0d0
            tabdfc(ki,kj,1) = 0.0d0
            tabdfc(kj,ki,1) = 0.0d0

c     the ( , ,ntab) element of these lookup tables is beyond dijmax and
c     must be here just in case r = dijmax exactly.

            atable(ki,kj,ntab) = 0.0d0
            atable(kj,ki,ntab) = atable(ki,kj,ntab)
            datable(ki,kj,ntab) = 0.0d0
            datable(kj,ki,ntab) = datable(ki,kj,ntab)

            rtable(ki,kj,ntab) = 0.0d0
            rtable(kj,ki,ntab) = rtable(ki,kj,ntab)
            drtable(ki,kj,ntab) = 0.0d0
            drtable(kj,ki,ntab) = drtable(ki,kj,ntab)

            tabfc(ki,kj,ntab) = 0.0d0
            tabfc(kj,ki,ntab) = tabfc(ki,kj,ntab)
            tabdfc(ki,kj,ntab) = 0.0d0
            tabdfc(kj,ki,ntab) = tabdfc(ki,kj,ntab)

c$$$            do 950 itab = 1, ntab
c$$$               iunit = 26 + ki * 4 + kj
c$$$               rc = (itab-1) * ddtab(ki,kj)
c$$$               write(iunit, '(5g20.10)') rc, rtable(ki,kj,itab),
c$$$     .              atable(ki,kj,itab), drtable(ki,kj,itab), 
c$$$     .              datable(ki,kj,itab)
c$$$ 950        continue

 900     continue
 1000 continue

      return
      end
c

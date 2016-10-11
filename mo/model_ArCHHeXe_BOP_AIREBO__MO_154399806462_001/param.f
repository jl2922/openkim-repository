c-----------------------------------------------------------------------
c     param sets lots of parameters for the REBO potential
c-----------------------------------------------------------------------

      subroutine param(rbuffr, iverb)

      include 'common_files.inc'
      include 'common_pots.inc'

c     iverb  = verbosity (0 = silent, 9 = debug)
c     rbuffr = buffer distance for pair list

      real*8  rbuffr
      integer iverb

c     local variables

      character messg*160,
     .     filnam*104
      integer itype, jtype

      dimension y(4), y1(4), y2(4), y12(4)
      real*8 xmin, xmax, ymin, ymax

      real*8 pirc(0:4,0:4,1:10), dpidi(0:4,0:4,1:10),
     .     dpidj(0:4,0:4,1:10), dpidn(0:4,0:4,1:10)

      real*8 tij(0:4,0:4,1:10), dtidi(0:4,0:4,1:10),
     .     dtidj(0:4,0:4,1:10), dtidn(0:4,0:4,1:10)

      dimension 
     .     f(0:1,0:1,0:1), fx(0:1,0:1,0:1), fy(0:1,0:1,0:1), 
     .     fz(0:1,0:1,0:1), fxy(0:1,0:1,0:1), fxz(0:1,0:1,0:1),
     .     fyz(0:1,0:1,0:1), fxyz(0:1,0:1,0:1)

c     Quintic splines g(cos theta) for the angular dependence of bij

c     g_C(cos theta)

      do 100 ig = 1, 4
         call qucof(spgcx(ig), spgcy(ig), spgcdy(ig), spgcdd(ig),
     .        spgc(1,ig))
 100  continue

c     gamma_C(cos theta)

c     caution:  the spgc*(5) elements are overwritten here and
c     irretrievably lost.  This is a slightly clumsy way of having
c     the g_C and gamma_C coefficients share one spline knot.

      spgcx(5) = spgcx(6)
      spgcy(5) = spgcy(6)
      spgcdy(5) = spgcdy(6)
      spgcdd(5) = spgcdd(6)
      call qucof(spgcx(4), spgcy(4), spgcdy(4), spgcdd(4),
     .     spgc(1,5))

c     g_H(cos theta)

      do 110 ig = 1, 3
         call qucof(spghx(ig), spghy(ig), spghdy(ig), spghdd(ig), 
     .        spgh(1,ig))
 110  continue

c     g_F(cos theta)

      do 120 ig = 1, 3
         call qucof(spgfx(ig), spgfy(ig), spgfdy(ig), spgfdd(ig),
     .        spgf(1,ig))
 120  continue

c     Zero bicubic spline coefficients

      do 230 i=1,ntypes
         do 220 l=1,10
            do 210 m=1,10
               do 200 j=1,16
                  clm(j,i,l,m) = 0.d0
 200           continue
 210        continue
 220     continue
 230  continue

c     bicubic spline

c     initialize bicubic spline values at knot points
 
      do 260 itype = 1, ntypes
         do 250 l = 1, 10
            do 240 m = 1, 10
               xh(itype,l,m) = 0.d0
 240        continue
 250     continue
 260  continue
   
c     This is a hack, for now. The eventual approach should be to
c     abandon the multidimensional P_ij cubic splines completely
c     (keeping a legacy version for use by HC REBO and AIREBO), 
c     and moving to the approach described in SJS' notebook on 11/11/08.

c     For now, P_CC is P_CC(N^C, N^H+N^F), which is a stopgap for HFC
c     systems.
      
      xh(icarb,1,3) = -0.027603d0
      xh(icarb,2,2) = -0.01096d0
      xh(icarb,2,3) =  0.00317953083d0
      xh(icarb,3,1) = -0.0005d0
      xh(icarb,3,2) =  0.00632624824d0
      xh(icarb,4,1) =  0.0161253646d0

c     These are P_CH values
c     Values for the following splines are in the following format:
c     xh(itype,#H+1,#C+1)
c     The smallest value is xh(itype,1,1).
c     Methane (CH4)  would have a spline of xh(ihyd,4,1) -MTK

      xh(ihyd,1,1)  =  0.d0
      xh(ihyd,2,1)  =  0.209336733d0
      xh(ihyd,3,1)  = -0.0644496154d0
      xh(ihyd,4,1)  = -0.303927546d0
      xh(ihyd,1,2)  =  0.01d0
      xh(ihyd,1,3)  = -0.122042146d0
      xh(ihyd,2,2)  = -0.125123401d0
      xh(ihyd,3,2)  = -0.298905246d0
      xh(ihyd,1,4)  = -0.307584705d0
      xh(ihyd,2,3)  = -0.300529172d0

c     P_CF values haven't been fit yet - Testing by MTK 6/08
c     Bonding types can be found in the REBO2 paper by Don Brenner 
c     Brenner, J. Phy.: Cond. Matter, 14, 783 (2002)
c     Hydrogen is replaced by fluorine - MTK

      xh(ifluor,1,1)  =  0.250d0
      xh(ifluor,2,1)  = 0.3d0
      xh(ifluor,4,1) =  -0.0450d0
      xh(ifluor,1,2)  =  1.2d0
      xh(ifluor,1,3)  =  1.5d0
      xh(ifluor,3,2)  =  0.070d0
      xh(ifluor,1,4)  =  0.0d0
      xh(ifluor,2,3)  =  0.30d0
      xh(ifluor,3,3)  =  0.0d0
      xh(ifluor,4,2)  =  0.10d0
      xh(ifluor,3,1)  =  0.425
      xh(ifluor,2,2)  =  0.825d0

c     these are set for consistency with DWB's REBO values
c     If a bond type is not present above, add it. 
c     The below bonding types do not affect the energy, but are
c     included for posterity.   -MTK 

      xh1(ihyd,3,1) = (xh(ihyd,4,1) - xh(ihyd,2,1)) / 2.0d0
      xh1(ihyd,2,2) = (xh(ihyd,3,2) - xh(ihyd,1,2)) / 2.0d0
      
      xh2(ihyd,2,2) = (xh(ihyd,2,3) - xh(ihyd,2,1)) / 2.0d0
      xh2(ihyd,1,3) = (xh(ihyd,1,4) - xh(ihyd,1,2)) / 2.0d0
 
      xh1(ifluor,3,1) = (xh(ifluor,4,1) - xh(ifluor,2,1)) / 2.0d0
      xh1(ifluor,2,2) = (xh(ifluor,3,2) - xh(ifluor,1,2)) / 2.0d0
      
      xh2(ifluor,2,2) = (xh(ifluor,2,3) - xh(ifluor,2,1)) / 2.0d0
      xh2(ifluor,1,3) = (xh(ifluor,1,4) - xh(ifluor,1,2)) / 2.0d0

c     calculate bicubic spline coefficients, using values at knots
      
      do 290 i = 1, ntypes
         do 280 l = 1, 9
            do 270 m = 1, 9
               xmin = l
               xmax = l + 1
               ymin = m
               ymax = m + 1
                  
               y(1) = xh(i,l,m)
               y(2) = xh(i,l+1,m)
               y(3) = xh(i,l+1,m+1)
               y(4) = xh(i,l,m+1)
               
               y1(1) = xh1(i,l,m)
               y1(2) = xh1(i,l+1,m)
               y1(3) = xh1(i,l+1,m+1)
               y1(4) = xh1(i,l,m+1)
               
               y2(1) = xh2(i,l,m)
               y2(2) = xh2(i,l+1,m)
               y2(3) = xh2(i,l+1,m+1)
               y2(4) = xh2(i,l,m+1)
               
               y12(1) = 0.d0
               y12(2) = 0.d0
               y12(3) = 0.d0
               y12(4) = 0.d0
               
               call bcucof(xmin, xmax, ymin, ymax, y, y1, y2, y12,
     .              clm(1,i,l,m))
               
 270        continue
 280     continue
 290  continue

c     zero tricubic spline coefficients

      do 340 j = 1, mxprty
         do 330 n = 1, 10
            do 320 m = 1, 10
               do 310 l = 1, 10
                  do 300 i = 1, 64
                     clmn(i,l,m,n,j) = 0.d0
 300              continue
 310           continue
 320        continue
 330     continue
 340  continue

c     build the spline coefficient arrays for the tricubic splines, 
c     from the values at the spline knots

c     choose a scheme for what to do with derivatives:

c     strictly according to DWB 2002 reference (and a few corrections
c     thereto)

c     1 => set all derivatives strictly according to DWB 2002
c                  (subject to a few corrections)
c     2 => zero all derivatives
c     3 => set all derivatives by finite difference

      iders = 1

      do 790 iprty = 1, mxprty

c     function and derivatives are zero unless explicitly stated
c     otherwise

         do 370 i = 0, 4
            do 360 j = 0, 4
               do 350 n = 1, 10
                  pirc(i,j,n) = 0.d0
                  dpidi(i,j,n) = 0.d0
                  dpidj(i,j,n) = 0.d0
                  dpidn(i,j,n) = 0.d0
 350           continue
 360        continue
 370     continue

         if (iprty .eq. ijcc) then

c           pi^rc_CC values at spline knots

c           values taken from DWB et al, JP:CM, 2002 and/or SJS et al, JCP,
c           2000.

c           function values:

            do 380 n = 3, 9
               pirc(0,0,n) = 0.0099172158d0
 380        continue
            pirc(1,0,1) = 0.04338699d0
            do 390 n = 2, 9
               pirc(1,0,n) = pirc(0,0,3)
 390        continue
            pirc(1,1,1) = 0.105d0
            pirc(1,1,2) = -0.0041775d0
            do 400 n = 3, 9
               pirc(1,1,n) = -0.0160856d0
 400        continue
            pirc(2,0,1) = 0.0493976637d0
            pirc(2,0,2) = -0.011942669d0
            do 410 n = 3, 9
               pirc(2,0,n) = pirc(0,0,3)
 410        continue
            pirc(2,1,1) = 0.0096495698d0
            pirc(2,1,2) = 0.03d0
            pirc(2,1,3) = -0.02d0
            pirc(2,1,4) = -0.0233778774d0
            pirc(2,1,5) = -0.0267557548d0
            do 420 n = 6, 9
               pirc(2,1,n) = -0.030133632d0
 420        continue
            pirc(2,2,1) = 0.09444957d0
            pirc(2,2,2) = 0.022d0
            pirc(2,2,3) = 0.03970587d0
            do 430 n = 4, 9
               pirc(2,2,n) = (9-n) / 6.d0 * pirc(2,2,3)
 430        continue
            pirc(3,0,1) = -0.119798935d0
            pirc(3,0,2) = pirc(3,0,1)
            do 440 n = 3, 9
               pirc(3,0,n) = pirc(0,0,3)
 440        continue
            do 450 n = 2, 9
               pirc(3,1,n) = -0.124836752d0
 450        continue
            do 460 n = 1, 9
               pirc(3,2,n) = -0.044709383d0
 460        continue

c           derivatives, from DWB '02
c           These may be overridden, depending on the value of iders
            
            dpidn(1,1,2) = -0.060543d0
            dpidi(2,1,1) = -0.0525d0
            dpidn(2,1,4) = -0.020044d0
            dpidn(2,1,5) = dpidn(2,1,4)
            do 470 n = 5, 9
               dpidi(2,1,n) = -0.054376d0
 470        continue
            do 480 n = 4, 8
               dpidn(2,2,n) = -0.006618d0
 480        continue
c           This value does not appear in DWB's paper, but should:
            dpidj(3,1,2) = 0.0375447764d0
            do 490 n = 2, 9
               dpidj(3,2,n) = 0.062418d0
 490        continue

         else if (iprty .eq. ijch) then

c           pi^rc_CH values at spline knots

c           values taken from SJS et al, JCP, 2000.  They differ somewhat from
c           the values in DWB et al, JP:CM, 2002.

            do 500 n = 1, 9
               pirc(1,1,n) = -0.1d0
 500        continue
            pirc(1,1,3) = -0.6d0
            do 510 n = 5, 9
               pirc(2,0,n) = -0.0090477875161288110d0
 510        continue
            pirc(2,1,2) = -0.5d00
            pirc(2,1,3) = pirc(2,1,2)
            do 520 n = 1, 9
               pirc(3,1,n) = -0.2d0
 520        continue
            pirc(3,1,2) = -0.25d0
            pirc(3,1,3) = pirc(3,1,2)

         else if (iprty .eq. ijhh) then

c           pi^rc_HH spline knots

c           value taken from DWB et al, JP:CM, 2002

            pirc(1,1,1) = 0.249831916d0

         else if (iprty .eq. ijcf) then

c           these values have not been fitted yet

         else if (iprty .eq. ijff) then

c           these values have not been fitted yet

         else if (iprty .eq. ijhf) then

c           these values have not been fitted yet

         else
            write(isterr, *) 'param: unrecognized pair type ', iprty
            call ioclos
            stop
         endif

c        make sure the function is flat at the top end, and symmetric

c        symmetry doesn't entirely make sense for CH and CF, but the
c        published parameters for CH are symmetric

         i = 4
         do 540 j = 0, 3
            do 530 n = 1, 9
               pirc(i,j,n) = pirc(i-1,j,n)
 530        continue
 540     continue

         do 570 i = 0, 3
            do 560 j = i+1, 4
               do 550 n = 1, 9
                  pirc(i,j,n) = pirc(j,i,n)
 550           continue
 560        continue
 570     continue

         i = 4
         j = 4
         do 580 n = 1, 9
            pirc(i,j,n) = pirc(i-1,j,n)
 580     continue

         n = 10
         do 600 i = 0, 4
            do 590 j = 0, 4
               pirc(i,j,n) = pirc(i,j,n-1)
 590        continue
 600     continue

c        The derivatives should also be symmetric

         do 630 i = 0, 2
            do 620 j = i+1, 3
               do 610 n = 1, 10
                  dpidi(i,j,n) = dpidj(j,i,n)
                  dpidj(i,j,n) = dpidi(j,i,n)
                  dpidn(i,j,n) = dpidn(j,i,n)
 610           continue
 620        continue
 630     continue

         if (iders .eq. 1) then

c           keep the values that we have set above

         else if (iders .eq. 2) then

c           zero the derivatives

            do 660 i = 0, 4
               do 650 j = 0, 4
                  do 640 n = 1, 10
                     dpidi(i,j,n) = 0.d0
                     dpidj(i,j,n) = 0.d0
                     dpidn(i,j,n) = 0.d0
 640              continue
 650           continue
 660        continue

         else if (iders .eq. 3) then
            
c           set the derivatives by finite difference

            do 690 i = 0, 3
               do 680 j = 0, 3
                  do 670 n = 1, 9
                     if (i .eq. 0) then
                        dpidi(i,j,n) = pirc(i+1,j,n) / 2.d0
                     else
                        dpidi(i,j,n) = (pirc(i+1,j,n) - pirc(i-1,j,n))
     .                       / 2.d0
                     endif
                     if (j .eq. 0) then
                        dpidj(i,j,n) = pirc(i,j+1,n) / 2.d0
                     else
                        dpidj(i,j,n) = (pirc(i,j+1,n) - pirc(i,j-1,n))
     .                       / 2.d0
                     endif
                     if (n .eq. 1) then
                        dpidn(i,j,n) = pirc(i,j,n+1) / 2.d0
                     else
                        dpidn(i,j,n) = (pirc(i,j,n+1) - pirc(i,j,n-1))
     .                       / 2.d0
                     endif
 670              continue
 680           continue
 690        continue

         else

            write(isterr, *) 'param: unknown treatment iders = ', iders,
     .           ' for treatment of pi^rc derivatives at spline knots'
            call ioclos
            stop
            
         endif

c       open a file to write out the spline coefficients

         if (iverb .ge. 9) then
            do 710 itype = 1, ntypes
               do 700 jtype = 1, ntypes
                  if (ijty(itype,jtype) .eq. iprty) then
                     if (symbol(itype)(1:1) .eq. ' ') then
                        if (symbol(jtype)(1:1) .eq. ' ') then
                           write(filnam, '(3a)') 'pirc.',
     .                          symbol(itype)(2:2), symbol(jtype)(2:2)
                        else
                           write(filnam, '(3a)') 'pirc.',
     .                          symbol(itype)(2:2), symbol(jtype)(1:2)
                        endif
                     else
                        if (symbol(jtype)(1:1) .eq. ' ') then
                           write(filnam, '(3a)') 'pirc.',
     .                          symbol(itype)(1:2), symbol(jtype)(2:2)
                        else
                           write(filnam, '(3a)') 'pirc.',
     .                          symbol(itype)(1:2), symbol(jtype)(1:2)
                        endif
                     endif
                     go to 720
                  endif
 700           continue
 710        continue
 720        open(itempf, file = filnam, err = 9990, iostat = ioval)
         endif

c        calculate the spline coefficients

         do 780 i = 0, 3
            xmin = i + 1
            xmax = i + 2
            do 770 j = 0, 3
               ymin = j + 1
               ymax = j + 2
               do 760 n = 1, 9
                  zmin = n
                  zmax = n + 1

c                 fill the function and derivative arrays

                  do 750 i2 = 0, 1
                     do 740 j2 = 0, 1
                        do 730 n2 = 0, 1
                           f(i2,j2,n2) = pirc(i+i2,j+j2,n+n2)
                           fx(i2,j2,n2) = dpidi(i+i2,j+j2,n+n2)
                           fy(i2,j2,n2) = dpidj(i+i2,j+j2,n+n2)
                           fz(i2,j2,n2) = dpidn(i+i2,j+j2,n+n2)
                           fxy(i2,j2,n2) = 0.d0
                           fxz(i2,j2,n2) = 0.d0
                           fyz(i2,j2,n2) = 0.d0
                           fxyz(i2,j2,n2) = 0.d0
 730                    continue
 740                 continue
 750              continue

c                 evaluate the coefficients

                  call tcucof(xmin, xmax, ymin, ymax, zmin, zmax,
     .                 f, fx, fy, fz, fxy, fxz, fyz, fxyz,
     .                 clmn(1,i+1,j+1,n,iprty))

c                 write the spline coefficients to a file  

                  if (iverb .ge. 9) then
                     write(itempf, *) i+1, j+1, n
                     write(itempf, '(3f24.10)')
     .                    (clmn(k,i+1,j+1,n,iprty),k=1,64)
                  endif

 760           continue
 770        continue
 780     continue

c        close the file containing the spline coefficients

         if (iverb .ge. 9) then
            close(itempf)
         endif

 790  continue

c$$$c     debugging code to measure a diagonal slice through pi^rc:
c$$$
c$$$      do 791 i = 0, 1000
c$$$         dr = (2.d0 - 1.2d0) / 1000
c$$$         r = 1.2d0 + i * dr
c$$$         call rebopr(6,6,r,fc,dfcdr,vr,dvrdr,va,dvadr)
c$$$         x = 2.d0 * fc + 1.d0
c$$$         yd = 2.d0 * fc + 1.d0
c$$$         z = 1.d0 + 8.d0 * fc * fc
c$$$         ix = int(x)
c$$$         iy = int(yd)
c$$$         iz = int(z)
c$$$         call tricub(x, yd, z, clmn(1,ix,iy,iz,ijcc), s, sx, sy, sz)
c$$$         write(istout, *) r, fc, s, sx, sy, sz
c$$$ 791  continue
c$$$      stop

c     tricubic spline pi^dh for rotation about multiple bonds

      do 830 n = 1, 10
         do 820 m = 1, 10
            do 810 l = 1, 10
               do 800 i = 1, 64
                  tlmn(i,l,m,n) = 0.d0
 800           continue
 810        continue
 820     continue
 830  continue

c     zero the function and its derivatives at the spline knots

      do 860 i = 0, 4
         do 850 j = 0, 4
            do 840 n = 1, 10
               tij(i,j,n) = 0.d0
               dtidi(i,j,n) = 0.d0
               dtidj(i,j,n) = 0.d0
               dtidn(i,j,n) = 0.d0
 840        continue
 850     continue
 860  continue

c     set the nonzero values
      
      tij(2,2,1) = -0.070280d0

      do 870 i = 2, 10
	 tij(2,2,i) = -0.008096d0
 870  continue

      do 930 i = 0, 3
	 xmin = i+1
	 xmax = i+2
	 do 920 j = 0, 3
	    ymin = j+1
	    ymax = j+2
	    do 910 n= 1, 9
	       zmin = n
	       zmax = n + 1

c              fill the function and derivative arrays

               do 900 i2 = 0, 1
                  do 890 j2 = 0, 1
                     do 880 n2 = 0, 1
                        f(i2,j2,n2) = tij(i+i2,j+j2,n+n2)
                        fx(i2,j2,n2) = dtidi(i+i2,j+j2,n+n2)
                        fy(i2,j2,n2) = dtidj(i+i2,j+j2,n+n2)
                        fz(i2,j2,n2) = dtidn(i+i2,j+j2,n+n2)
                        fxy(i2,j2,n2) = 0.d0
                        fxz(i2,j2,n2) = 0.d0
                        fyz(i2,j2,n2) = 0.d0
                        fxyz(i2,j2,n2) = 0.d0
 880                 continue
 890              continue
 900           continue

c              evaluate  the coefficients

               call tcucof(xmin, xmax, ymin, ymax, zmin, zmax,
     .              f, fx, fy, fz, fxy, fxz, fyz, fxyz,
     .              tlmn(1, i+1, j+1, n))
 910        continue
 920     continue
 930  continue

c      do 960 i = 0, 3
c	 xx = i + 1.d0
c	 do 950 j = 0, 3
c	    yy = j + 1.d0
c	    do 940 n= 1, 9
c	       zz = n 
c	       lindex = int(i+1)
c	       mindex = int(j+1)
c	       nindex = int(n)
c	       call tricub(xx, yy, zz, 
c     .              clmn(1,lindex, mindex, nindex, ijcc), 
c     .              rad, dradi, dradj, drdc)
c	       print *, i, j, n, 0, " rad", rad
c	       print *, i, j, n, 0, " dradi", dradi
c	       print *, i, j, n, 0, " dradj", dradj
c	       print *, i, j, n, 0, " drdc", drdc
c 940        continue 
c 950     continue 
c 960  continue 

      pibydq = pi / (bgqmax - bgqmin)

c     REBO parameters for covalent bonds

      do 1020 i = 1, ntypes
         do 1010 j = 1, ntypes
            do 1000 k = 1, ntypes
               ell(i,j,k) = 0.0d0
 1000       enddo
 1010    enddo
 1020 enddo

      pidt = pi/0.30d0

      do 1120 i=1,ntypes
         do 1110 j=1,ntypes
            do 1100 k=1,ntypes
               reg(i,j,k) = 1.0d0
 1100       continue
 1110    continue
 1120 continue

      rhh = 0.7415886997d0
      rch = 1.09d0
      ell(ihyd,ihyd,ihyd) = 4.0d0

      ell(ihyd,icarb,ihyd) = 4.0d0
      ell(ihyd,ihyd,icarb) = 4.0d0

      ell(ihyd,icarb,icarb) = 4.0d0

      reg(ihyd,icarb,ihyd) = exp(ell(ihyd,icarb,ihyd)*(rhh-rch))
      reg(ihyd,ihyd,icarb) = exp(ell(ihyd,ihyd,icarb)*(rch-rhh))

      prtrig = (0.5d0 * rbuffr) ** 2

c     these are parameters used in the torsional potential,
c        V(phi) = tor0 + tor10 * cos(phi) ** 10
c     if tor0 = -1/10 V and tor10 = 256/405 V, then an X3CCX3 molecule
c     (such as ethane) will have a torsional barrier of V.

c     default values

      do 1210 itype = 1, ntypes
         if (lrebot(itype)) then
            tor0(itype,itype) = 0.d0
            tor10(itype,itype) = 0.d0
            do 1200 jtype = itype+1, ntypes
               if (lrebot(jtype)) then
                  tor0(itype,jtype) = 0.d0
                  tor10(itype,jtype) = 0.d0

                  tor0(jtype,itype) = tor0(itype,jtype)
                  tor10(jtype,itype) = tor0(itype,jtype)
               endif
 1200       continue
         endif
 1210 continue

c     for H-C-C-H torsions, the potential reproduces the torsional
c     barrier in ethane

      tor0(ihyd,ihyd) = -hhtors / 10.d0
      tor10(ihyd,ihyd) = hhtors * 256.d0 / 405.d0

c     for C-C-C-H torsions, the potential reproduces the torsional
c     barrier in propane

      tor0(ihyd,icarb) = -chtors / 10.d0
      tor10(ihyd,icarb) = chtors * 256.d0 / 405.d0

      tor0(icarb,ihyd) = tor0(ihyd,icarb)
      tor10(icarb,ihyd) = tor10(ihyd,icarb) 

c     for C-C-C-C torsions, the potential reproduces the energy
c     difference between anti and gauche butane

      tor0(icarb,icarb) = -cctors / 10.d0
      tor10(icarb,icarb) = cctors * 256.d0 / 405.d0

c     for F-C-C-F torsions, the value should be fit to reproduce
c     the torsional barrier in perfluoroethane

c     for C-C-C-F torsions, the value should be fit to reproduce
c     the torsional barrier in perfluoropropane

c     these are the bounds on the switching function on cos(theta) 
c     which turns off the torsion interaction for nearly linear
c     bond angles.

c     torsion interaction begins to switch off at arccos(tthmax).
c     (this is 174.3 deg for tthmax = -0.995)
 
      tthmin = -1.d0
      tthmax = -0.995d0

c     the torsion interaction MUST be turned off completely for any
c     dihedral angles involving linear bonds, in order to avoid
c     singularities.

      if (tthmin .lt. -1.d0) then
         write(isterr, *) 'param: tthmin = ', tthmin
         write(isterr, *) '       set tthmin >= -1 and recompile'
      endif

      return

 9990 write(messg, '(2a)') 'ioopen: error opening file ',
     .     filnam(1:index(filnam, ' ')-1)
      call perror(messg)
      call ioclos
      stop

      end
c
c-----------------------------------------------------------------------
c     cuspa1 reads custom parameter values from files, when they
c     exist, for a specified parameter that is indexed by a single  atom
c     type....sjs 4/6/09
c-----------------------------------------------------------------------
c
      subroutine cuspa1(parray, name, unit)

      include 'common_files.inc'
      include 'common_pots.inc'

c     parray = array in which to store parameter values
c     name   = name of parameter file (prefix)
c     unit   = description of units assumed for user's parameter value

      real*8 parray(ntypes)
      character name*20
      character unit*20

c     local variables

      character filnam*104,
     .     dirnam*80,
     .     atname*3

      call getdir(dirnam)
      do 100 itype = 1, ntypes
         write(atname, '(a)') symbol(itype)(index(symbol(itype),' ')+1:)
         write(filnam, '(5a)') dirnam(1:index(dirnam,' ')-1),
     .        'param/', name(1:index(name,' ')-1), '.',
     .        atname(1:index(atname,' '))
         open(itempf, file=filnam, status='old', err=100)
         read(itempf, *) parray(itype)
         write(istout, *) 'Using custom parameter ',
     .        name(1:index(name,' ')-1), '_',
     .        atname(1:index(atname,' ')-1), ' = ', parray(itype), 
     .        unit
         close(itempf)
 100  continue
      return
      end
c
c-----------------------------------------------------------------------
c     cuspa2 reads custom parameter values from files, when they exist,
c     for a specified parameter that is indexed by two atom types.
c     ...sjs 4/7/09
c-----------------------------------------------------------------------
c
      subroutine cuspa2(parray, name, unit, symmet)

      include 'common_files.inc'
      include 'common_pots.inc'

c     parray = 2D array in which to store parameter values
c     name   = prefix for parameter file
c     unit   = description of units assumed for user's parameter value
c     symmet = whether the parameter is symmetric wrt exchange of indices

      real*8 parray(ntypes, ntypes)
      character name*20
      character unit*20
      logical symmet

c     local variables

      character filnam*104,
     .     dirnam*80,
     .     prname*5

      call getdir(dirnam)
      do 110 itype = 1, ntypes
         do 100 jtype = 1, ntypes
            write(prname, '(2a)')
     .           symbol(itype)(index(symbol(itype),' ')+1:),
     .           symbol(jtype)(index(symbol(jtype),' ')+1:)
            write(filnam, '(5a)') dirnam(1:index(dirnam,' ')-1),
     .           'param/', name(1:index(name,' ')-1), '.',
     .           prname(1:index(prname,' '))
            open(itempf, file=filnam, status='old', err=100)
            read(itempf, *) parray(itype,jtype)
            if (symmet) then
               parray(jtype,itype) = parray(itype,jtype)
            endif
            write(istout, *) 'Using custom parameter ',
     .           name(1:index(name,' ')-1), '_',
     .           prname(1:index(prname,' ')-1), ' = ',
     .           parray(itype,jtype), unit
            close(itempf)
 100     continue
 110  continue
      return
      end
c
c-----------------------------------------------------------------------
c     cuspa3 reads custom parameter values from files, when they exist,
c     for a specified parameter that is indexed by three atom types.
c     ...sjs 2/9/11
c-----------------------------------------------------------------------
c
      subroutine cuspa3(parray, name, unit, symmet)

      include 'common_files.inc'
      include 'common_pots.inc'

c     parray(r,s,t) = parameter value for types r,s,t (output)
c     name          = prefix for parameter file (input)
c     unit          = description of units assumed for user's parameter value
c     symmet        = whether parameter is symmetric wrt exchange of indices

      real*8 parray(ntypes, ntypes, ntypes)
      character name*20
      character unit*20
      logical symmet

c     local variables

      character filnam*104,
     .     dirnam*80,
     .     trname*7

      call getdir(dirnam)
      do 120 itype = 1, ntypes
         do 110 jtype = 1, ntypes
            do 100 ktype = 1, ntypes
               write(trname, '(3a)')
     .              symbol(itype)(index(symbol(itype),' ')+1:),
     .              symbol(jtype)(index(symbol(jtype),' ')+1:),
     .              symbol(ktype)(index(symbol(ktype),' ')+1:)
               write(filnam, '(5a)') dirnam(1:index(dirnam,' ')-1),
     .              'param/', name(1:index(name,' ')-1), '.',
     .              trname(1:index(trname,' '))
               open(itempf, file=filnam, status='old', err=100)
               read(itempf, *) parray(itype,jtype,ktype)
               if (symmet) then
                  parray(ktype,jtype,itype) = parray(itype,jtype,ktype)
               endif
               write(istout, *) 'Using custom parameter ',
     .              name(1:index(name,' ')-1), '_',
     .              trname(1:index(trname,' ')-1), ' = ',
     .              parray(itype,jtype,ktype), unit
               close(itempf)
 100        continue
 110     continue
 120  continue
      return
      end
c

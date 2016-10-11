c-----------------------------------------------------------------------
c     blkdat initializes some data that lives in common blocks. some
c     compilers object when this is done in the common block statement
c     itself. Only the common blocks in common_pots.inc are initialized
c     here. The common blocks in common_wrappers.inc are initialized in
c     initwrapper()...sjs 6/4/07
c-----------------------------------------------------------------------
c
      block data blkdat

      include 'common_files.inc'
      include 'common_pots.inc'

c     conversion from atomic number to atom type and vice versa

      data kt/ihyd,ihel,
     .     3*0,icarb,0,ioxy,ifluor,
     .     4*0,isili,3*0,iargon,
     .     13*0,igerm,4*0,
     .     17*0,ixenon,
     .     46*0/

c     some compilers (pgf77) choke on overlapping initializations, so
c     this line is not used
c     data kt2/100*0/

      data kt2(ihyd)/1/
      data kt2(ihel)/2/
      data kt2(icarb)/6/
      data kt2(ioxy)/8/
      data kt2(ifluor)/9/
      data kt2(isili)/14/
      data kt2(iargon)/18/
      data kt2(igerm)/32/
      data kt2(ixenon)/54/

c     whether an atom type uses the REBO potential

      data lrebot(ihyd)/.true./
      data lrebot(ihel)/.false./
      data lrebot(icarb)/.true./
      data lrebot(ioxy)/.true./
      data lrebot(ifluor)/.true./
      data lrebot(isili)/.true./
      data lrebot(iargon)/.false./
      data lrebot(igerm)/.true./
      data lrebot(ixenon)/.false./

c     symbols for elements

      data symbol(ihyd)/' H'/
      data symbol(ihel)/'He'/
      data symbol(icarb)/' C'/
      data symbol(ioxy)/' O'/
      data symbol(ifluor)/' F'/
      data symbol(isili)/'Si'/
      data symbol(iargon)/'Ar'/
      data symbol(igerm)/'Ge'/
      data symbol(ixenon)/'Xe'/

c     atom pair types

c     some compilers (pgf77) choke on overlapping initializations, so
c     this line is not used
c     data ijty/64*0.d0/

      data ijty(icarb,icarb)/ijcc/
      data ijty(icarb,ihyd)/ijch/
      data ijty(ihyd,icarb)/ijch/
      data ijty(ihyd,ihyd)/ijhh/
      data ijty(icarb,ifluor)/ijcf/
      data ijty(ifluor,icarb)/ijcf/
      data ijty(ifluor,ifluor)/ijff/
      data ijty(ihyd,ifluor)/ijhf/
      data ijty(ifluor,ihyd)/ijhf/

c     maximum valences for different atom types

      data maxval(ihyd)/1/
      data maxval(ihel)/0/
      data maxval(icarb)/4/
      data maxval(ioxy)/2/
      data maxval(ifluor)/1/
      data maxval(isili)/4/
      data maxval(iargon)/0/
      data maxval(igerm)/4/
      data maxval(ixenon)/0/

c     number of lone pairs for different atom types

      data numlps(ihyd)/0/
      data numlps(ihel)/1/
      data numlps(icarb)/0/
      data numlps(ioxy)/2/
      data numlps(ifluor)/3/
      data numlps(isili)/0/
      data numlps(iargon)/4/
      data numlps(igerm)/0/
      data numlps(ixenon)/0/

c     overpolarization penalty. set to 2 (rough draft) because this is
c     sufficient to make H + MH_n -> MH_n + H rather than MH_n+1 for
c     M = H, C, Si (with n = maxval(M))

      data opcnst(ihyd)/0.d0/
      data opcnst(ihel)/0.d0/
      data opcnst(icarb)/0.d0/
      data opcnst(ioxy)/2.d0/
      data opcnst(ifluor)/2.d0/
      data opcnst(isili)/2.d0/
      data opcnst(iargon)/0.d0/
      data opcnst(igerm)/2.d0/
      data opcnst(ixenon)/0.d0/

c     quintic spline coefficients for G_C(cos theta) and 
c     gamma_C(cos theta)

      data spgcx/-1.d0, -0.66666666666667d0, -0.5d0, -0.3333333333d0,
     .     1.d0, 1.d0/
      data spgcy/-0.01d0, 0.028207d0, 0.052804d0, 0.097321d0, 8.d0, 
     .     1.d0/
      data spgcdy/0.104d0, 0.131443d0, 0.17d0, 0.4d0,
     .     20.2436d0, 2.83457d0/
      data spgcdd/0.d0, 0.140229d0, 0.37d0, 1.98d0,
     .     43.9336d0, 10.2647d0/

c     parameters for ZBL universal screening function
c     (Ziegler et al, Stopping and Range of Ions in Matter, Pergamon,
c     1977).

      data zbla/0.1818d0, 0.5099d0, 0.2802d0, 0.02817d0/
      data zblb/3.2d0, 0.9423d0, 0.4028d0, 0.2016d0/
      
c     lookup table for spline region based on three-body angle
c     about a carbon atom

      data igc/16*4,2*3,2*2,5*1/

c     6th-order spline coefficients for G_H(cos theta)

      data spghx/-1.d0, -0.8333333333d0, -0.5d0, 1.d0/
      data spghy/11.2357, 12.5953d0, 16.8111d0, 19.9918d0/
      data spghdy/0.d0, 13.8543d0, 8.64123d0, 0.333013d0/
      data spghdd/115.115d0, 32.3618d0, -25.0617d0, -0.474189d0/

c     6th-order spline coefficients for G_F(cos theta)

      data spgfx/-1.d0, -0.8333333333d0, -0.5d0, 1.d0/
      data spgfy/11.2357, 12.5953d0, 16.8111d0, 19.9918d0/
      data spgfdy/0.d0, 13.8543d0, 8.64123d0, 0.333013d0/
      data spgfdd/115.115d0, 32.3618d0, -25.0617d0, -0.474189d0/

c     lookup table for spline region based on three-body angle
c     about a hydrogen atom

      data igh/18*3,4*2,3*1/

c     lookup table for spline region based on three-body angle
c     about a fluorine atom

      data igf/18*3,4*2,3*1/

c     initial zeroing of P_ij lookup table and its derivatives

      data xh/isizxh*0.0/,
     .     xh1/isizxh*0.0/,
     .     xh2/isizxh*0.0/

      end

      include 'parameters_both.inc'
      include 'parameters_pots.inc'

c     variables having to do with electronegativity equalization:
c     chi0(t)      = electronegativity of atom type t
c     qvec(i)      = linear coefficient of charge i (electroneg, E field)
c     drj(s,t)     = bin spacing for J(r) lookup table for atom types s and t
c     iqi2at(i)    = atom corresponding to ith independent charge
c     iseety(t)    = whether atom type t has EE parameters
c     ivalsh(t)    = valence shell for atom type t
c     J0(t)        = hardness of atom type t
c     jbh(s,t,b)   = bond hardness for a bond of type b between atom types s,t
c     ljdir        = whether to calculate J(r) directly (vs lookup table)
c     nqind        = number of independent charges
c     bcqmat(a,b)  = electrostatic interaction between bond charges a and b
c     qmat(i,j)    = electrostatic interaction between atoms i and j
c     rjmax(s,t)   = maximum r for J(r) lookup table for atoms types s and t
c     tabdj(s,t,k) = dJ(r)/dr at r=k*drj(s,t) for atom types s and t
c     tabj(s,t,k)  = J(r) at r=k*drj(s,t) for atom types s and t
c     zeta(t)      = Slater zeta orbital exponent for atom type t, 1/A

      logical iseety, ljdir
      real*8  J0, J0min, jbh

      common/eleceq/bcqmat(nqmat,nqmat), qmat(nqmat,nqmat),
     .     tabj(ntypes,ntypes,0:njtab),
     .     tabdj(ntypes,ntypes,0:njtab),
     .     bcqvec(nbcmax),
     .     jbh(ntypes,ntypes,nbctyp),
     .     drj(ntypes,ntypes), rjmax(ntypes,ntypes),
     .     chi0(ntypes), J0(ntypes), zeta(ntypes),
     .     J0min,
     .     iqi2at(npmax), qvec(npmax),
     .     ivalsh(ntypes),
     .     nqind,
     .     iseety(ntypes),
     .     ljdir

c     Ewald summation:
c     coskx(i,n,d) = cos^n(2 pi n (r_i)_d / L_d)
c     ewlkap = Ewald kappa (1/A)
c     ewlkmx = largest |k| for Ewald k-space sum, where
c              k_d = 2 pi n_d / L_d, d=x,y,z
c     ewlkp2 = kappa squared, for Ewald real-space sum
c     nlamax(d) = maximum lattice vector along dimension d in k space
c     numkvc = number of k vectors
c     nvec   = lattice vectors associated with k vectors
c     r2kbrp = 2 kappa / sqrt(pi)
c     r2pb3v = 2 pi / 3V
c     r4pb3v = 4 pi / 3V
c     rkbrp  = kappa / sqrt(pi)
c     rkvec  = k vectors
c     rk2exp = 1 / pi / V * 4 pi^2 / k^2 * exp(-k^2 / 4 / kappa^2)
c     rk2fac = k^2 / 4 / kappa^2
c     sinkx(i,n,d) = sin^n(2 pi n (r_i)_d / L_d)

      common/ewald/coskx(npmax,-maxk1d:maxk1d,ndim),
     .     sinkx(npmax,-maxk1d:maxk1d,ndim),
     .     rkvec(maxkvc,ndim),
     .     rk2exp(maxkvc), rk2fac(maxkvc),
     .     rntok(ndim),
     .     ewlkap, ewlkmx, ewlkp2, r2kbrp, r2pb3v, r4pb3v, rkbrp,
     .     nvec(maxkvc,ndim),
     .     nlamax(ndim),
     .     numkvc

c     variables having to do with molecular membership lists
c     inext    = next atom in molecule linked list
c     lmolol   = whether molecule list is old
c     molec(i) = molecule to which atom i belongs, also head of linked list
c     molno(i) = molecule number of ith molecule
c     natmol   = number of atoms in specified molecule
c     nmolec   = number of molecules
c     nmonat   = number of monatomic molecules
c     npolat   = number of atoms in polyatomic (non-monatomic) molecules
c     npolmo   = number of polyatomic (non-monatomic) molecules

      integer inext, molec, molno,
     .     natmol, nmolec, nmonat, npolat, npolmo
      logical lmolol

      common/molecp/
     .     inext(npmax), molec(npmax),
     .     molno(maxmol), natmol(maxmol),
     .     nmolec, nmonat, npolat, npolmo,
     .     lmolol

C Common blocks for empirical hydrocarbon and fluorocarbon potentials

      common/kxis2/igf(25), igh(25)
      common/xis12/xh(ntypes,10,10),xh1(ntypes,10,10),xh2(ntypes,10,10)

c     potential variables:

c     glamda(s,t,u) = BOP angular coeff for angle of type s-t-u
c     lhvang(s,t,u) = whether angular params exist for angle of type s-t-u
c     lhvcov(s,t)   = whether covalent bond params exist for types s and t
c     opcnst(t)     = overcoordination penalty constant for type t

      logical lhvang, lhvcov
      real*8  glamda, opcnst

      common/bondor/glamda(ntypes,ntypes,ntypes),
     .     b1(ntypes,ntypes), beta1(ntypes,ntypes),
     .     b2(ntypes,ntypes), beta2(ntypes,ntypes),
     .     b3(ntypes,ntypes),
     .     beta3(ntypes,ntypes),
     .     capa(ntypes,ntypes), alfa(ntypes,ntypes),
     .     capq(ntypes,ntypes), pibyd(ntypes,ntypes),
     .     dijmin(ntypes,ntypes), dijmax(ntypes,ntypes),
     .     r2rbmx(ntypes,ntypes),
     .     dijmx2(ntypes,ntypes),
     .     spgc(6,5),
     .     spgf(6,3), spgh(6,3),
     .     adb(ntypes), cdb(ntypes), cdb2(ntypes),
     .     ddb(ntypes), ddb2(ntypes), hdb(ntypes),
     .     opcnst(ntypes), xtn1(ntypes), xtn2(ntypes),
     .     spgcdd(6), spgcdy(6), spgcx(6), spgcy(6),
     .     spgfdd(4), spgfdy(4), spgfx(4), spgfy(4),
     .     spghdd(4), spghdy(4), spghx(4), spghy(4),
     .     pibydq,
     .     lhvang(ntypes,ntypes,ntypes),
     .     lhvcov(ntypes,ntypes)

      common/carbi/igc(25)

c     variables encoding some basic chemistry:

c     maxval(t)     = maximum valence for atom type t
c     numlps(t)     = number of lone pairs for atom type t
c     symbol(t) = chemical symbol of atom type t

      character symbol*2
      integer maxval, numlps

      common/chmstry/
     .     maxval(ntypes), numlps(ntypes),
     .     symbol(ntypes)


      common/barr/ell(ntypes,ntypes,ntypes),reg(ntypes,ntypes,ntypes),
     .     xxdb
      common/interp/clmn(64,10,10,10,mxprty),
     .     clm(16,ntypes,10,10),
     .     tlmn(64,10,10,10),
     .     pidt
      common/tors/tthmin, tthmax

c     lookup tables for the full REBO pair list (including buffer distance)
c     dexx1(p)  = -1/2 1/r dV^A/dr for pair p
c     drepel(p) = -1/r dV^R/dr for pair p
c     dww(p)    = df_c/dr for pair p
c     exx1(p)   = attractive term -1/2 V^A for pair p
c     fij(p)    = cutoff function f_c for pair p
c     ihalf(p)  = first of two atoms participating in pair p
c     jhalf(p)  = second of two atoms participating in pair p
c     repel(p)  = repulsive pairwise term V^R for pair p
c     lrebop(p) = whether pair p is a bona fide interacting REBO pair right now

      integer ihalf, jhalf
      logical lrebop
      real*8  dexx1, drepel, dww, exx1, fij, repel
      common/frplst/dexx1(nlmax), drepel(nlmax), dww(nlmax),
     .     exx1(nlmax), fij(nlmax), repel(nlmax),
     .     ihalf(nlmax), jhalf(nlmax),
     .     lrebop(nlmax)

      common/bond/rnpls1(npmax,ntypes)

c     variables having to do with REBO lookup tables:
c     atable(s,t,k)  = -1/2 V^A(r) at r=k*ddtab(s,t) for atom types s and t
c     datable(s,t,k) = -1/2 1/r dV^A/dr at r=k*ddtab(s,t) for types s and t
c     ddtab(s,t)     = bin spacing for lookup tables for atom types s and t
c     drtable(s,t,k) = -1/r dV^R/dr at r=k*ddtab(s,t) for atom types s and t
c     lljdir         = whether to calculate lj terms directly
c     lpdir          = whether to calculate pairwise terms directly
c     rtable(s,t,k)  = V^R(r) at r=k*ddtab(s,t) for atom types s and t
c     tabdfc(s,t,k)  = df_c/dr at r=k*ddtab(s,t) for atom types s and t
c     tabfc(s,t,k)   = f_c(r) at r=k*ddtab(s,t) for atom types s and t

      logical lljdir, lpdir

      common/chlist/atable(ntypes,ntypes,ntab),
     .     datable(ntypes,ntypes,ntab), drtable(ntypes,ntypes,ntab),
     .     rtable(ntypes,ntypes,ntab),
     .     tabdfc(ntypes,ntypes,ntab), tabfc(ntypes,ntypes,ntab),
     .     ddtab(ntypes,ntypes),
     .	   lljdir, lpdir

      common/zblrep/repmax(ntypes,ntypes), repmin(ntypes,ntypes)

c from common_lj_new.inc

      common/eatint/iv(nmabig), jv(nmabig),
     .     iljwv(nwlmax), jljwv(nwlmax),
     .     nljprs, nljwpr
      COMMON/csurf/epss(ntypes,ntypes), sigs(ntypes,ntypes),
     .     tau(ntypes),
     .     surf
      COMMON/ncsurf/ndir 
      common/ljlook/vljtb(maxbin,ntypes,ntypes),
     .     dvljtb(maxbin,ntypes,ntypes),
     .     dellj

c     LJ parameters

c     r2mxlj(s,t) = square of outer edge of LJ switch for atom types s & t
c     rmaxlj(s,t) = outer edge of LJ switch for atom types s & t
c     rminlj(s,t) = inner edge of LJ switch for atom types s & t
c     rslj(s,t)   = square of LJ neighbor list cutoff for types s & t
c     sig(s,t)    = square of LJ sigma for atom types s & t
c     sigcut      = outer edge of LJ switch, in sigma

      common/eatd_pot/bijmax(ntypes,ntypes), bijmin(ntypes,ntypes),
     .     eps(ntypes,ntypes),
     .     rlj0(ntypes,ntypes), rlj1(ntypes,ntypes),
     .     r2mxlj(ntypes,ntypes), rmaxlj(ntypes,ntypes),
     .     rminlj(ntypes,ntypes), rslj(ntypes,ntypes),
     .     sig(ntypes,ntypes),
     .     pvdw, sigcut

c from common_md.inc

      common/c0l/r0l(npmax,3)

c     number of things:
c     nnra   = number of non-REBO atoms
c     npairs = number of pairs in the REBO pair list
c     nrba   = number of REBO atoms

      common/number_pot/nnra, npairs, nrba

c     pair list info:
c     lcovol = whether covalent paiir list has expired
c     lposol = whether atoms have moved
c     lvdwol = whether van der Waals pair list has expired
c     rijv   = rij vector for a given atom pair
c     rcor   = rij distance for a given atom pair
c     ichch  = list of pairs which are actually REBO neighbors
c     i2chch = double-length list of pairs which are REBO neighbors
c     nabors = indexing of main REBO pair list by 1st half atom
c     n2strt = indexing of i2chch pair list by 1st half atom
c     prtrig = square of 1/2 buffer distance, used for pairlist updates
c     nchch  = number of actual REBO neighbors
c     n2chch = twice number of actual REBO neighbors

      logical lcovol, lposol, lvdwol
      common/prlist_pot/rijv(nlmax,3),
     .     rcor(nlmax),
     .     ichch(nlmax), i2chch(n2lmax),
     .     nabors(npmax+1), n2strt(npmax+1),
     .     prtrig,
     .     nchch, n2chch,
     .     lcovol, lposol, lvdwol

c     vlrc   = long-range LJ correction to energy, eV

      common/state/vlrc

      common/consts/pi

      common/lists_pots/inrlst(npmax), irblst(npmax), 
     .     noa(ntypes)

c     torsional potential parameters
c     V  = tor0 + tor10 * cos(phi)**10

      common/torsn_pot/tor0(ntypes,ntypes), tor10(ntypes,ntypes)

      common/xref_pot/ijty(ntypes,ntypes)

c     cross-reference variables

c     iat2ty(a) = atom type of atom a
c     kt(n)     = atom type for atomic number n
c     kt2(t)    = atomic number for atom type t
c     lrebot(t) = whether or not atom type t uses REBO
c     pad       = padding used to align variables on 8-byte boundaries (4-byte
c                 boundaries may be good enough?) and silence a compiler
c                 warning; size will change when ntypes changes
c     pad2      = another padding variable, used to make common block size
c                 a multiple of 16 bytes after rounding up to a multiple of
c                 4 bytes, in order to silence a gfortran compiler warning;
c                 size will change when ntypes changes

      integer iat2ty, kt, kt2
      character pad*4, pad2*4
      logical lrebot

      common/xref_pots/iat2ty(npmax),
     .     kt(mxatno), kt2(ntypes),
     .     pad,
     .     lrebot(ntypes),
     .     pad2

c    ZBL universal screening potential

      common/zbl/zbla(4), zblb(4)

c     Change these parameters to alter the number of atoms that can
c     be simulated.  Suggested values work well together.
c     nmabig = maximum # of pairs on LJ safe list (525k or 1875k or 4e6)

      parameter (nmabig = 1875000)

c     n2lmax = maximum # pairs in double-length REBO neighbor list
c     ntab   = array size for potential table look-up
c     nwlmax = maximum # of pairs on LJ watch list
c              (no greater than 4/3 pi (dijmax(C,C)+rbuffr)**3 * dens * N / 2
c              where dens = max. number density.  effectively <= 150 npmax)

      parameter (n2lmax = 2*nlmax)
      parameter (ntab = 10000)
      parameter (nwlmax = 500000)

c     # of bins in LJ lookup table

      parameter (maxbin = 13000)

c     size of xh arrays

      parameter (isizxh = ntypes * 10 * 10)

c     atom pair types

      parameter (ijcc   = 1)
      parameter (ijch   = 2)
      parameter (ijhh   = 3)
      parameter (ijcf   = 4)
      parameter (ijff   = 5)
      parameter (ijhf   = 6)
      parameter (mxprty = 6)

c     atom pair categories
c     hhc = hydro/halocarbon

      parameter (ihhc   = 1)
      parameter (isige  = 2)
      parameter (imixty = 3)

c     pair lists

      parameter (iljsaf = 2)
	    parameter (iljwat = 3)

c     switching region in coordination number-space for switch between
c     G_C and gamma_C

      parameter (bgqmin = 3.2d0)
      parameter (bgqmax = 3.7d0)

c     switching region in coordination number-space for conjugation
c     of carbons

      parameter (fmn = 2.d0)
      parameter (fmx = 3.d0)

c     flags:

c     REBO potential values:
c     hhtors = torsional barrier in ethane = 2882 cal/mol.
c     chtors = fictional torsional barrier for H3C-C(CR)3, chosen to
c              reproduce barrier to dihedral rotation in propane
c     cctors = fictional torsional barrier for (CR)3C-C(CR)3, chosen
c              to reproduce gauche/trans energy difference in butane

      parameter (hhtors = 2882.d0 * f2kilo * fkc2eV)

      parameter (chtors = 4120.d0 * f2kilo * fkc2eV)
      parameter (cctors = 7100.d0 * f2kilo * fkc2eV)

c     parameters for ZBL universal screening potential
c     (Ziegler, Biersack and Littmark, Stopping and Range of Ions in
c     Matter, Pergamon, 1977):

      parameter (zblcof = 0.8854d0)
      parameter (zblexp = 0.23d0)

c     parameters for J(r) integration
c     rjkmax = upper limit of integral, Angstroms
c     njk    = number of points in integral

      parameter (rjkmax = 100)
      parameter (njk    = 500)

c     parameters for J(r) lookup table
c     njtab = number of points in lookup table

      parameter (njtab  = 10000)

c     fundamental constants:

c     Bohr radius, '97 IUPAC, Angstroms
      parameter (bohr = 0.529177249d0)

c     speed of light, Angstroms/fs)
      parameter (clight = 299792458.d0 * fm2A / f2fem)

c     prefactor for electrostatic interacions,
c     1/(4 pi eps_0), in eV Angstrom/e^2, where eps_0 is the permittivity
c     of free space. The MKS value is c^2 * 1e-7 N/A^2
      parameter (epsinv = 1.d-7 * clight * clight * fJ2eV / fC2e / fC2e
     .                    * f2fem * f2fem / fm2A)

c     unit conversions...

c     ...from Angstroms to atomic units (Bohr)
      parameter (fA2au  = 1.d0 / bohr)

c     ...from cm to Anstroms
      parameter (fcm2A = 1.d8)

c     ...from erg to J
      parameter (ferg2J = 1.d-7)

c     Ewald parameters:

c     maximum number of lattice vectors in one Cartesian direction
      parameter (maxk1d = 15)
      
c     maximum number of k vectors
      parameter (maxkvc = 2000)

c     tolerance for error in Coulomb interaction in lookup tables (eV/e^2)
      parameter (eserrt = 1.d-6)

c     dtttol = tolerance for "zero" value of dS/dt for the dihedral angle
c     switching funcction
      parameter (dtttol = 1.d-8)

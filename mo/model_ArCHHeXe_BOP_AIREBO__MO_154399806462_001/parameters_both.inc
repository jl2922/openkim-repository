      include 'parameters_bothkim.inc'

c     maxmol = naximum # of molecules
      parameter (maxmol = npmax)

c     nqmat = maximum # of charges that can be optimized by matrix methods

      parameter (nqmat = 2500)

c     unit numbers:

      parameter (itempf = 99)

c     sizes of things:

c     kbyte  = number of bytes in a kilobyte (1024)
c     mbyte  = number of bytes in a megabyte
c     mxatno = largest atomic number

      parameter (maxcod = 10)
      parameter (mxatno = 100)
      parameter (kbyte  = 2 ** 10)
      parameter (mbyte  = kbyte * kbyte)

c     number of types of bond charges (sigma, pi)

      parameter (nbctyp = 2)

c     fundamental constants:

c     Avogadro's number, mol^-1
      parameter (avonum = 6.02214d23)

c     gas constant, R, J / mol / K
      parameter (gascon = 8.3145d0)

c     Boltzmann's constant, J / K
      parameter (bolzsi = gascon / avonum)

c     Planck's constant, J s
      parameter (plnksi = 6.626d-34)

c     unit conversions...

c     ...from X to gigaX
      parameter (f2giga = 1.d-9)

c     ...from X to kiloX
      parameter (f2kilo = 1.d-3)
c     ...and back
      parameter (fkilo2 = 1.d3)

c     ...from X to nanoX
      parameter (f2nano = 1.d9)

c     ...from X to picoX
      parameter (f2pico = 1.d12)
c     ...and back
      parameter (fpico2 = 1.d-12)

c     ...from X to femtoX
      parameter (f2fem = 1.d15)
c     ...and back
      parameter (ffem2 = 1.d-15)

c     ...from picoX to femtoX
      parameter (fp2f = fpico2 * f2fem)

c     ...from meters to Angstroms
      parameter (fm2A = 1.d10)
c     ...and back
      parameter (fA2m = 1.d-10)
      parameter (fA2cm = 1.d-8)

c     ...from electron charge to Coulomb
      parameter (fe2C = 1.602177d-19)
c     ...and back
      parameter (fC2e = 1.d0 / fe2C)

c     ...from electron charge to esu
      parameter (fe2esu = 4.803250d-10)

c     ...from esu-Angstroms to Debye
      parameter (fesA2D = 1.d10)

c     ...from electron-Angstroms to Debye
      parameter (feA2D = fe2esu * fesA2D)

c     ...from amu to grams
      parameter (famu2g = 1.d0 / avonum)
c     ... and back
      parameter (fg2amu = avonum)

c     ...from eV/A to nanoNewtons (force units)
      parameter (feV2nN = fe2C * fm2A * f2nano)
c     ...and back
      parameter (fnN2eV = 1.d0 / feV2nN)

c     ...from cal to J
      parameter (fcal2J = 4.184d0)

c     ...from eV to J
      parameter (feV2J = fe2C)
c     ...and back
      parameter (fJ2eV = 1.d0 / feV2J)

c     ...from J to erg
      parameter (fJ2erg = 1.d7)

c     ...from kcal/mol to eV
      parameter (fkc2eV = fkilo2 / avonum * fcal2J * fJ2eV)

c     ...from program energy units (amu A^2 / fs^2) to eV
      parameter (fE2eV = famu2g * f2kilo * fA2m * fA2m * f2fem * f2fem
     .     * fJ2eV)
c     ...and back
      parameter (feV2E = 1.d0 / fE2eV)

c     ...from J to K
      parameter (fJ2K = 1.d0 / bolzsi)

c     ...from eV to K
      parameter (feV2K = feV2J * fJ2K)
c     ...and back
      parameter (fK2eV = 1.d0 / feV2K)

c     Lorentz-Berthelot combining rules:

      logical lorber
      parameter (lorber = .true.)


c-----------------------------------------------------------------------
c     stvlrc is a wrapper that sets the long-range correction to the
c     potential energy on the potential side. written by yl.
c-----------------------------------------------------------------------
c
      subroutine stvlrc(vlrc1)

      include 'common_files.inc' 
      include 'common_pots.inc' 
	 
      vlrc = vlrc1

      end
c
c-----------------------------------------------------------------------
c     addtyp updates some variables and lists to know about the type
c     of a newly added atom....sjs 10/14/04
c-----------------------------------------------------------------------
c
      subroutine addtyp(i,iatnum)

      include 'common_files.inc'
      include 'common_pots.inc'

c     i      = number of new atom
c     iatnum = atomic number of new atom

      integer i, iatnum

c     keep track of the internal atom type

      iat2ty(i) = kt(iatnum)

c     complain if we've never heard of this atom

      if (iat2ty(i) .eq. 0) then
         write(isterr, *) 'addtyp: unknown atom type for atom ', i
         call ioclos
         stop
      endif

c     update the count of how many atoms of each type we have
      
      noa(iat2ty(i)) = noa(iat2ty(i)) + 1

c     update the lists of atoms with REBO type potentials and those
c     without

      if (lrebot(iat2ty(i))) then
         nrba = nrba + 1
         irblst(nrba) = i
      else
         nnra = nnra + 1
         inrlst(nnra) = i
      endif

c     mark the pair lists as expired

      lcovol = .true.
      lvdwol = .true.

c     note that some atoms have moved

      return
      end
c
c-----------------------------------------------------------------------
c     deltyp updates some variables and lists to know that an atom has
c     recently been deleted....sjs 10/15/04
c-----------------------------------------------------------------------
c
      subroutine deltyp(iatom,np)

      include 'common_files.inc'
      include 'common_pots.inc'

c     iatom = number of deleted atom
c     np    = number of atoms

      integer iatom, np

c     ... and then there were N-1

      noa(iat2ty(iatom)) = noa(iat2ty(iatom)) - 1
      
c     update the lists of atoms with REBO type potentials and those
c     without. make sure to use the new atom numbers.

      if (lrebot(iat2ty(iatom))) then
         nrba = nrba - 1
         do 10 irba = 1, nrba
            if (irblst(irba) .ge. iatom) then
               irblst(irba) = irblst(irba+1)
               if (irblst(irba) .gt. iatom) then
                  irblst(irba) = irblst(irba) - 1
               endif
            endif
 10      continue
         do 15 inra = 1, nnra
            if (inrlst(inra) .gt. iatom) then
               inrlst(inra) = inrlst(inra) - 1
            endif
 15      continue
      else
         nnra = nnra - 1
         do 20 inra = 1, nnra
            if (inrlst(inra) .ge. iatom) then
               inrlst(inra) = inrlst(inra+1)
               if (inrlst(inra) .gt. iatom) then
                  inrlst(inra) = inrlst(inra) - 1
               endif
            endif
 20      continue
         do 25 irba = 1, nrba
            if (irblst(irba) .gt. iatom) then
               irblst(irba) = irblst(irba) - 1
            endif
 25      continue
      endif

c     update the list of atom types

      do 100 jatom = iatom, np
         iat2ty(jatom) = iat2ty(jatom+1)
 100  continue

c     mark the pair lists as expired

      lcovol = .true.
      lvdwol = .true.

      return
      end
c
c-----------------------------------------------------------------------
c     expmol expires the molecular membership list...sjs 11/25/09
c-----------------------------------------------------------------------
c
      subroutine expmol

      include 'common_files.inc' 
      include 'common_pots.inc' 

      lmolol = .true.

      return
      end
c
c-----------------------------------------------------------------------
c     exppos marks the atom positions as having moved since the last
c     time the pair list was checked....sjs 12/17/09
c-----------------------------------------------------------------------
c
      subroutine exppos

      include 'common_files.inc'
      include 'common_pots.inc'

      lposol = .true.
      lmolol = .true.

      return
      end
c

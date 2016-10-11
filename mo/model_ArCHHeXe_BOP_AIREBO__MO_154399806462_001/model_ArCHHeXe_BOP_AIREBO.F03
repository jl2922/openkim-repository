!
! CDDL HEADER START
!
! The contents of this file are subject to the terms of the Common Development
! and Distribution License Version 1.0 (the "License").
!
! You can obtain a copy of the license at
! http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
! specific language governing permissions and limitations under the License.
!
! When distributing Covered Code, include this CDDL HEADER in each file and
! include the License file in a prominent location with the name LICENSE.CDDL.
! If applicable, add the following below this CDDL HEADER, with the fields
! enclosed by brackets "[]" replaced with your own identifying information:
!
! Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
!
! CDDL HEADER END
!

!
! Copyright (c) 2012, Clemson University.  All rights reserved.
!
! Contributors:
!    Steven J. Stuart
!    Ellad B. Tadmor
!
!-------------------------------------------------------------------------------
! 
! Module model_ArCHHeXe_BOP_AIREBO
!
! Compute energy and forces on a hydrocarbon system using the AIREBO
! potential.[1]
!
! [1] S.J. Stuart, A.B. Tutein and J.A. Harrison, "A Reactive Potential for
!     Hydrocarbons with Intermolecular Interactions", J. Chem. Phys. 112,
!     6472-6486 (2000).
!
! Language: Fortran 2003
!
!-------------------------------------------------------------------------------

#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))

module model_ArCHHeXe_BOP_AIREBO

  use, intrinsic :: iso_c_binding
  use KIM_API_F03
  implicit none

  save
  private
  public Compute_Energy_Forces
  public model_cutoff
  public Destroy
  public cd

  ! parameters

  integer(c_int), parameter :: cd = c_double ! used for literal constants

  real(c_double), parameter :: model_cutoff = 10.4_cd ! C-C cutoff [A]
  ! neighbor list styles
  integer(c_int), parameter          :: inone = 0
  integer(c_int), parameter          :: ihalf = 1
  integer(c_int), parameter          :: ifull = 2
  ! iterator styles for neighbor list
  integer(c_int), parameter          :: iterator = 1
  integer(c_int), parameter          :: ilocator = 2
  ! another definition of iterator modes for neighbor list
  integer(c_int), parameter          :: iteratorMode = 0
  integer(c_int), parameter          :: locatorMode  = 1
  ! requests for iterator-based neighbor lists
  integer(c_int), parameter          :: neighReset = 0
  integer(c_int), parameter          :: neighIncrement = 1

contains

!-------------------------------------------------------------------------------
! Compute energy and forces for a configuration
!-------------------------------------------------------------------------------

  integer(c_int) function Compute_Energy_Forces(pkim) bind(c)
    implicit none
    include 'genpar.inc'
    include 'parameters_bothkim.inc'

    ! Arguments
    type(c_ptr), intent(in)  :: pkim

    ! Local variables
    integer(c_int) :: ier

    ! maximum number of neighbors of an atom, using the published cutoff
    integer(c_int), parameter :: prkim = 100

    character(len=80) :: errormsg
    character(len=25) :: modelname = 'model_ArCHHeXe_BOP_AIREBO'
    integer(c_int) :: pairs(nlmax,2)
    integer(c_int) :: neighborsOf(npmax,prkim)
    integer(c_int) :: nneighbors(npmax)
    integer(c_int) :: atnum(ntypes)
    integer(c_int) :: atom2, comp_energy, comp_force, &
               i, iatom, idum, ineigh, j, jatom, jindex, &
               NLaccess, NLstyle, nneigh, npairs
    logical :: useNL

    ! AIREBO variables
    integer(c_int) iatno(npmax)
    integer(c_int) atomof(nbcmax,2)
    integer(c_int) igfunc, ipot, nbc, isperiodic
    logical lpbc(ndim)
    logical lbc, les, lewald, lewsrf, llj, lnones, lqsolv, ltors, lvor
    real(c_double) fext(npmax,ndim), fint(npmax,ndim), r0(npmax,ndim)
    real(c_double) chi(npmax), q0(npmax)
    real(c_double) bcchi(nbcmax), bctype(nbcmax), bq0(nbcmax)
    real(c_double) uu(ndim,ndim)
    real(c_double) cube(ndim), Efield(ndim)
    real(c_double) ewlkpl, ewkmxw, tote

    ! KIM variables
    integer(c_int), pointer :: particleSpecies(:); type(c_ptr) :: pparticleSpecies ! atom types
    real(c_double), pointer :: boxlen(:);          type(c_ptr) :: pboxlen          ! box size
    real(c_double), pointer :: coord(:,:);         type(c_ptr) :: pcoord           ! coordinates
    real(c_double), pointer :: Rij_dummy(:,:);     type(c_ptr) :: pRij_dummy       ! dummy
    real(c_double), pointer :: energy;             type(c_ptr) :: penergy          ! PE
    real(c_double), pointer :: force(:,:);         type(c_ptr) :: pforce           ! forces
    integer(c_int), pointer :: nonGhost;           type(c_ptr) :: pnonGhost        ! # non-ghost particles
    integer(c_int), pointer :: np;                 type(c_ptr) :: pnp              ! # particles
    integer(c_int), pointer :: nASpecs;            type(c_ptr) :: pnASpecs         ! # atom species
    character(len=KIM_KEY_STRING_LENGTH) NBC_Method                                ! boundary conditions
    integer(c_int), pointer :: neighbors(:);       type(c_ptr) :: pneighbors       ! neighbor list

    isperiodic = 0

    ! Check to determine boundary conditions and neighbor list model

    ier = kim_api_get_nbc_method(pkim, NBC_Method)
    if (ier .lt. KIM_STATUS_OK) then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_get_nbc_method_f", ier)
       goto 42
    endif

    if (index(NBC_Method, "CLUSTER") .eq. 1) then 
       useNL = .false.       ! no neighbor list is available
       NLstyle = inone
       NLaccess = inone
       do i = 1, ndim
          lpbc(i) = .false.  ! no periodic boundary conditions should be applied
       end do
    else
       ! check the neighbor list access mode
       NLaccess = kim_api_get_neigh_mode(pkim, ier)
       if (ier .lt. KIM_STATUS_OK) then
          idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                      "kim_api_get_neigh_mode", ier)
          goto 42
       endif
       if (NLaccess .eq. 3) then  ! both access modes are supported by Test...
          NLaccess = ilocator     ! ... so choose easier Locator mode
       endif
       if (NLaccess .ne. iterator .and. NLaccess .ne. ilocator) then
          ier = KIM_STATUS_FAIL
          write(errormsg, *) "unknown neighbor list access mode ", NLaccess
          idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                      errormsg, ier)
          goto 42
       endif
       if (index(NBC_Method, "MI_OPBC_F") .eq. 1) then
          useNL = .true.        ! neighbor list is available
          NLstyle = ifull       ! full (redundant) neighbor list
          do i = 1, ndim
             lpbc(i) = .true.   ! periodic boundary conditions must be applied
          end do
          isperiodic = 1
       else if (index(NBC_Method, "MI_OPBC_H") .eq. 1) then
          useNL = .true.        ! neighbor list is available
          NLstyle = ihalf       ! half (non-redundant) neighbor list
          do i = 1, ndim
             lpbc(i) = .true.   ! periodic boundary conditions must be applied
          end do
          isperiodic = 1
       else if (index(NBC_Method, "NEIGH_PURE_F") .eq. 1) then
          useNL = .true.        ! neighbor list is available
          NLstyle = ifull       ! full (redundant) neighbor list
          do i = 1, ndim
             lpbc(i) = .false.  ! periodic boundary conditions should not be applied
          end do
       else if (index(NBC_Method, "NEIGH_PURE_H") .eq. 1) then
          useNL = .true.        ! neighbor list is available
          NLstyle = ihalf       ! half (redundant) neighbor list
          do i = 1, ndim
             lpbc(i) = .false.  ! periodic boundary conditions should not be applied
          end do
       else
          write(errormsg, *) "cannot support neighbor list method ", &
               NBC_Method(1:12)
          ier = KIM_STATUS_FAIL
          idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                      errormsg, ier)
          goto 42
       endif
    endif

    ! Check to see if we need energy and/or forces

    call kim_api_getm_compute(pkim, ier, &
                              "energy", comp_energy, 1, &
                              "forces", comp_force,  1)

    if (ier .lt. KIM_STATUS_OK) then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_getm_compute", &
                                   ier)
       goto 42
    endif

    ! Unpack data from KIM object

    call kim_api_getm_data(pkim, ier, &
         "numberOfParticles", pnp,             1,                           &
         "numberOfSpecies",   pnASpecs,        1,                           &
         "particleSpecies",   pparticleSpecies,1,                           &
         "coordinates",       pcoord,          1,                           &
         "boxSideLengths",    pboxlen,         TRUEFALSE(isperiodic.eq.1),  &
         "energy",            penergy,         TRUEFALSE(comp_energy.eq.1), &
         "forces",            pforce,          TRUEFALSE(comp_force.eq.1))
    if (ier .lt. KIM_STATUS_OK) then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_getm_data", ier)
       goto 42
    endif

    call c_f_pointer(pnp,              np)
    call c_f_pointer(pnASpecs,         nASpecs)
    call c_f_pointer(pparticleSpecies, particleSpecies, [np])
    call c_f_pointer(pcoord,           coord,           [ndim,np])
    if (isperiodic.eq.1)  call c_f_pointer(pboxlen, boxlen, [ndim])
    if (comp_energy.eq.1) call c_f_pointer(penergy, energy)
    if (comp_force.eq.1)  call c_f_pointer(pforce,  force, [ndim,np])

    ! Make sure the (statically allocated) AIREBO code can handle this many atoms
    if (np .gt. npmax) then
       write(isterr, *) modelname, ': numberOfParticles = ', np, &
            ' greater than npmax = ', npmax
       write(isterr, *) '     increase npmax and recompile'
       ier = KIM_STATUS_FAIL
       goto 42
    endif

    ! Check how many ghost atoms we have

    ! repack coordinates into form needed by AIREBO routines
    do i = 1, np
       do j = 1, ndim
          r0(i,j) = coord(j,i)
       end do
    end do

    ! repack box side lengths into form needed by AIREBO routines
    if (isperiodic .eq. 1) then
       do i = 1, ndim
          cube(i) = boxlen(i)
       end do
    endif

    ! kludge
    atnum(ihyd) = 1
    atnum(ihel) = 2
    atnum(icarb) = 6
    atnum(iargon) = 18
    atnum(ixenon) = 54
    ! repack atom types into AIREBO array
    do i = 1, np
       iatno(i) = atnum(particleSpecies(i))
       call addtyp(i, iatno(i))
    end do

    ! we don't have charges; set them to zero
    do i = 1, np
       q0(i) = 0.0_cd
    end do

    ! no bond charges, either
    nbc = 0
    do i = 1, nbc
       bq0(i) = 0.0_cd
       atomof(i,1) = 0
       atomof(i,2) = 0
       bctype(i) = isigma
    end do

    ! no external electric field
    do i = 1, ndim
       Efield(i) = 0.0_cd
    end do

    ! we can only handle cluster (non-periodic) boundary conditions for now, so
    ! just define a big box
    do i = 1, ndim
       cube(i) = 9.e24_cd
       lpbc(i) = .false.
    end do

    ! use a REBO-style potential
    ipot = irebo

    ! use spline-based angular functions
    igfunc = isplin

    ! use Lennard-Jones
    llj = .true.

    ! use torsions
    ltors = .true.

    ! no electrostatics
    les = .false.
    lqsolv = .false.
    lnones = .true.
    lewald = .false.
    lewsrf = .false.
    ewlkpl = 6.0_cd  ! value is irrelevant
    ewkmxw = 1.6_cd  ! value is irrelevant

    ! don't calculate Voronoi volumes
    lvor = .false.

    if (useNL) then

       ! set up some AIREBO state variables

       call bldstt(np, iatno, cube, lpbc, lewald, ewlkpl, ewkmxw)

       if (NLstyle .eq. ihalf) then
          
          ! check whether there are any ghost atoms

          pnonGhost = kim_api_get_data(pkim, "numberContributingParticles",ier)
          if (ier .lt. KIM_STATUS_OK) then
             write(errormsg, *) "kim_api_get_data failed"
             idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         errormsg, ier)
             goto 42
          endif
          call c_f_pointer(pnonGhost, nonGhost)

          if (nonGhost .ne. np) then
             ier = KIM_STATUS_FAIL
             write(errormsg, *) "Can't handle ghost atoms yet"
             idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         errormsg, ier)
             goto 42
          endif
          
       endif

       ! get the neighbor list info from the KIM API

       if (NLaccess .eq. iterator) then

          ! store neighbors in an array where they can be accessed directly
          
          do iatom = 1, np
             nneighbors(iatom) = 0
          end do

          ier = kim_api_get_neigh(pkim, iteratorMode, neighReset, atom2, &
                                  nneigh, pneighbors, pRij_dummy)
          if (ier .ne. KIM_STATUS_NEIGH_ITER_INIT_OK) then
             write(errormsg, *) "kim_api_get_neigh failed"
             idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                         errormsg, ier)
             goto 42
          endif

          do 
             ier = kim_api_get_neigh(pkim, iteratorMode, neighIncrement, &
                                     atom2, nneigh, pneighbors, pRij_dummy)
             if (ier .eq. KIM_STATUS_NEIGH_ITER_PAST_END) then
                exit
             else if (ier .ne. KIM_STATUS_OK) then
                write(errormsg, *) "kim_api_get_neigh failed"
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                            errormsg, ier)
                goto 42
             endif

             iatom = atom2
             nneighbors(iatom) = nneigh
             if (nneighbors(iatom) .gt. prkim) then
                write(errormsg, *) "atom ", iatom, " has ", nneighbors(iatom), &
                     " neighbors; only prkim = ", prkim, " allowed"
                ier = KIM_STATUS_FAIL
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                            errormsg, ier)
                goto 42
             endif

             call c_f_pointer(pneighbors, neighbors, [nneigh])

             do ineigh = 1, nneighbors(iatom)
                neighborsOf(iatom,ineigh) = neighbors(ineigh)
             end do

          end do

          ! now repack them into a non-redundant, ordered pair list

          npairs = 0
          do iatom = 1, np
             do ineigh = 1, nneighbors(iatom)
                if (NLstyle .eq. ihalf &
                     .or. (NLstyle .eq. ifull .and. iatom .le. jatom)) then
                   npairs = npairs + 1
                   pairs(npairs,1) = iatom
                   pairs(npairs,2) = neighborsOf(iatom,ineigh)
                endif
             end do
          end do

       else

          npairs = 0
          do iatom = 1, np
             ier = kim_api_get_neigh(pkim, locatorMode, iatom, atom2, nneigh, &
                                     pneighbors, pRij_dummy)
             if (ier .lt. KIM_STATUS_OK) then
                write(errormsg, *) "kim_api_get_neigh failed"
                idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                            errormsg, ier)
                goto 42
             endif

             call c_f_pointer(pneighbors, neighbors, [nneigh])

                ! repack the non-redundant pair list provided by KIM into 
                ! an array

                do jindex = 1, nneigh
                   jatom = neighbors(jindex)
                   if (NLstyle .eq. ihalf &
                        .or. (NLstyle .eq. ifull .and. iatom .le. jatom)) then
                      npairs = npairs + 1
                      pairs(npairs,1) = iatom
                      pairs(npairs,2) = jatom
                   endif
                end do
                
          end do

          ! send the non-redundant pair list to AIREBO for use in generating
          ! several AIREBO pair lists

          call pakpls(npairs, pairs, np, cube, r0, lpbc)

       endif

    else
       ! no pair list is available; build it ourself

       call build(np, r0, iatno, cube, lpbc, lewald, ewlkpl, ewkmxw, lvor)
    endif

    ! for charge-based calculations we would need to initialize charges next
    ! (see energymain.f) but that is not needed here

    ! Compute energy and forces

    call calcforce(np, r0, q0, lbc, nbc, bq0, atomof, bctype, Efield, lpbc, &
                   cube, ipot, igfunc, llj, ltors, les, lqsolv, lnones, &
                   lewald, lewsrf, fint, chi, bcchi, uu, tote)

    if (comp_energy .eq. 1) then
       energy = tote
    endif
    if (comp_force .eq. 1) then
       do i = 1, np
          do j = 1, ndim
             force(j,i) = fint(i,j)
          end do
       end do
    endif
    
    ! Finished happily
    ier = KIM_STATUS_OK

42  continue
    Compute_Energy_Forces = ier
    return
  end function Compute_Energy_Forces

  !-------------------------------------------------------------------------------
  ! Model destroy routine
  !-------------------------------------------------------------------------------

  integer(c_int) function Destroy(pkim) bind(c)

    implicit none

    ! Arguments
    type(c_ptr), intent(in) :: pkim

    ! No published parameters, so nothing to deallocate

    Destroy = KIM_STATUS_OK
    return

  end function Destroy

end module model_ArCHHeXe_BOP_AIREBO

!-------------------------------------------------------------------------------
! Model initialization routine (REQUIRED)
!-------------------------------------------------------------------------------

integer(c_int) function model_init(pkim) bind(c)

  use, intrinsic :: iso_c_binding
  use model_ArCHHeXe_BOP_AIREBO
  use KIM_API_F03

  implicit none
  include 'parameters_bothkim.inc'

  ! Argument
  type(c_ptr), intent(in) :: pkim

  ! KIM variables

  real(c_double), pointer :: cutoff; type(c_ptr) :: pcutoff

  ! Local variables
  integer(c_int), parameter :: one=1
  integer(c_int) ityc3p(maxc3t,3)
  integer(c_int) ityc2p(maxc2t,2)
  integer(c_int) ityc1p(maxc1t,1)
  integer(c_int) idum, ier, igfunc, ipot, iverb, ncpair, nctrip, nctype
!  integer(c_int), allocatable, target :: ityc3p(:,:), itypc2p(:,:), ityc1p(:)
  logical usec3p(maxc3t,maxc3p)
  logical usec2p(maxc2t,maxc2p)
  logical usec1p(maxc1t,maxc1p)
  logical lbc, lckfil, les, ljdir, llj, lljdir, lpdir, lvor
  real(c_double) valc3p(maxc3t,maxc3p)
  real(c_double) valc2p(maxc2t,maxc2p)
  real(c_double) valc1p(maxc1t,maxc1p)
  real(c_double) rbuffr


  ! store pointer to compute function in KIM object
  ier = kim_api_set_method(pkim, "compute", one, c_funloc(Compute_Energy_Forces))
  if (ier .lt. KIM_STATUS_OK) then
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_set_method", ier)
     goto 42
  endif

  ! store pointer to destroy subroutine in KIM object
  ier = kim_api_set_method(pkim, "destroy", one, c_funloc(Destroy))
  if (ier .lt. KIM_STATUS_OK) then
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_set_method", ier)
     goto 42
  endif
  
  ! store model cutoff in KIM object
  pcutoff = kim_api_get_data(pkim, "cutoff", ier)
  if (ier .lt. KIM_STATUS_OK) then
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                 "kim_api_get_data", ier)
     goto 42
  endif
  call c_f_pointer(pcutoff, cutoff)
  cutoff = model_cutoff
  
  ! initialize the AIREBO model

  ! this initializes everything that must be initialized before the system state
  ! is known
  nctype = 0
  ncpair = 0
  nctrip = 0
  lckfil = .false. ! parameters will come from KIM, not files
  call preint(nctype, ityc1p, usec1p, valc1p, &
              ncpair, ityc2p, usec2p, valc2p, &
              nctrip, ityc3p, usec3p, valc3p, lckfil)

  ! this initializes everything that depends on run conditions (irrelevant for
  ! for KIM)

  igfunc = isplin ! use splines for the angular function (defaut AIREBO)
  call runint(igfunc)

  ! assume direct evaluation of pair functions for now,
  ! because (1) we don't yet know how many atoms the system will have or how
  ! many evaluations are needed, and (2) it is more conservative and (3) it
  ! avoids having to think about how to store the lookup table in the KIM API

  ljdir = .true.  ! direct evaluation of Coulomb interactions
  lpdir = .true.  ! direct evaluation of pair functions
  lljdir = .true. ! direct evaluation of Lennard-Jones interactions

  ! AIREBO class of potential:
  lbc = .false.  ! no bond charges
  ipot = irebo   ! REBO family of potential
  llj = .true.   ! use LJ
  les = .false.  ! no elecrostatics

  ! some flags that affect what gets reported:
  rbuffr = 0.0_cd  ! buffer (skin) width for pair list [A] (test handles pair
                   ! list, so we don't need a skin on pair list)
  lvor = .false. ! don't calculate Voronoi volumes
  iverb = 3      ! verbosity (0 = silent, 9 = debug)

  ! this is where most of the potential parameters are initialized
  call init(lbc, ipot, llj, rbuffr, les, lvor, ljdir, lpdir, lljdir, iverb)

  ier = KIM_STATUS_OK
42 continue
  model_init = ier
  return
  
end function model_init


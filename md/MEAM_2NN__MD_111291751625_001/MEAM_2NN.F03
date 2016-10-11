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
! Copyright (c) 2012--2014, Regents of the University of Minnesota.
! All rights reserved.
!
! Contributors:
!    Ryan S. Elliott
!    Ellad B. Tadmor
!    Valeriu Smirichinski
!    Feng Gao
!    Eunkoo Lee
!    Bohumir Jelinek
!

!****************************************************************************
!**
!**  MODULE kcc_meam
!**
!**  kcc_meam pair potential model for <FILL element name>
!**
!**  Reference: <FILL>
!**
!**  Language: Fortran 90
!**
!**  Release: 
!**
!****************************************************************************


#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))

!
! All MEAM parameters stored in this module
!
module meam_parameters
use, intrinsic :: iso_c_binding
implicit none

integer(c_int), parameter :: cd = c_double  ! used for literal constants

save
public

type meam_params

   integer(c_int) :: max_comp
   integer(c_int) :: LineTable
   ! lj_1
   real(c_double), pointer :: PhiTab(:,:,:)
   real(c_double), pointer :: DPhiTab(:,:,:)
   ! lj_2
   real(c_double)          :: Rmin
   real(c_double)          :: Rsqmin
   real(c_double)          :: Rcutoff
   real(c_double)          :: OverDeltaRsq
   ! eam1
   real(c_double), pointer :: Ec(:,:)
   real(c_double), pointer :: Re(:)
   real(c_double), pointer :: B(:,:)
   real(c_double), pointer :: A(:)
   real(c_double), pointer :: beta(:,:)
   real(c_double), pointer :: t(:,:)
   ! eam2
   real(c_double), pointer :: RhoZ(:)
   real(c_double), pointer :: Rho0_Bar(:)
   real(c_double), pointer :: Z(:)
   real(c_double), pointer :: Z2(:)
   real(c_double), pointer :: A2(:)
   ! eam3
   real(c_double), pointer :: Cmin(:,:,:)
   real(c_double), pointer :: Cmax(:,:,:)
   real(c_double)          :: Dcutoff
   real(c_double)          :: OverDeltaR
   real(c_double), pointer :: Eu_d(:)
   ! eam4
   real(c_double), pointer :: Re_Alloy(:,:)
   real(c_double), pointer :: alpha(:,:)
   real(c_double), pointer :: omega(:,:)
   ! ener (only amass)
   real(c_double), pointer :: Amass(:)
   ! scrn
   real(c_double)          :: scr_2nn
   !
   integer(c_int) :: Ncomp
   character(len=2), pointer :: Icomp(:)

end type meam_params

contains

   subroutine allocate_meam_parameters(mp,max_comp,LineTable)
   implicit none

!--Transferred variables
      type(meam_params), intent(inout) :: mp
      integer(c_int), intent(in) :: max_comp
      integer(c_int), intent(in) :: LineTable
      mp%max_comp  = max_comp
      mp%LineTable = LineTable
      allocate(mp%PhiTab(LineTable+1,max_comp,max_comp))
      allocate(mp%DPhiTab(LineTable+1,max_comp,max_comp))
      allocate(mp%Ec(max_comp,max_comp))
      allocate(mp%Re(max_comp))
      allocate(mp%B(max_comp,max_comp))
      allocate(mp%A(max_comp))
      allocate(mp%beta(max_comp,4))
      allocate(mp%t(max_comp,3))
      allocate(mp%RhoZ(max_comp))
      allocate(mp%Rho0_Bar(max_comp))
      allocate(mp%Z(max_comp))
      allocate(mp%Z2(max_comp))
      allocate(mp%A2(max_comp))
      allocate(mp%Cmin(max_comp,max_comp,max_comp))
      allocate(mp%Cmax(max_comp,max_comp,max_comp))
      allocate(mp%Eu_d(max_comp))
      allocate(mp%Re_Alloy(max_comp,max_comp))
      allocate(mp%alpha(max_comp,max_comp))
      allocate(mp%omega(max_comp,max_comp))
      allocate(mp%Amass(max_comp))
      allocate(mp%Icomp(max_comp))
    end subroutine allocate_meam_parameters

    subroutine deallocate_meam_parameters(mp)
      implicit none
      type(meam_params), intent(inout) :: mp
      deallocate(mp%PhiTab)
      deallocate(mp%DPhiTab)
      deallocate(mp%Ec)
      deallocate(mp%Re)
      deallocate(mp%B)
      deallocate(mp%A)
      deallocate(mp%beta)
      deallocate(mp%t)
      deallocate(mp%RhoZ)
      deallocate(mp%Rho0_Bar)
      deallocate(mp%Z)
      deallocate(mp%Z2)
      deallocate(mp%A2)
      deallocate(mp%Cmin)
      deallocate(mp%Cmax)
      deallocate(mp%Eu_d)
      deallocate(mp%Re_Alloy)
      deallocate(mp%alpha)
      deallocate(mp%omega)
      deallocate(mp%Amass)
      deallocate(mp%Icomp)
    end subroutine deallocate_meam_parameters

end module meam_parameters

module kcc_meam

use, intrinsic :: iso_c_binding
use KIM_API_F03
implicit none

save
private
public Compute_Energy_Forces, &
       ReInit,                &
       Destroy,               &
       DIM

! Below are the definitions and values of all Model parameters
integer(c_int), parameter          :: DIM=3          ! dimensionality of space

contains

!-------------------------------------------------------------------------------
!
! Compute energy and forces on atoms from the positions.
!
!-------------------------------------------------------------------------------
  integer(c_int) function Compute_Energy_Forces(pkim) bind(c)
    use, intrinsic :: iso_c_binding
    use meam_parameters
    implicit none

    !-- Transferred variables
    type(c_ptr), intent(in)  :: pkim

    !-- Local variables
    integer(c_int) ier
    integer(c_int) i,idum
    type(meam_params), pointer :: mp; type(c_ptr) :: a_mp

    !-- KIM variables
    integer(c_int), pointer :: numberOfParticles;  type(c_ptr) :: pnAtoms
    integer(c_int), pointer :: nparticleSpecies;   type(c_ptr) :: pnparticleSpecies
    integer(c_int), pointer :: particleSpecies(:); type(c_ptr) :: pptypes
    real(c_double), pointer :: cutoff;             type(c_ptr) :: pcutoff
    real(c_double), pointer :: energy;             type(c_ptr) :: penergy
    real(c_double), pointer :: coor(:,:);          type(c_ptr) :: pcoor
    real(c_double), pointer :: force(:,:);         type(c_ptr) :: pforce
    real(c_double), pointer :: ene_pot(:);         type(c_ptr) :: penepot
    real(c_double), pointer :: virial_global(:);   type(c_ptr) :: pvirial
    integer(c_int) comp_energy, comp_force, comp_enepot, comp_virial
    real(c_double), pointer :: boxSideLengths(:);  type(c_ptr) :: pboxSideLengths

    ! Get model buffer from KIM object
    !
    a_mp = kim_api_get_model_buffer(pkim, ier)
    if (ier.lt.KIM_STATUS_OK) then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_get_model_buffer", ier)
       goto 42
    endif
    call c_f_pointer(a_mp, mp)

    ! Check to see if we have been asked to compute the energy, forces, energyperatom,
    ! and virial
    !
    call kim_api_getm_compute(pkim, ier,   &
         "energy",         comp_energy, 1, &
         "forces",         comp_force,  1, &
         "particleEnergy", comp_enepot, 1, &
         "virial",         comp_virial, 1)
    if (ier.lt.KIM_STATUS_OK) then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_getm_compute", ier)
       goto 42
    endif
    if (comp_enepot.eq.1) then
       ier = KIM_STATUS_FAIL
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
         "ERROR: Energy per atom not supported by MEAM yet.", ier)
       goto 42
    endif

    ! Unpack data from KIM object
    !
    call kim_api_getm_data(pkim, ier,                                 &
         "numberOfParticles",           pnAtoms,         1,           &
         "numberOfSpecies",             pnparticleSpecies, 1,         &
         "particleSpecies",             pptypes,         1,           &
         "cutoff",                      pcutoff,         1,           &
         "coordinates",                 pcoor,           1,           &
         "boxSideLengths",              pboxSideLengths, 1,           &
         "energy",                      penergy,         comp_energy, &
         "forces",                      pforce,          comp_force,  &
         "particleEnergy",              penepot,         comp_enepot, &
         "virial",                      pvirial,         comp_virial)
    if (ier.lt.KIM_STATUS_OK) then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                   "kim_api_getm_data", ier)
       goto 42
    endif
    call c_f_pointer(pnAtoms, numberOfParticles)
    call c_f_pointer(pptypes, particleSpecies, [numberOfParticles])
    call c_f_pointer(pnparticleSpecies, nparticleSpecies)
    call c_f_pointer(pcutoff, cutoff)
    call c_f_pointer(pcoor,   coor,            [DIM,numberOfParticles])
    call c_f_pointer(pboxSideLengths, boxSideLengths, [DIM])
    if (comp_energy.eq.1) &
       call c_f_pointer(penergy, energy)
    if (comp_force.eq.1) &
       call c_f_pointer(pforce,  force,   [DIM,numberOfParticles])
    if (comp_enepot.eq.1) &
       call c_f_pointer(penepot, ene_pot, [numberOfParticles])
    if (comp_virial.eq.1) &
       call c_f_pointer(pvirial, virial_global, [6])

    ! Check that the number of atoms does not exceed the hard limit
    ! coded in MEAM.F
    if (numberOfParticles>20000) then
       ier = KIM_STATUS_FAIL ! assume an error
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
            "Number of atoms exceeds MEAM hard limit of 20000", ier)
       goto 42
    endif

    ! Check to be sure that the atom types are correct by comparing
    ! the provided species codes to the value given here (which should
    ! be the same as that given in the .kim file).
    !
    ier = KIM_STATUS_FAIL ! assume an error
    do i = 1,numberOfParticles
       if (particleSpecies(i) .gt. mp%max_comp) then
          idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
            "Wrong Atom Species", ier)
          goto 42
       endif
    enddo
    ier = KIM_STATUS_OK ! everything is ok

    ! Initialize potential energies, forces, virial term
    !
    if (comp_enepot.eq.1) ene_pot(1:numberOfParticles) = 0.0_cd
    if (comp_energy.eq.1) energy = 0.0_cd
    if (comp_force.eq.1)  force(1:3,1:numberOfParticles) = 0.0_cd
    if (comp_virial.eq.1) virial_global = 0.0_cd

    ! Call MEAM potential
    !
    call MEAM_potential(pkim,mp,numberOfParticles,coor,particleSpecies, &
            boxSideLengths,energy,ene_pot,force,virial_global,          &
            comp_energy,comp_enepot,comp_force,comp_virial,ier)

    ier = KIM_STATUS_OK
42  continue
    Compute_Energy_Forces = ier
    return

  end function Compute_Energy_Forces

!-------------------------------------------------------------------------------
!
! Model reinitialization routine
!
!-------------------------------------------------------------------------------
  integer(c_int) function ReInit(pkim) bind(c)
    implicit none

    !-- Transferred variables
    type(c_ptr), intent(in) :: pkim

    ! assume all is well
    !
    ReInit = KIM_STATUS_OK

  end function ReInit

!-------------------------------------------------------------------------------
!
! Model destroy routine
!
!-------------------------------------------------------------------------------
  integer(c_int) function Destroy(pkim) bind(c)
    use KIM_API_F03
    use meam_parameters
    implicit none

    !-- Transferred variables
    type(c_ptr), intent(in)  :: pkim

    !-- Local variables
    integer(c_int) ier
    integer(c_int) idum
    type(meam_params), pointer :: mp; type(c_ptr) :: a_mp

    ! get model buffer from KIM object
    a_mp = kim_api_get_model_buffer(pkim, ier)
    if (ier.lt.KIM_STATUS_OK) then
       idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                     "kim_api_get_model_buffer", ier)
       goto 42
    endif
    call c_f_pointer(a_mp, mp)

    ! free storage
    call deallocate_meam_parameters(mp)

    ier = KIM_STATUS_OK
42  continue
    Destroy = ier
    return

  end function Destroy

end module kcc_meam

!-------------------------------------------------------------------------------
!
! Model initialization routine (REQUIRED)
!
!-------------------------------------------------------------------------------
integer(c_int) function model_driver_init(pkim, pparamfile, nmstrlen, numparamfiles) bind(c)
use kcc_meam
use meam_parameters
use KIM_API_F03
implicit none

!-- Transferred variables
type(c_ptr),        intent(in) :: pkim
type(c_ptr), value, intent(in) :: pparamfile
integer(c_int),     intent(in) :: nmstrlen
integer(c_int),     intent(in) :: numparamfiles
character(len=nmstrlen), pointer :: paramfile(:)

!-- Local variables
integer(c_int), parameter :: one=1
integer(c_int) ier, idum, return_error
integer(c_int) i, j, strlen

!-- KIM variables
real(c_double), pointer :: cutoff; type(c_ptr) :: pcutoff
character(len=KIM_KEY_STRING_LENGTH) species;
integer(c_int) num_species
type(meam_params), pointer :: mp

! assume all is well
!
return_error = KIM_STATUS_OK

call c_f_pointer(pparamfile, paramfile, [numparamfiles])

! find first null character and write spaces everywhere afterward
do i = 1, numparamfiles
  j = index(paramfile(i),char(0))
  paramfile(i)(j:)=" "
end do

! store pointers in KIM object
!
call kim_api_setm_method(pkim, ier, &
     "compute", one, c_funloc(Compute_Energy_Forces), 1, &
     "reinit",  one, c_funloc(ReInit),                1, &
     "destroy", one, c_funloc(Destroy),               1)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_setm_method", ier)
   return_error = ier
   goto 1000
endif

! get number and names of species from KIM API object
!
ier = kim_api_get_num_sim_species(pkim, num_species, strlen)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_get_num_sim_species",ier)
   return_error = ier
   goto 1000
endif

! allocate storage for MEAM parameter tables
!
allocate(mp)
call allocate_meam_parameters(mp,num_species,350001)

do i=1,mp%max_comp
   ier = kim_api_get_sim_species(pkim, i, species)
   if (ier.lt.KIM_STATUS_OK) then
     idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
       "kim_api_get_sim_species", ier)
     return_error = ier
     goto 1000
   endif
   mp%Icomp(i) = trim(species)
!@
print *,i,trim(mp%Icomp(i))
enddo                                                ! which signals end of string

! set codes for species read in
do i=1,mp%max_comp
   call kim_api_set_species_code(pkim,mp%Icomp(i),i,ier)
   if (ier.lt.KIM_STATUS_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "kim_api_set_species_code",ier)
      return_error = ier
      goto 1000
   endif
enddo
                                                     ! in API
! Read in model parameters from parameter file
!
call Define_EAM_Potential(mp,paramfile(1),nmstrlen,mp%max_comp,mp%Icomp,ier)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "Unable to read MEAM parameters, kimerror = ",KIM_STATUS_FAIL)
   return_error = ier
   goto 1000
endif

! store pointer to parameters in KIM API object
!
call kim_api_set_model_buffer(pkim,c_loc(mp),ier)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_set_model_buffer",ier)
   return_error = ier
   goto 1000
endif

! store model cutoff in KIM object
pcutoff =  kim_api_get_data(pkim,"cutoff",ier)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_get_data", ier)
   return_error = ier
   goto 1000
endif
call c_f_pointer(pcutoff, cutoff)
cutoff = mp%Rcutoff

1000 continue
model_driver_init = return_error
return

end function model_driver_init

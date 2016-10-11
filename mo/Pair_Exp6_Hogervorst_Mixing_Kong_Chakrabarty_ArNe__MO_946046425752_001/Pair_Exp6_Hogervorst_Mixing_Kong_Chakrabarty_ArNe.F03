!****************************************************************************
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
! Copyright (c) 2012, Regents of the University of Minnesota.  All rights reserved.
!
!****************************************************************************
!
! MODULE Pair_Exp6_Hogervorst
!
!
! Compute energy and forces on an isolated cluster of atoms using 
! the exp-6 potential shifted to have zero energy at cutoff. 
! Parameters due to Hogervorst [1] and mixing rule between species 
! due to Kong and Chakrabarty [2].
!
! Author : E. B. Tadmor (10-MAR-12), updated to kim api v1.6 4-AUG-14.
!
! Support atom species:
!
!  1 = "Ar"
!  2 = "Ne"
!
!  References:
!
!  [1] W. Hogervorst, Physica, vol 51, 77-89 (1971). (Like interactions.)
!
!  [2] C. L. Kong and M. R. Chakrabarty, J. Phys. Chem., Vol. 77, 2668-2670
!      (1973). (Mixing rule.)
!
!****************************************************************************

#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__

module Pair_Exp6_Hogervorst

use, intrinsic :: iso_c_binding
use KIM_API_F03
implicit none

save
private
public Compute_Energy_Forces,cutoff

integer(c_int), parameter :: cd = c_double ! used for literal constants

!
! Define potential parameters
!
real(c_double), parameter :: kb     = 8.6173e-5_cd ! Boltzmann's constant [eV/K]
real(c_double), parameter :: cutoff = 8.15_cd      ! cutoff radius [A] (arbitrary value)
real(c_double), parameter :: eps11  = 138.0_cd*kb  ! Ar-Ar epsilon parameter [eV]
real(c_double), parameter :: rm11   = 3.77_cd      ! Ar-Ar rm parameter [A]
real(c_double), parameter :: alf11  = 14.8_cd      ! Ar-Ar alpha parameter [-]
real(c_double), parameter :: eps22  = 43.0_cd*kb   ! Ne-Ne epsilon parameter [eV]
real(c_double), parameter :: rm22   = 3.03_cd      ! Ne-Ne rm parameter [A]
real(c_double), parameter :: alf22  = 16.0_cd      ! Ne-Ne alpha parameter [-]
real(c_double), parameter :: eps12  = 68.89_cd*kb  ! Ar-Ne epsilon parameter [eV]
real(c_double), parameter :: rm12   = 3.447_cd     ! Ar-Ne rm parameter [A]
real(c_double), parameter :: alf12  = 15.52_cd     ! Ar-Ne alpha parameter [-]

contains

!-------------------------------------------------------------------------------
!
!  Calculate exp-6 pair potential phi(r) and its derivative dphi(r)
!  Note: potential is shifted to have zero energy at r=cutoff
!
!-------------------------------------------------------------------------------
subroutine calc_phi_dphi(r,phi,dphi,eps,rm,alf,cutoff)
implicit none

!-- Transferred variables
real(c_double), intent(in)  :: r
real(c_double), intent(in)  :: eps,rm,alf,cutoff
real(c_double), intent(out) :: phi,dphi

!-- Local variables
real(c_double) soa,amp,ror,cor,phicut

if (r .gt. cutoff) then
   ! Argument exceeds cutoff radius
   phi    = 0.0_cd
   dphi   = 0.0_cd
else
   soa = 6.0_cd/alf
   amp = eps/(1.0_cd-soa)
   cor = cutoff/rm
   ror = r/rm
   phicut  = amp*(soa*exp(alf*(1.0_cd-cor))-cor**(-6))
   phi     = amp*(soa*exp(alf*(1.0_cd-ror))-ror**(-6)) - phicut
   dphi    = -6.0_cd*amp/rm*(exp(alf*(1.0_cd-ror))-ror**(-7))
endif

return

end subroutine calc_phi_dphi

!-------------------------------------------------------------------------------
!
! Compute energy and forces on atoms from the positions.
!
!-------------------------------------------------------------------------------
integer(c_int) function Compute_Energy_Forces(pkim) bind(c)
implicit none

!--Transferred Variables
type(c_ptr), intent(in)  :: pkim

!--Local Variables
integer(c_int), parameter :: DIM=3
integer(c_int) ier
integer(c_int) i,j
real(c_double) Rij(DIM),Rsqij,cutoffsq,r,phi,dphi

!-- KIM variables
integer(c_int), pointer :: numberOfParticles;  type(c_ptr) :: pnAtoms
real(c_double), pointer :: energy;             type(c_ptr) :: penergy
real(c_double), pointer :: coor(:,:);          type(c_ptr) :: pcoor
real(c_double), pointer :: force(:,:);         type(c_ptr) :: pforce
integer(c_int), pointer :: particleSpecies(:); type(c_ptr) :: pparticleSpecies
integer(c_int) comp_energy,comp_force
integer(c_int) idum

!
! Check to see if we have been asked to compute the energy and forces
!
call kim_api_getm_compute(pkim, ier,   &
     "energy",         comp_energy, 1, &
     "forces",         comp_force,  1)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_getm_compute",ier)
   goto 42
endif

!
! Unpack data from KIM object
!
call kim_api_getm_data(pkim, ier, &
     "numberOfParticles", pnAtoms,           1,           &
     "particleSpecies",   pparticleSpecies,  1,           &
     "coordinates",       pcoor,             1,           &
     "energy",            penergy,           comp_energy, &
     "forces",            pforce,            comp_force)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_getm_data_f", ier)
   goto 42
endif

call c_f_pointer(pnAtoms,          numberOfParticles)
call c_f_pointer(pparticleSpecies, particleSpecies, [numberOfParticles])
call c_f_pointer(pcoor,            coor,            [DIM,numberOfParticles])
if (comp_energy.eq.1) call c_f_pointer(penergy, energy)
if (comp_force.eq.1)  call c_f_pointer(pforce,  force,  &
                                                    [DIM,numberOfParticles])
!
! Verify input
!
ier = KIM_STATUS_OK ! assume all is OK
do i=1,numberOfParticles
   if (particleSpecies(i).ne.1 .and. particleSpecies(i).ne.2) ier = KIM_STATUS_FAIL
enddo
if (ier.ne.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "Unexpected species detected", ier)
   goto 42
endif

!
! Compute energy and forces (use half-list approach)
!
if (comp_energy.eq.1) energy = 0.0_cd
if (comp_force.eq.1)  force  = 0.0_cd
cutoffsq = cutoff**2
do i=1,numberOfParticles-1
   do j=i+1,numberOfParticles
      Rij(:) = coor(:,j) - coor(:,i)
      Rsqij = dot_product(Rij,Rij)
      if ( Rsqij.lt.cutoffsq ) then
         r = sqrt(Rsqij)
         if (particleSpecies(i).eq.1 .and. particleSpecies(j).eq.1) then
            call calc_phi_dphi(r,phi,dphi,eps11,rm11,alf11,cutoff)
         elseif (particleSpecies(i).eq.2 .and. particleSpecies(j).eq.2) then
            call calc_phi_dphi(r,phi,dphi,eps22,rm22,alf22,cutoff)
         else
            call calc_phi_dphi(r,phi,dphi,eps12,rm12,alf12,cutoff)
         endif
         if (comp_energy.eq.1) energy = energy + phi
         if (comp_force.eq.1) then
            force(:,i) = force(:,i) + dphi*Rij/r
            force(:,j) = force(:,j) - dphi*Rij/r
         endif
      endif
   enddo
enddo

ier = KIM_STATUS_OK
42 continue
Compute_Energy_Forces = ier
return

end function Compute_Energy_Forces

end module Pair_Exp6_Hogervorst

!-------------------------------------------------------------------------------
!
! Model initialization routine (REQUIRED)
!
!-------------------------------------------------------------------------------
integer(c_int) function model_init(pkim) bind(c)
use, intrinsic :: iso_c_binding
use Pair_Exp6_Hogervorst
use KIM_API_F03
implicit none

!-- Transferred variables
type(c_ptr), intent(in) :: pkim

!-- Local variables
integer(c_int), parameter :: one=1
integer(c_int) ier, idum

!-- KIM variables
real(c_double), pointer :: rcut;  type(c_ptr) :: prcut

! store pointer to compute function in KIM object
ier = kim_api_set_method(pkim,"compute",one,c_funloc(Compute_Energy_Forces))
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_set_method", ier)
   goto 42
endif

! store model cutoff in KIM object
prcut =  kim_api_get_data(pkim,"cutoff",ier)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_get_data", ier)
   goto 42
endif
call c_f_pointer(prcut, rcut)
rcut = cutoff

ier = KIM_STATUS_OK
42 continue
model_init = ier
return

end function model_init

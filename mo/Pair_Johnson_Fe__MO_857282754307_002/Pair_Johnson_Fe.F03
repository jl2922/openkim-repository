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
! Copyright (c) 2013--2014, Regents of the University of Minnesota.
! All rights reserved.
!
! Contributors:
!    Ryan S. Elliott
!    Ellad B. Tadmor
!    Sung-Yup Kim
!

!****************************************************************************
!**
!**  MODULE iron_pair_potential
!**
!**  Reference: J. R. Morris, R. S. Aga, V. Levashov and T. Egami, 
!**             "Many-body effects in bcc metals: An embedded atom model 
!**             extension of the modified Johnson pair potential for iron"
!**             Phys. Rev. B 77, 174201 (2008).
!**
!**  Language: Fortran 2003
!**
!**  Release: This file is part of the kim-api.git repository.
!**
!****************************************************************************


#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))

module iron_pair_potential

use, intrinsic :: iso_c_binding
use KIM_API_F03
implicit none

save
private
public Compute_Energy_Forces, &
       model_cutoff

! Below are the definitions and values of all Model parameters
integer(c_int), parameter :: cd = c_double ! used for literal constants
integer(c_int), parameter :: DIM=3          ! dimensionality of space
integer(c_int), parameter :: speccode = 1   ! internal species code
real(c_double), parameter :: model_cutoff  =  3.44_cd
real(c_double), parameter :: model_cutsq   = model_cutoff**2

!-------------------------------------------------------------------------------
! Below are the definitions and values of all additional model parameters
!
! Recall that the Fortran 2003 format for declaring parameters is as follows:
!
! integer(c_int), parameter :: parname = value   ! This defines an integer
!                                                ! parameter called `parname'
!                                                ! with a value equal to
!                                                ! `value' (a number)
!
! real(c_double), parameter :: parname = value   ! This defines a real(c_double)
!                                                ! parameter called `parname'
!                                                ! with a value equal to
!                                                ! `value' (a number)
!-------------------------------------------------------------------------------
real(c_double), parameter :: a0  = -0.2002108_cd
real(c_double), parameter :: a1  = -0.5047747_cd
real(c_double), parameter :: a2  =    1.372738_cd
real(c_double), parameter :: a3  =  -15.09618_cd
real(c_double), parameter :: a4  =  -12.90021_cd
real(c_double), parameter :: b0  =   -1.581570_cd
real(c_double), parameter :: b1  =    0.477871_cd
real(c_double), parameter :: b2  =   -0.639230_cd
real(c_double), parameter :: c0  =   -0.1469636_cd
real(c_double), parameter :: c1  =    0.4521426_cd
real(c_double), parameter :: c2  =    0.2221241_cd
real(c_double), parameter :: c3  =    1.725326_cd
real(c_double), parameter :: c4  = - 12.91063_cd
real(c_double), parameter :: c5  =   14.67111_cd
real(c_double), parameter :: aa0 = 8752.934_cd
real(c_double), parameter :: B   =    4.572488_cd
contains

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi(r,phi)
implicit none

!-- Transferred variables
real(c_double), intent(in)  :: r
real(c_double), intent(out) :: phi

!-- Local variables
! <FILL place any local variable definitions here>

if (r <  1.90_cd) then

   phi = aa0*exp(-B*r)

else if ( r < 2.40_cd) then

   phi = a0 +  a1*(r-2.4_cd) + a2*(r-2.4_cd)**2 + a3*(r-2.4_cd)**3 + a4*(r-2.4_cd)**4

else if (  r < 3.00_cd) then

   phi = b2*(r-3.115829_cd)**3 +  b1*r +  b0

else if ( r <  3.44_cd) then

   phi = c0 +  c1*(r-3.0_cd) + c2*(r-3.0_cd)**2 + c3*(r-3.0_cd)**3 + &
         c4*(r-3.0_cd)**4 + c5*(r-3.0_cd)**5

else

  phi = 0.0_cd

endif

end subroutine calc_phi

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r) and its derivative dphi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi_dphi(r,phi,dphi)
implicit none

!-- Transferred variables
real(c_double), intent(in)  :: r
real(c_double), intent(out) :: phi,dphi

!-- Local variables
! <FILL place any local variable definitions here>

if (r <  1.90_cd) then

   phi = aa0*exp(-B*r)
   dphi = -B*aa0*exp(-B*r)

else if ( r < 2.40_cd) then

   Phi  = a0 +  a1*(r-2.4) + a2*(r-2.4_cd)**2 + a3*(r-2.4_cd)**3 + a4*(r-2.4_cd)**4
   dphi = a1 + 2.0_cd*a2*(r-2.4) + 3.0_cd*a3*(r-2.4_cd)**2 + 4.0_cd*a4*(r-2.4_cd)**3

else if (  r < 3.00_cd) then

   phi  = b2*(r-3.115829_cd)**3 +  b1*r +  b0
   dphi = 3.0_cd*b2*(r-3.115829_cd)**2 + b1

else if ( r <  3.44_cd) then

   phi  = c0 + c1*(r-3.0) + c2*(r-3.0_cd)**2 + c3*(r-3.0_cd)**3 + &
          c4*(r-3.0_cd)**4 + c5*(r-3.0_cd)**5
   dphi = c1 + 2.0_cd*c2*(r-3.0_cd) + 3.0_cd*c3*(r-3.0_cd)**2 +   &
          4.0_cd*c4*(r-3.0_cd)**3 + 5.0_cd*c5*(r-3.0_cd)**4

else

  phi  = 0.0_cd
  dphi = 0.0_cd

endif


end subroutine calc_phi_dphi

!-------------------------------------------------------------------------------
!
! Compute energy and forces on atoms from the positions.
!
!-------------------------------------------------------------------------------
integer(c_int) function Compute_Energy_Forces(pkim) bind(c)
implicit none

!-- Transferred variables
type(c_ptr), intent(in)  :: pkim

!-- Local variables
real(c_double) :: Rij(DIM)
real(c_double) :: r,Rsqij,phi,dphi,dEidr = 0.0_cd
integer(c_int) :: i,j,jj,numnei,atom_ret,comp_force,comp_enepot,comp_virial, &
                  comp_energy
integer(c_int) :: ier
character(len=80) :: error_message

!-- KIM variables
integer(c_int), pointer :: N;                 type(c_ptr) :: pN
real(c_double), pointer :: energy;            type(c_ptr) :: penergy
real(c_double), pointer :: coor(:,:);         type(c_ptr) :: pcoor
real(c_double), pointer :: force(:,:);        type(c_ptr) :: pforce
real(c_double), pointer :: enepot(:);         type(c_ptr) :: penepot
real(c_double), pointer :: boxSideLengths(:); type(c_ptr) :: pboxSideLengths
real(c_double), pointer :: Rij_list(:,:);     type(c_ptr) :: pRij_list
integer(c_int), pointer :: numContrib;        type(c_ptr) :: pnumContrib
integer(c_int), pointer :: nei1atom(:);       type(c_ptr) :: pnei1atom
integer(c_int), pointer :: particleSpecies(:);type(c_ptr) :: pparticleSpecies
real(c_double), pointer :: virial(:);         type(c_ptr) :: pvirial
character(len=KIM_KEY_STRING_LENGTH) NBC_Method
integer(c_int) IterOrLoca
integer(c_int) HalfOrFull
integer(c_int) NBC
integer(c_int) numberContrib
integer(c_int) idum

numberContrib = 0 ! initialize

! Determine neighbor list boundary condition (NBC)
! and half versus full mode:
! *****************************
! * HalfOrFull = 1 -- Half
! *            = 2 -- Full
! *****************************
!
!
ier = kim_api_get_nbc_method(pkim, NBC_Method)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_get_nbc_method", ier)
   goto 42
endif
if (index(NBC_Method,"CLUSTER").eq.1) then
   NBC = 0
   HalfOrFull = 1
elseif (index(NBC_Method,"MI_OPBC_H").eq.1) then
   NBC = 1
   HalfOrFull = 1
elseif (index(NBC_Method,"MI_OPBC_F").eq.1) then
   NBC = 1
   HalfOrFull = 2
elseif (index(NBC_Method,"NEIGH_PURE_H").eq.1) then
   NBC = 2
   HalfOrFull = 1
elseif (index(NBC_Method,"NEIGH_PURE_F").eq.1) then
   NBC = 2
   HalfOrFull = 2
elseif (index(NBC_Method,"NEIGH_RVEC_F").eq.1) then
   NBC = 3
   HalfOrFull = 2
else
   ier = KIM_STATUS_FAIL
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "Unknown NBC method", ier)
   goto 42
endif

! Determine neighbor list handling mode
!
if (NBC.ne.0) then
   !*****************************
   !* IterOrLoca = 1 -- Iterator
   !*            = 2 -- Locator
   !*****************************
   IterOrLoca = kim_api_get_neigh_mode(pkim, ier)
   if (ier.lt.KIM_STATUS_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "kim_api_get_neigh_mode", ier)
      goto 42
   endif
   if (IterOrLoca.ne.1 .and. IterOrLoca.ne.2) then
      ier = KIM_STATUS_FAIL
      write(error_message,'(a,i1)') &
         'Unsupported IterOrLoca mode = ',IterOrLoca
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  error_message, ier)
      goto 42
   endif
else
   IterOrLoca = 2   ! for CLUSTER NBC
endif

! Check to see if we have been asked to compute the forces, energyperatom,
! energy and virial
!
call kim_api_getm_compute(pkim, ier, &
     "energy",         comp_energy, 1, &
     "forces",         comp_force,  1, &
     "particleEnergy", comp_enepot, 1, &
     "virial",         comp_virial, 1)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_getm_compute", ier)
   goto 42
endif

! Unpack data from KIM object
!
call kim_api_getm_data(pkim, ier, &
     "numberOfParticles",           pN,              1,                           &
     "particleSpecies",             pparticleSpecies,1,                           &
     "coordinates",                 pcoor,           1,                           &
     "numberContributingParticles", pnumContrib,     TRUEFALSE(HalfOrFull.eq.1),  &
     "boxSideLengths",              pboxSideLengths, TRUEFALSE(NBC.eq.1),         &
     "energy",                      penergy,         TRUEFALSE(comp_energy.eq.1), &
     "forces",                      pforce,          TRUEFALSE(comp_force.eq.1),  &
     "particleEnergy",              penepot,         TRUEFALSE(comp_enepot.eq.1), &
     "virial",                      pvirial,         TRUEFALSE(comp_virial.eq.1))
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_getm_data_f", ier)
   goto 42
endif

call c_f_pointer(pN,               N)
call c_f_pointer(pparticleSpecies, particleSpecies, [N])
call c_f_pointer(pcoor,            coor,            [DIM,N])
if (HalfOrFull.eq.1)  call c_f_pointer(pnumContrib,     numContrib)
if (NBC.eq.1)         call c_f_pointer(pboxSideLengths, boxSideLengths, [DIM])
if (comp_energy.eq.1) call c_f_pointer(penergy,         energy)
if (comp_force.eq.1)  call c_f_pointer(pforce,          force,          [DIM,N])
if (comp_enepot.eq.1) call c_f_pointer(penepot,         enepot,         [N])
if (comp_virial.eq.1) call c_f_pointer(pvirial,         virial,         [6])

if (HalfOrFull.eq.1) then
   if (NBC.ne.0) then ! non-CLUSTER cases
      numberContrib = numContrib
   else               ! CLUSTER case
      numberContrib = N
   endif
endif

! Check to be sure that the atom types are correct
!
ier = KIM_STATUS_FAIL ! assume an error
do i = 1,N
   if (particleSpecies(i).ne.speccode) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "Unexpected species detected", ier)
      goto 42
   endif
enddo
ier = KIM_STATUS_OK ! everything is ok

! Initialize potential energies, forces, virial term
!
if (comp_enepot.eq.1) enepot = 0.0_cd
if (comp_energy.eq.1) energy = 0.0_cd
if (comp_force.eq.1)  force  = 0.0_cd
if (comp_virial.eq.1) virial = 0.0_cd

! Initialize neighbor handling for CLUSTER NBC
!
if (NBC.eq.0) then
   allocate( nei1atom(N) )
endif

! Initialize neighbor handling for Iterator mode
!
if (IterOrLoca.eq.1) then
   ier = kim_api_get_neigh(pkim,0,0,atom_ret,numnei, &
                           pnei1atom,pRij_list)
   ! check for successful initialization
   if (ier.ne.KIM_STATUS_NEIGH_ITER_INIT_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "kim_api_get_neigh", ier)
      ier = KIM_STATUS_FAIL
      goto 42
   endif
endif

!
!  Compute energy and forces
!

!  Loop over particles and compute energy and forces
!
i = 0
do

   ! Set up neighbor list for next atom for all NBC methods
   !
   if (IterOrLoca.eq.1) then    ! ITERATOR mode
      ier = kim_api_get_neigh(pkim,0,1,atom_ret,numnei, &
                              pnei1atom,pRij_list)
      if (ier.eq.KIM_STATUS_NEIGH_ITER_PAST_END) exit
                                ! incremented past the end of the list,
                                ! terminate loop
      if (ier.lt.KIM_STATUS_OK) then     ! some sort of problem, exit
         idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                     "kim_api_get_neigh", ier)
         goto 42
      endif

      i = atom_ret

   else                         ! LOCATOR mode
      i = i + 1
      if (i.gt.N) exit          ! incremented past end of list,
                                ! terminate loop
      if (NBC.eq.0) then     ! CLUSTER NBC method
         numnei = N - i      ! number of neighbors in list i+1, ..., N
         nei1atom(1:numnei) = (/ (i+jj, jj = 1,numnei) /)
         ier = KIM_STATUS_OK
      else
         ier = kim_api_get_neigh(pkim,1,i,atom_ret,numnei, &
                                 pnei1atom,pRij_list)
         if (ier.ne.KIM_STATUS_OK) then ! some sort of problem, exit
            idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                        "kim_api_get_neigh", ier)
            ier = KIM_STATUS_FAIL
            goto 42
         endif
      endif
   endif

   if (NBC.ne.0) call c_f_pointer(pnei1atom, nei1atom, [numnei])
   if (NBC.eq.3) call c_f_pointer(pRij_list, Rij_list, [DIM,numnei])

   ! Loop over the neighbors of atom i
   !
   do jj = 1, numnei

      j = nei1atom(jj)                            ! get neighbor ID

      ! compute relative position vector
      !
      if (NBC.ne.3) then                          ! all methods except NEIGH_RVEC_F
         Rij(:) = coor(:,j) - coor(:,i)           ! distance vector between i j
      else
         Rij(:) = Rij_list(:,jj)
      endif

      ! apply periodic boundary conditions if required
      !
      if (NBC.eq.1) then
         where ( abs(Rij) .gt. 0.5_cd*boxSideLengths )
                                                       ! periodic boundary conditions
            Rij = Rij - sign(boxSideLengths,Rij)       ! applied where needed.
         end where
      endif

      ! compute energy and forces
      !
      Rsqij = dot_product(Rij,Rij)                ! compute square distance
      if ( Rsqij .lt. model_cutsq ) then          ! particles are interacting?

         r = sqrt(Rsqij)                          ! compute distance
         if (comp_force.eq.1.or.comp_virial.eq.1) then
            call calc_phi_dphi(r,phi,dphi)        ! compute pair potential
                                                  !   and it derivative
            if ((HalfOrFull.eq.1) .and. &
                (j .le. numberContrib)) then      ! HALF mode
               dEidr = dphi                       !      double contribution
            else                                  ! FULL mode
               dEidr = 0.5_cd*dphi                 !      regular contribution
            endif
         else
            call calc_phi(r,phi)                  ! compute just pair potential
         endif

         ! contribution to energy
         !
         if (comp_enepot.eq.1) then
            enepot(i) = enepot(i) + 0.5_cd*phi   ! accumulate energy
            if ((HalfOrFull.eq.1) .and. &
                (j .le. numberContrib)) &        ! HALF mode
               enepot(j) = enepot(j) + 0.5_cd*phi! (i and j share it)
         endif
         if (comp_energy.eq.1) then
            if ((HalfOrFull.eq.1) .and. &
                (j .le. numberContrib)) then     ! HALF mode
               energy = energy + phi             !      add v to total energy
            else                                 ! FULL mode
               energy = energy + 0.5_cd*phi       !      add half v to total energy
            endif
         endif

         ! contribution to virial tensor, virial(i,j)=r(i)*r(j)*(dV/dr)/r
         !
         if (comp_virial.eq.1) then
            virial(1) = virial(1) + Rij(1)*Rij(1)*dEidr/r
            virial(2) = virial(2) + Rij(2)*Rij(2)*dEidr/r
            virial(3) = virial(3) + Rij(3)*Rij(3)*dEidr/r
            virial(4) = virial(4) + Rij(2)*Rij(3)*dEidr/r
            virial(5) = virial(5) + Rij(1)*Rij(3)*dEidr/r
            virial(6) = virial(6) + Rij(1)*Rij(2)*dEidr/r
         endif

         ! contribution to forces
         !
         if (comp_force.eq.1) then
            force(:,i) = force(:,i) + dEidr*Rij/r ! accumulate force on atom i
            force(:,j) = force(:,j) - dEidr*Rij/r ! accumulate force on atom j
         endif

      endif

   enddo  ! loop on jj

enddo  ! infinite do loop (terminated by exit statements above)

! Free temporary storage
!
if (NBC.eq.0) deallocate( nei1atom )

! Everything is great
!
ier = KIM_STATUS_OK
42 continue
Compute_Energy_Forces = ier
return

end function Compute_Energy_Forces

end module iron_pair_potential

!-------------------------------------------------------------------------------
!
! Model initialization routine (REQUIRED)
!
!-------------------------------------------------------------------------------
integer(c_int) function model_init(pkim) bind(c)
use, intrinsic :: iso_c_binding
use iron_pair_potential
use KIM_API_F03
implicit none

!-- Transferred variables
type(c_ptr), intent(in) :: pkim

!-- Local variables
integer(c_int), parameter :: one=1
integer(c_int) ier, idum

!-- KIM variables
real(c_double), pointer :: cutoff; type(c_ptr) :: pcutoff

! store pointer to compute function in KIM object
ier = kim_api_set_method(pkim,"compute",one,c_funloc(Compute_Energy_Forces))
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_set_method", ier)
   goto 42
endif

! store model cutoff in KIM object
pcutoff =  kim_api_get_data(pkim,"cutoff",ier)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_get_data", ier)
   goto 42
endif
call c_f_pointer(pcutoff, cutoff)
cutoff = model_cutoff

ier = KIM_STATUS_OK
42 continue
model_init = ier
return

end function model_init

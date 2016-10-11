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
!
! CDDL HEADER END
!

!
! Copyright (c) 2012, Mark R. Gilbert, CCFE Fusion Association.
! All rights reserved.
!
! Contributors:
!    Mark R. Gilbert
!    Ryan S. Elliott
!

!****************************************************************************
!**
!**  MODULE model_driver_PF_cubic_splines
!**
!**  EAM-like potential with cubic splines representing knot functions
!**  magnetic ability also available via B parameter
!**
!**  Language: Fortran 2003
!**
!**
!**
!****************************************************************************


#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))

module model_driver_pf_cubic_splines

use, intrinsic :: iso_c_binding
use KIM_API_F03
implicit none

save
private
public Compute_Energy_Forces, &
       Destroy, &
       model_parameters

type model_parameters
  real(c_double) :: model_cutoff  ! cutoff radius in angstroms
  real(c_double) :: model_cutsq
  ! new variables for general potential - to be read in - 12/10/06
  real(c_double)       :: A, B, Z,r1,r2,a_inter(7),a_rho
  integer(c_int)       :: nknotv,nknotp,interpolate_num
  logical              :: linear_interpolate
  real(c_double) , pointer :: vknotcoeff(:), vknotpoint(:), &
    pknotcoeff(:), pknotpoint(:)
end type model_parameters

! Below are the definitions and values of all Model parameters
integer(c_int), parameter          :: cd = c_double  ! used for literal constants
integer(c_int), parameter          :: DIM=3          ! dimensionality of space
integer(c_int), parameter          :: speccode = 1   ! internal species code

    real(c_double) , parameter :: pc=1.0_cd, pi=3.141592653589793_cd
! electric constant and electron charge ( units of
    real(c_double) , parameter :: E_0=8.854187817e-12_cd, Ec=1.60217653e-19_cd
!bohr constant in angstroms
   real(c_double) , parameter :: a_B=0.5291772108_cd

! constants of biersack-ziegler coulomb potential
    real(c_double) , parameter :: p1=0.1818_cd, p2=0.5099_cd, p3=0.2802_cd, p4=0.02817_cd
    real(c_double) , parameter :: Bt1=-3.2_cd, Bt2=-0.9423_cd, Bt3=-0.4029_cd, Bt4=-0.2016_cd
!

contains

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi(r,params,phi)
implicit none

!-- Transferred variables
real(c_double),         intent(in)  :: r
type(model_parameters), intent(in)  :: params
real(c_double),         intent(out) :: phi

!-- Local variables
real(c_double) :: a_s,phi_r,B_Zpre
integer(c_int) :: i

if (r .gt. params%model_cutoff) then
   ! Argument exceeds cutoff radius
   phi = 0.0_cd
else
    IF (r>=params%r2) THEN
      phi=0.0_cd
      DO i=1,params%nknotv
        IF (r<params%vknotpoint(i)) THEN
         phi=phi+params%vknotcoeff(i)*(params%vknotpoint(i)-r)**3
        END IF
      END DO
    ELSE IF (r<=params%r1) THEN
     a_s=0.88534_cd*a_B/(sqrt(2.0_cd)*(params%Z**(1.0_cd/3.0_cd)))
     B_Zpre=params%Z**2*Ec*1.0e10_cd/(4.0_cd*pi*E_0)
     phi_r=p1*exp(Bt1*r/a_s)+p2*exp(Bt2*r/a_s)+p3*exp(Bt3*r/a_s)+ &
               p4*exp(Bt4*r/a_s)
     phi=B_Zpre*phi_r/r
    ELSE
     IF (params%linear_interpolate) THEN
      phi=params%a_inter(1)+params%a_inter(2)*r+params%a_inter(3)*(r**2)+ &
       params%a_inter(4)*(r**3)+params%a_inter(5)*(r**4)+params%a_inter(6)*(r**5)+params%a_inter(7)*(r**6)
     ELSE
      phi=exp(params%a_inter(1)+params%a_inter(2)*r+params%a_inter(3)*(r**2)+ &
       params%a_inter(4)*(r**3)+params%a_inter(5)*(r**4)+params%a_inter(6)*(r**5)+params%a_inter(7)*(r**6))
     END IF
    END IF
endif

end subroutine calc_phi

!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(r) and its derivative dphi(r)
!
!-------------------------------------------------------------------------------
subroutine calc_phi_dphi(r,params,phi,dphi)
implicit none

!-- Transferred variables
real(c_double),         intent(in)  :: r
type(model_parameters), intent(in)  :: params
real(c_double),         intent(out) :: phi,dphi

!-- Local variables
real(c_double) :: a_s,phi_r,B_Zpre,dphi_r
integer(c_int) :: i

if (r .gt. params%model_cutoff) then
   ! Argument exceeds cutoff radius
   phi    = 0.0_cd
   dphi   = 0.0_cd
else
    IF (r>=params%r2) THEN
      dphi = 0.0_cd
      phi  = 0.0_cd
      DO i=1,params%nknotv
        IF (r<params%vknotpoint(i)) THEN
         dphi=dphi-3*params%vknotcoeff(i)*(params%vknotpoint(i)-r)**2
         phi=phi+params%vknotcoeff(i)*(params%vknotpoint(i)-r)**3
        END IF
      END DO
    ELSE IF (r<=params%r1) THEN
     a_s=0.88534_cd*a_B/(sqrt(2.0_cd)*(params%Z**(1.0_cd/3.0_cd)))
     B_Zpre=params%Z**2*Ec*1.0e10_cd/(4.0_cd*pi*E_0)
     phi_r=p1*exp(Bt1*r/a_s)+p2*exp(Bt2*r/a_s)+p3*exp(Bt3*r/a_s)+ &
               p4*exp(Bt4*r/a_s)
     dphi_r=(p1*Bt1*exp(Bt1*r/a_s)+p2*Bt2*exp(Bt2*r/a_s)+p3*Bt3*exp(Bt3*r/a_s)+ &
               p4*Bt4*exp(Bt4*r/a_s))/a_s
     dphi=-B_Zpre*phi_r/(r**2)+B_Zpre*dphi_r/r
     phi=B_Zpre*phi_r/r
    ELSE
     IF (params%linear_interpolate) THEN
      dphi=params%a_inter(2)+2.0_cd*params%a_inter(3)*r+3.0_cd*params%a_inter(4)*(r**2)+ &
                      4.0_cd*params%a_inter(5)*(r**3)+5.0_cd*params%a_inter(6)*(r**4)+6.0_cd*params%a_inter(7)*(r**5)
      phi=params%a_inter(1)+params%a_inter(2)*r+params%a_inter(3)*(r**2)+ &
          params%a_inter(4)*(r**3)+params%a_inter(5)*(r**4)+params%a_inter(6)*(r**5)+params%a_inter(7)*(r**6)
     ELSE
       phi_r=exp(params%a_inter(1)+params%a_inter(2)*r+params%a_inter(3)*(r**2)+ &
        params%a_inter(4)*(r**3)+params%a_inter(5)*(r**4)+params%a_inter(6)*(r**5)+params%a_inter(7)*(r**6))
       dphi_r=(params%a_inter(2)+2.0_cd*params%a_inter(3)*r+3.0_cd*params%a_inter(4)*(r**2)+ &
                      4.0_cd*params%a_inter(5)*(r**3)+5.0_cd*params%a_inter(6)*(r**4)+6.0_cd*params%a_inter(7)*(r**5))
       dphi=phi_r*dphi_r
       phi=exp(params%a_inter(1)+params%a_inter(2)*r+params%a_inter(3)*(r**2)+ &
          params%a_inter(4)*(r**3)+params%a_inter(5)*(r**4)+params%a_inter(6)*(r**5)+params%a_inter(7)*(r**6))
     END IF

    END IF
endif

end subroutine calc_phi_dphi

!-------------------------------------------------------------------------------
!
!  Calculate electron density g(r)
!
!-------------------------------------------------------------------------------
subroutine calc_g(r,params,g)
implicit none

!-- Transferred variables
real(c_double),         intent(in)  :: r
type(model_parameters), intent(in)  :: params
real(c_double),         intent(out) :: g

!-- Local variables
integer(c_int) :: i

if (r .gt. params%model_cutoff) then
   ! Argument exceeds cutoff radius
   g = 0.0_cd
else
   g = 0.0_cd
    DO i=1,params%nknotp
      IF (r<params%pknotpoint(i)) THEN
       g=g+params%pknotcoeff(i)*(params%pknotpoint(i)-r)**3
      END IF
    END DO
endif

end subroutine calc_g

!-------------------------------------------------------------------------------
!
!  Calculate electron density derivative dg(r)
!
!-------------------------------------------------------------------------------
subroutine calc_dg(r,params,dg)
implicit none

!-- Transferred variables
real(c_double),         intent(in)  :: r
type(model_parameters), intent(in)  :: params
real(c_double),         intent(out) :: dg

!-- Local variables
integer(c_int) :: i

if (r .gt. params%model_cutoff) then
   ! Argument exceeds cutoff radius
   dg = 0.0_cd
else
   dg = 0.0_cd
    DO i=1,params%nknotp
      IF (r<params%pknotpoint(i)) THEN
       dg=dg-3.0_cd*params%pknotcoeff(i)*(params%pknotpoint(i)-r)**2
      END IF
    END DO
endif

end subroutine calc_dg

!-------------------------------------------------------------------------------
!
!  Calculate embedding function U(rho)
!
!-------------------------------------------------------------------------------
subroutine calc_U(rho,params,U)
implicit none

!-- Transferred variables
real(c_double),         intent(in)  :: rho
type(model_parameters), intent(in)  :: params
real(c_double),         intent(out) :: U

    IF (rho<pc) THEN
     U=-params%A*sqrt(rho)-(params%B/log(2.0_cd))*(1.0_cd-sqrt(rho/pc))*log(2.0_cd-rho/pc)+ &
                     params%a_rho*(rho**2)
    ELSE
      U=-params%A*sqrt(rho)+params%a_rho*(rho**2)
    END IF

end subroutine calc_U

!-------------------------------------------------------------------------------
!
!  Calculate embedding function U(rho) and first derivative dU(rho)
!
!-------------------------------------------------------------------------------
subroutine calc_U_dU(rho,params,U,dU)
implicit none

!-- Transferred variables
real(c_double),         intent(in)  :: rho
type(model_parameters), intent(in)  :: params
real(c_double),         intent(out) :: U,dU

 IF(rho.le.1e-10_cd) THEN
    dU = 0.0_cd
    U  = 0.0_cd
 ELSE
    IF (rho<pc) THEN
     U=-params%A*sqrt(rho)-(params%B/log(2.0_cd))*(1.0_cd-sqrt(rho/pc))*log(2.0_cd-rho/pc)+ &
                     params%a_rho*(rho**2)
     dU=-params%A/(2.0_cd*sqrt(rho))-(params%B/log(2.0_cd))*( &
               (-1.0_cd/(2.0_cd*sqrt(rho*pc)))*log(2.0_cd-rho/pc)+ &
               (1.0_cd-sqrt(rho/pc))*(-1.0_cd/(2.0_cd*pc-rho)))+2.0_cd*params%a_rho*rho
    ELSE
      U=-params%A*sqrt(rho)+params%a_rho*(rho**2)
      dU=-params%A/(2.0_cd*sqrt(rho))+2.0_cd*params%a_rho*rho
    END IF
 END IF




end subroutine calc_U_dU


!-------------------------------------------------------------------------------
!
!  Calculate embedding function U(rho) and first +2nd derivative dU(rho)
!
!-------------------------------------------------------------------------------
subroutine calc_U_dU_ddU(rho,params,U,dU,ddU)
implicit none

!-- Transferred variables
real(c_double),         intent(in)  :: rho
type(model_parameters), intent(in)  :: params
real(c_double),         intent(out) :: U,dU,ddU

!-- Local variables
! <FILL place any local variable definitions here>
 IF(rho.le.1e-10_cd) THEN
    dU  = 0.0_cd
    U   = 0.0_cd
    ddU = 0.0_cd
 ELSE
    IF (rho<pc) THEN
     U=-params%A*sqrt(rho)-(params%B/log(2.0_cd))*(1.0_cd-sqrt(rho/pc))*log(2.0_cd-rho/pc)+ &
                     params%a_rho*(rho**2)
     dU=-params%A/(2.0_cd*sqrt(rho))-(params%B/log(2.0_cd))*( &
               (-1.0_cd/(2.0_cd*sqrt(rho*pc)))*log(2.0_cd-rho/pc)+ &
               (1.0_cd-sqrt(rho/pc))*(-1.0_cd/(2.0_cd*pc-rho)))+2.0_cd*params%a_rho*rho
     ddU=params%A/(4.0_cd*(rho**(3.0_cd/2.0_cd)))-(params%B/log(2.0_cd))*( &
                (rho**(-3.0_cd/2.0_cd)/(4.0_cd*sqrt(pc)))*log(2.0_cd-rho/pc)+ &
                 (1.0_cd/sqrt(pc*rho))*(1.0_cd/(2.0_cd*pc-rho))- &
                  (1.0_cd-sqrt(rho/pc))*(1.0_cd/((2.0_cd*pc-rho)**2)))+2.0_cd*params%a_rho
    ELSE
      U=-params%A*sqrt(rho)+params%a_rho*(rho**2)
      dU=-params%A/(2.0_cd*sqrt(rho))+2.0_cd*params%a_rho*rho
      ddU=params%A/(4.0_cd*(rho**(3.0_cd/2.0_cd)))+2.0_cd*params%a_rho
    END IF
 END IF




end subroutine calc_U_dU_ddU




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
real(c_double) :: r,Rsqij,phi,dphi,g,dg,dU,U,dphieff
real(c_double) :: dphii,dUi,Ei,dphij,dUj,Ej
integer(c_int) :: ier
integer(c_int) :: i,j,jj,numnei,comp_force,comp_enepot,comp_virial,comp_energy
real(c_double), allocatable :: rho(:),derU(:)
integer(c_int), allocatable, target :: nei1atom_substitute(:)
character(len=80) :: error_message

!-- KIM variables
integer(c_int), pointer :: N;                  type(c_ptr) :: pN
real(c_double), pointer :: energy;             type(c_ptr) :: penergy
real(c_double), pointer :: coor(:,:);          type(c_ptr) :: pcoor
real(c_double), pointer :: force(:,:);         type(c_ptr) :: pforce
real(c_double), pointer :: ene_pot(:);         type(c_ptr) :: penepot
real(c_double), pointer :: boxSideLengths(:);  type(c_ptr) :: pboxSideLengths
real(c_double), pointer :: Rij_list(:,:);      type(c_ptr) :: pRij_list
integer(c_int), pointer :: numContrib;         type(c_ptr) :: pnumContrib
integer(c_int), pointer :: nei1atom(:);        type(c_ptr) :: pnei1atom
integer(c_int), pointer :: particleSpecies(:); type(c_ptr) :: pparticleSpecies
real(c_double), pointer :: virial_global(:);   type(c_ptr) :: pvirial
character(len=KIM_KEY_STRING_LENGTH) NBC_Method
integer(c_int) IterOrLoca
integer(c_int) HalfOrFull
integer(c_int) NBC
integer(c_int) numberContrib
integer(c_int) idum
integer(c_int) atom_ret

type(model_parameters), pointer :: params; type(c_ptr) :: pparams

pparams = kim_api_get_model_buffer(pkim, ier)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_get_model_buffer", ier)
   goto 42
endif
call c_f_pointer(pparams, params)

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
     "virial",                      pvirial,         TRUEFALSE(comp_virial.eq.1)  &
     )
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_getm_data", ier)
   goto 42
endif

call c_f_pointer(pN, N)
call c_f_pointer(pparticleSpecies, particleSpecies, [N])
call c_f_pointer(pcoor, coor, [DIM,N])
if (HalfOrFull.eq.1)  call c_f_pointer(pnumContrib, numContrib)
if (NBC.eq.1)         call c_f_pointer(pboxSideLengths, boxSideLengths, [DIM])
if (comp_energy.eq.1) call c_f_pointer(penergy, energy)
if (comp_force.eq.1)  call c_f_pointer(pforce,force, [DIM,N])
if (comp_enepot.eq.1) call c_f_pointer(penepot,ene_pot, [N])
if (comp_virial.eq.1) call c_f_pointer(pvirial,virial_global, [6])
if (HalfOrFull.eq.1) then
   if (NBC.ne.0) then ! non-CLUSTER cases
      numberContrib = numContrib
   else               ! CLUSTER cases
      numberContrib = N
   endif
endif

! Check to be sure that the atom types are correct
!
ier = KIM_STATUS_FAIL ! assume an error
do i = 1,N
   if (particleSpecies(i).ne.speccode) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "Unexpected species type detected", ier)
      goto 42
   endif
enddo
ier = KIM_STATUS_OK ! everything is ok

! Initialize potential energies, forces, virial term, electron density
!
! Note: that the variable `ene_pot' does not need to be initialized
!       because it's initial value is set during the embedding energy
!       calculation.
!
if (comp_energy.eq.1) energy           = 0.0_cd
if (comp_force.eq.1)  force(1:DIM,1:N) = 0.0_cd
if (comp_virial.eq.1) virial_global    = 0.0_cd
allocate( rho(N) )  ! pair functional electron density
rho(1:N) = 0.0_cd
if (comp_force.eq.1.or.comp_virial.eq.1) allocate( derU(N) )  ! EAM embedded energy deriv

! Initialize neighbor handling for CLUSTER NBC
!
if (NBC.eq.0) then
   allocate( nei1atom_substitute(N) )
   nei1atom => nei1atom_substitute
endif

!
!  Compute energy and forces
!

! Reset iterator if one is being used
!
if (IterOrLoca.eq.1) then
   ier = kim_api_get_neigh(pkim,0,0,atom_ret,numnei,pnei1atom,pRij_list)
   if (ier.lt.KIM_STATUS_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "kim_api_get_neigh", ier)
      goto 42
   endif
endif

!  Loop over particles in the neighbor list a first time,
!  to compute electron density (=coordination)
!
i = 0
do

   ! Set up neighbor list for next atom for all NBC methods
   !
   call get_current_atom_neighbors(IterOrLoca,NBC,N,pkim, &
                                   i,numnei,nei1atom,Rij_list,ier)
   if (ier.eq.KIM_STATUS_NEIGH_ITER_PAST_END) exit  ! atom counter incremented past end of list
   if (ier.lt.KIM_STATUS_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "get_current_atom_neighbors", ier)
      goto 42
   endif

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
         where ( abs(Rij) .gt. 0.5_cd*boxSideLengths )  ! periodic boundary conditions
            Rij = Rij - sign(boxSideLengths,Rij)       ! applied where needed.
         end where                                !
      endif

      ! compute contribution to electron density
      !
      Rsqij = dot_product(Rij,Rij)                ! compute square distance
      if ( Rsqij .lt. params%model_cutsq ) then   ! particles are interacting?
         r = sqrt(Rsqij)                          ! compute distance
         call calc_g(r,params,g)                  ! compute electron density
         rho(i) = rho(i) + g                      ! accumulate electron density
         if ((HalfOrFull.eq.1) .and. &
             (j .le. numberContrib)) &            ! HALF mode
            rho(j) = rho(j) + g                   !      (add contrib to j)
      endif

   enddo  ! loop on jj

enddo  ! infinite do loop (terminated by exit statements above)

!  Now that we know the electron densities, calculate embedding part of energy
!  U and its derivative U' (derU)
!
do i = 1,N
   if (comp_force.eq.1.or.comp_virial.eq.1) then
      call calc_U_dU(rho(i),params,U,dU)          ! compute embedding energy
                                                  !   and its derivative
      derU(i) = dU                                ! store dU for later use
   else
      call calc_U(rho(i),params,U)                ! compute just embedding energy
   endif

   ! accumulate the embedding energy contribution
   !
   ! Assuming U(rho=0) = 0.0_cd
   !
   if (comp_enepot.eq.1) then                     ! accumulate embedding energy contribution
      ene_pot(i) = U
   endif
   if (comp_energy.eq.1) then
      energy = energy + U
   endif

   if ((HalfOrFull.eq.1) .and. (i .gt. numberContrib)) exit
enddo

!  Loop over particles in the neighbor list a second time, to compute
!  the forces and complete energy calculation
!

! Reset iterator if one is being used
!
if (IterOrLoca.eq.1) then
   ier = kim_api_get_neigh(pkim,0,0,atom_ret,numnei,pnei1atom,pRij_list)
   if (ier.lt.KIM_STATUS_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "kim_api_get_neigh", ier)
      goto 42
   endif
endif

i = 0
do

   ! Set up neighbor list for next atom for all NBC methods
   !
   call get_current_atom_neighbors(IterOrLoca,NBC,N,pkim,      &
                                   i,numnei,nei1atom,Rij_list,ier)

   if (ier.eq.KIM_STATUS_NEIGH_ITER_PAST_END) exit  ! atom counter incremented past end of list
   if (ier.lt.KIM_STATUS_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "get_current_atom_neighbors", ier)
      goto 42
   endif


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
         where ( abs(Rij) .gt. 0.5_cd*boxSideLengths )  ! periodic boundary conditions
            Rij = Rij - sign(boxSideLengths,Rij)       ! applied where needed.
         end where                                !
      endif

      ! compute energy and forces
      !
      Rsqij = dot_product(Rij,Rij)                ! compute square distance
      if ( Rsqij .lt. params%model_cutsq ) then   ! particles are interacting?

         r = sqrt(Rsqij)                          ! compute distance
         if (comp_force.eq.1.or.comp_virial.eq.1) then
            call calc_phi_dphi(r,params,phi,dphi) ! compute pair potential
                                                  !   and it derivative
            call calc_dg(r,params,dg)             ! compute elect dens first deriv
            if ((HalfOrFull.eq.1) .and. &
                (j .le. numberContrib)) then      ! HALF mode
               dphii  = 0.5_cd*dphi
               dphij  = 0.5_cd*dphi
               dUi    = derU(i)*dg
               dUj    = derU(j)*dg
            else                                  ! FULL mode
               dphii  = 0.5_cd*dphi
               dphij  = 0.0_cd
               dUi    = derU(i)*dg
               dUj    = 0.0_cd
            endif
            dphieff = dphii + dphij + dUi + dUj
         else
            call calc_phi(r,params,phi)           ! compute just pair potential
         endif
         if ((HalfOrFull.eq.1) .and. &
             (j .le. numberContrib)) then         ! HALF mode
            Ei     = 0.5_cd*phi
            Ej     = 0.5_cd*phi
         else                                  ! FULL mode
            Ei     = 0.5_cd*phi
            Ej     = 0.0_cd
         endif

         ! contribution to energy
         !
         if (comp_enepot.eq.1) then
            ene_pot(i) = ene_pot(i) + Ei          ! accumulate energy Ei
            ene_pot(j) = ene_pot(j) + Ej          ! accumulate energy Ej
         endif
         if (comp_energy.eq.1) then
            energy = energy + Ei                  ! accumulate energy
            energy = energy + Ej                  ! accumulate energy
         endif

         ! contribution to virial tensor
         !
         if (comp_virial.eq.1) then
            virial_global(1) = virial_global(1) + Rij(1)*Rij(1)*dphieff/r
            virial_global(2) = virial_global(2) + Rij(2)*Rij(2)*dphieff/r
            virial_global(3) = virial_global(3) + Rij(3)*Rij(3)*dphieff/r
            virial_global(4) = virial_global(4) + Rij(2)*Rij(3)*dphieff/r
            virial_global(5) = virial_global(5) + Rij(1)*Rij(3)*dphieff/r
            virial_global(6) = virial_global(6) + Rij(1)*Rij(2)*dphieff/r
         endif

         ! contribution to forces
         !
         if (comp_force.eq.1) then                        ! Ei contribution
            force(:,i) = force(:,i) + dphieff*Rij/r ! accumulate force on atom i
            force(:,j) = force(:,j) - dphieff*Rij/r ! accumulate force on atom j
         endif

      endif

   enddo  ! loop on jj

enddo  ! infinite do loop (terminated by exit statements above)

! Free temporary storage
!
if (NBC.eq.0) deallocate( nei1atom_substitute )
deallocate( rho )
if (comp_force.eq.1.or.comp_virial.eq.1) deallocate( derU )

! Everything is great
!
ier = KIM_STATUS_OK
42 continue
Compute_Energy_Forces = ier
return

end function Compute_Energy_Forces

!-------------------------------------------------------------------------------
!
! Get list of neighbors for current atom using all NBC methods
!
!-------------------------------------------------------------------------------
subroutine get_current_atom_neighbors(IterOrLoca,NBC,N,pkim,      &
                                      atom,numnei,nei1atom,Rij_list,ier)
implicit none

!-- Transferred variables
integer(c_int),                  intent(in)    :: IterOrLoca
integer(c_int),                  intent(in)    :: NBC
integer(c_int),                  intent(in)    :: N
type(c_ptr),                     intent(in)    :: pkim
integer(c_int),                  intent(inout) :: atom
integer(c_int),                  intent(out)   :: numnei
integer(c_int),                  intent(out)   :: ier
integer(c_int), pointer,         intent(inout) :: nei1atom(:);   type(c_ptr) pnei1atom
real(c_double), pointer,         intent(inout) :: Rij_list(:,:); type(c_ptr) pRij_list

!-- Local variables
integer(c_int) atom_ret, jj
integer(c_int) idum

! Set up neighbor list for next atom for all NBC methods
!
if (IterOrLoca.eq.1) then    ! ITERATOR mode
   ier = kim_api_get_neigh(pkim,0,1,atom_ret,numnei, &
                           pnei1atom,pRij_list)
   if (ier.eq.KIM_STATUS_NEIGH_ITER_PAST_END) then
                          ! past end of the list, terminate loop in
      return              ! calling routine
   endif
   if (ier.lt.KIM_STATUS_OK) then     ! some sort of problem, exit
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "kim_api_get_neigh", ier)
      return
   endif
   call c_f_pointer(pnei1atom, nei1atom, [numnei])
   call c_f_pointer(pRij_list, Rij_list, [DIM,numnei])
   atom = atom_ret
else                         ! LOCATOR mode

   atom = atom + 1
   if (atom.gt.N) then                     ! incremented past end of list,
      ier = KIM_STATUS_NEIGH_ITER_PAST_END ! terminate loop in calling routine
      return
   endif

   if (NBC.eq.0) then ! CLUSTER NBC method
      numnei = N - atom   ! number of neighbors in list atom+1, ..., N
      nei1atom(1:numnei) = (/ (atom+jj, jj = 1,numnei) /)
      ier = KIM_STATUS_OK
   else
      ier = kim_api_get_neigh(pkim,1,atom,atom_ret,numnei,pnei1atom,pRij_list)
      if (ier.ne.KIM_STATUS_OK) then ! some sort of problem, exit
         idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                     "kim_api_get_neigh", ier)
         ier = KIM_STATUS_FAIL
         return
      endif
      call c_f_pointer(pnei1atom, nei1atom, [numnei])
      call c_f_pointer(pRij_list, Rij_list, [DIM,numnei])
   endif
endif

return
end subroutine get_current_atom_neighbors
!-------------------------------------------------------------------------------
!
! Model driver destroy routine
!
!-------------------------------------------------------------------------------
integer(c_int) function destroy(pkim) bind(c)
use KIM_API_F03
implicit none

!-- Transferred variables
type(c_ptr), intent(in) :: pkim

!-- Local variables
integer(c_int) :: ier, idum
type(model_parameters), pointer :: params; type(c_ptr) :: pparams

pparams = kim_api_get_model_buffer(pkim, ier)
if (ier.ne.KIM_STATUS_OK) then
  idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
    "kim_api_get_model_buffer", ier)
  destroy = ier
  return
endif
call c_f_pointer(pparams, params)

deallocate( params%vknotcoeff, params%vknotpoint )
deallocate( params%pknotcoeff, params%pknotpoint )

deallocate( params )

destroy = KIM_STATUS_OK
return

end function destroy

end module model_driver_pf_cubic_splines

!-------------------------------------------------------------------------------
!
! Model driver initialization routine (REQUIRED)
!
!-------------------------------------------------------------------------------
integer(c_int) function model_driver_init(pkim, pparamfile, nmstrlen, numparamfiles) &
  bind(c)
use, intrinsic :: iso_c_binding
use model_driver_pf_cubic_splines
use KIM_API_F03
implicit none

!-- Transferred variables
type(c_ptr),                     intent(in) :: pkim
type(c_ptr), value,              intent(in) :: pparamfile
integer(c_int),                  intent(in) :: nmstrlen
integer(c_int),                  intent(in) :: numparamfiles
character(len=nmstrlen), pointer            :: paramfile(:)

!-- Local variables
integer(c_int), parameter :: cd = c_double  ! used for literal constants
integer(c_int), parameter :: one=1
integer(c_int) i, j, ier, idum, return_error
character (LEN=55) :: text
logical ldum
type(model_parameters), pointer :: params

! define variables for all model parameters to be read in

!-- KIM variables
real(c_double), pointer :: cutoff; type(c_ptr) :: pcutoff

! assume all is well
!
return_error = KIM_STATUS_OK

call c_f_pointer(pparamfile, paramfile, [numparamfiles])

! find first null character and write spaces everywhere afterward
do i = 1, numparamfiles
  j = index(paramfile(1),char(0))
  paramfile(i)(j:)=" "
end do

! store function pointers in KIM object
call kim_api_setm_method(pkim, ier, &
     "compute", one, c_funloc(Compute_Energy_Forces), 1, &
     "destroy", one, c_funloc(destroy),               1)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_setm_method", ier)
   return_error = ier
   goto 1000
endif

! allocate parameters and register in KIM object
allocate( params )
call kim_api_set_model_buffer(pkim, c_loc(params), ier)
if (ier.ne.KIM_STATUS_OK) then
  idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
    "kim_api_set_model_buffer", ier)
  return_error = ier
  goto 1000
endif

! Read in model parameters from parameter file
open(10,file=paramfile(1),status="old")
    read(10,*,iostat=ier,err=100) params%model_cutoff
    read(10,*,iostat=ier,err=100) params%Z
    read(10,*,iostat=ier,err=100) params%r1, params%r2, &
      params%linear_interpolate, params%interpolate_num, params%a_rho
    read(10,*,iostat=ier,err=100) params%nknotp, params%nknotv
    !
    ! for A,B, and knot functions there are numbers and logicals to skip
    read(10,*,iostat=ier,err=100) idum, params%A, ldum
    !
    !
    read(10,*,iostat=ier,err=100) idum, params%B, ldum

    allocate (params%pknotcoeff(params%nknotp),params%pknotpoint(params%nknotp))
    allocate (params%vknotcoeff(params%nknotv),params%vknotpoint(params%nknotv))

    i=1
    DO WHILE ((ier==0).AND.(i.LE.params%nknotp))
     read(10,*,iostat=ier,err=100) idum, params%pknotcoeff(i), ldum, &
       params%pknotpoint(i)
     i=i+1
    END DO
    i=1
    DO WHILE ((ier==0).AND.(i.LE.params%nknotv))
     read(10,*,iostat=ier,err=100) idum, params%vknotcoeff(i), ldum, &
       params%vknotpoint(i)
     i=i+1
    END DO
    i=1
    params%a_inter=0.0_cd
    DO WHILE ((ier==0).AND.(i.LE.params%interpolate_num))
     read(10,*,iostat=ier,err=100) params%a_inter(i)
     i=i+1
    END DO
    !PRINT *,Z,r1,r2,A,B,a_rho,nknotp,nknotv,a_inter,model_cutoff,interpolate_num
    !PRINT *,vknotpoint,vknotcoeff,pknotpoint,pknotcoeff
    read(10,'(A55)',iostat=ier,err=100) text
close(10)


    write (*,*) 'Using ',TRIM(ADJUSTL(text))
    write (*,*) 'atomic number is: ',params%Z

goto 200
100 continue
! reading parameters failed
idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                            "Unable to read potential parameters, kimerror = ",KIM_STATUS_FAIL)
return_error = KIM_STATUS_FAIL
goto 1000
!PRINT *,Z,r1,r2,A,B,a_rho,nknotp,nknotv,a_inter,model_cutoff,interpolate_num
!PRINT *,vknotpoint,vknotcoeff,pknotpoint,pknotcoeff

200 continue

! store model cutoff in KIM object
pcutoff =  kim_api_get_data(pkim,"cutoff",ier)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_get_data", ier)
   return_error = ier
   goto 1000
endif
call c_f_pointer(pcutoff, cutoff)
cutoff = params%model_cutoff
params%model_cutsq   = params%model_cutoff**2

1000 continue
model_driver_init = return_error
return

end function model_driver_init

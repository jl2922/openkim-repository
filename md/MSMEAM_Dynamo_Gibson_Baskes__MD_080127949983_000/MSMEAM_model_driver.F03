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
!    Valeriu Smirichinski
!    Stephen M. Whalen
!

!
! Portions Copyright (c) 2015, Joshua S. Gibson.
! All rights reserved.
!
! Contributors:
!    Joshua S. Gibson
!    Srinivasan G. Srivilliputhur
!    Michael I. Baskes
!    Ronald E. Miller
!    Angela K. Wilson
!

!****************************************************************************
!**
!**  MODULE MSMEAM
!**
!**  MSMEAM model driver
!**
!**  Reference: [PENDING]
!**
!**  Language: Fortran 2003
!**
!****************************************************************************


#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__
#define TRUEFALSE(TRUTH) merge(1,0,(TRUTH))

module MSMEAM

use, intrinsic :: iso_c_binding
use KIM_API_F03
implicit none

save
private
public Compute_Energy_Forces, &
       model_parameters,      &
       destroy,               &
       elecden,               &
       rscrn,                 &
       get_c,                 &
       calc_phi,              &
       rhof,                  &
       drhof,                 &
       dphif,                 &
       frhoi,                 &
       dfrhoi,                &
       get_s,                 &
       dcmij,                 &
       lagrangefxn,           &
       lagrangedrv

type model_parameters
   integer(c_int)                                    :: NBC
   integer(c_int)                                    :: HalfOrFull
   integer(c_int)                                    :: IterOrLoca
   integer(c_int)                                    :: nel
   integer(c_int)                                    :: ntab
   real(c_double)                                    :: rhotabmin
   real(c_double)                                    :: rhotabmax
   real(c_double)                                    :: rtabmax
   real(c_double)                                    :: rrtabmax
   real(c_double)                                    :: model_cutoff
   real(c_double)                                    :: model_cutsq
   real(c_double), pointer, dimension(:)             :: asubs
   real(c_double), pointer, dimension(:)             :: bsubs
   real(c_double), pointer, dimension(:,:)           :: esubs
   real(c_double), pointer, dimension(:,:)           :: res
   real(c_double), pointer, dimension(:,:)           :: embedtab
   real(c_double), pointer, dimension(:,:,:)         :: phitab
   real(c_double), pointer, dimension(:,:,:,:)       :: rhotab
   real(c_double), pointer, dimension(:,:,:,:)       :: scrp
   real(c_double), pointer, dimension(:,:,:,:)       :: scr0
   real(c_double), pointer, dimension(:,:,:,:)       :: scr1
   real(c_double), pointer, dimension(:,:,:,:)       :: scr2
   real(c_double), pointer, dimension(:,:,:,:)       :: scr3
   real(c_double), pointer, dimension(:,:,:,:)       :: scr54
end type model_parameters


! Below are the definitions and values of all model parameters
integer(c_int), parameter :: cd = c_double  ! used for literal constants
integer(c_int), parameter :: DIM = 3        ! dimensionality of space
integer(c_int), parameter :: speccode = 1   ! internal species code

!-------------------------------------------------------------------------------
! Below are the definitions and values of all additional model parameters
!
! Recall that the Fortran 2003 format for declaring parameters is as follows:
!
! integer(c_int), parameter :: parname = value ! This defines an
!                                              ! integer(c_int) parameter
!                                              ! called `parname' with a value
!                                              ! equal to `value' (a number)
!
! real(c_double), parameter :: parname = value ! This defines a real(c_double)
!                                              ! parameter called `parname'
!                                              ! with a value equal to
!                                              ! `value' (a number)
!-------------------------------------------------------------------------------
integer(c_int), parameter ::  natmax=10000
integer(c_int), parameter ::  nelmax=2
integer(c_int), parameter ::  ntable=10001
integer(c_int), parameter ::  ifree=0
integer(c_int), parameter ::  noscr=0
real(c_double), parameter ::  scnres=0.0_cd
real(c_double), parameter ::  legend00=3./5._cd
real(c_double), parameter ::  legend01=10./9._cd
real(c_double), parameter ::  legend02=5./21._cd

contains

!-------------------------------------------------------------------------------
!
!  Calculate lagrange interpolation
!
!-------------------------------------------------------------------------------
real(c_double) function lagrangefxn ( p, yd )
implicit none

integer(c_int), parameter                    :: n = 4
real(c_double)                               :: lazy
real(c_double), intent(in)                   :: p
real(c_double), intent(in), dimension(1:n)   :: yd
real(c_double), dimension(1:n)               :: cof
integer(c_int)                               :: j

lazy=0.0_cd

cof(1) = -0.166666667_cd*(p)*(p-1.0_cd)*(p-2.0_cd)
cof(2) = 0.5_cd*(p+1.0_cd)*(p-1.0_cd)*(p-2.0_cd)
cof(3) = -0.5_cd*(p+1.0_cd)*(p)*(p-2.0_cd)
cof(4) = 0.166666667_cd*(p)*(p+1.0_cd)*(p-1.0_cd)

do j=1,4
   lazy=lazy+cof(j)*yd(j)
enddo

lagrangefxn=lazy

return
end function lagrangefxn

!-------------------------------------------------------------------------------
!
!  Calculate derivative of lagrange interpolation
!
!-------------------------------------------------------------------------------
real(c_double) function lagrangedrv ( xt, yd )
implicit none

integer(c_int), parameter                    :: ione = 1
integer(c_int), parameter                    :: n = 4
real(c_double), intent(in)                   :: xt
real(c_double)                               :: lazy
real(c_double), intent(in), dimension(1:n)   :: yd
real(c_double), parameter, dimension(1:n)    :: xd = &
 & [-1,0,1,2]
real(c_double), parameter, dimension(1:n)    :: q = &
 & [-6,2,-2,6]

real(c_double), dimension(1:n)               :: d,pa,pb,dpa,dpb,dw
integer(c_int)                               :: j

pa(1)=1.0_cd
dpa(1)=0.0_cd
do j=1,n-ione
   d(j)=xt-xd(j)
   pa(j+ione)=pa (j)*d(j)
   dpa(j+ione)=dpa(j)*d(j)+pa(j)
enddo
d(n)=xt-xd(n)

pb(n)=1.0_cd
dpb(n)=0.0_cd
do j=n,2,-1
   pb (j-ione)=pb(j)*d(j)
   dpb(j-ione)=dpb(j)*d(j)+pb(j)
enddo

lazy=0.0_cd

do j=1,n
   dw(j)=(pa(j)*dpb(j)+dpa(j)*pb(j))/q(j)
   lazy=lazy+dw(j)*yd(j)
enddo

lagrangedrv=lazy

return
end function lagrangedrv

!-------------------------------------------------------------------------------
!
!  Calculate screening functions
!
!-------------------------------------------------------------------------------
subroutine get_s(cx,it,jt,kt,meatable,sp,s0,s1,s2,s3,s54)

real(c_double), intent(in)          :: cx
integer(c_int), intent(in)          :: it,jt,kt
type(model_parameters), intent(in)  :: meatable

real(c_double), intent(out)         :: sp,s0,s1,s2,s3,s54

real(c_double)                      :: p
real(c_double), parameter           :: ctabmax=3.0_cd
integer(c_int)                      :: k,ctr
integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double),dimension(1:Np)      :: sparray,s0array,s1array
real(c_double),dimension(1:Np)      :: s2array,s3array,s5array

if(cx.ge.ctabmax) then
     sp=1.0_cd
     s0=1.0_cd
     s1=1.0_cd
     s2=1.0_cd
     s3=1.0_cd
     s54=1.0_cd
     return
endif

p = cx * ( dble(meatable%ntab) - 1.0_cd ) / ctabmax + 1.0_cd
k = int(p)

if ( k .gt. meatable%ntab - 1 ) then
   sp = 0.0_cd
   s0 = 0.0_cd
   s1 = 0.0_cd
   s2 = 0.0_cd
   s3 = 0.0_cd
   s54 = 0.0_cd
   return
endif

k = min( k , meatable%ntab - Pr )
k = max( k , Pr )
p = p - dble( k )

!   make sure that p is less than Pr
!   then if cx is out of range, Pr and c = last value of screening
!

p = min( p ,dble(Pr) )
      
do ctr=1,Np
   sparray(ctr)=meatable%scrp( k - Pr + ctr , it, jt, kt)
   s0array(ctr)=meatable%scr0( k - Pr + ctr , it, jt, kt)
   s1array(ctr)=meatable%scr1( k - Pr + ctr , it, jt, kt)
   s2array(ctr)=meatable%scr2( k - Pr + ctr , it, jt, kt)
   s3array(ctr)=meatable%scr3( k - Pr + ctr , it, jt, kt)
   s5array(ctr)=meatable%scr54( k - Pr + ctr , it, jt, kt)
enddo

sp=lagrangefxn(p,sparray)
s0=lagrangefxn(p,s0array)
s1=lagrangefxn(p,s1array)
s2=lagrangefxn(p,s2array)
s3=lagrangefxn(p,s3array)
s54=lagrangefxn(p,s5array)

sp = max( 0.0_cd , sp )
sp = min( 1.0_cd , sp )

s0 = max( 0.0_cd , s0 )
s0 = min( 1.0_cd , s0 )

s1 = max( 0.0_cd , s1 )
s1 = min( 1.0_cd , s1 )

s2 = max( 0.0_cd , s2 )
s2 = min( 1.0_cd , s2 )

s3 = max( 0.0_cd , s3 )
s3 = min( 1.0_cd , s3 )

s54 = max( 0.0_cd , s54 )
s54 = min( 1.0_cd , s54 )

return
end subroutine get_s
!-------------------------------------------------------------------------------
!
!  Calculate background electron density
!
!-------------------------------------------------------------------------------
subroutine elecden(rho0,AA,dang1,dang2,rhobar)
implicit none

real(c_double)                 :: temp
real(c_double),intent(in)      :: rho0,AA
real(c_double),intent(out)     :: dang1,dang2,rhobar

!   rhobar^2 = rho0^2+AA
!   dang is derivative of rhobar wrt rho0 and  AA
    temp = rho0**2 + AA
if(temp.le.0.0_cd) then
    rhobar=1d-20
    dang1=0.0_cd
    dang2=0.0_cd
  else
    rhobar=sqrt(temp)
    dang1= rho0/rhobar
    dang2= 0.5_cd/rhobar
endif
return
end subroutine elecden
!-------------------------------------------------------------------------------
!
!  Calculate radial screening
!
!-------------------------------------------------------------------------------
subroutine rscrn(rij,meatable,scrn)
!   calculate r cut-off function
implicit none

type(model_parameters), intent(in)  :: meatable

real(c_double), intent(in) :: rij

real(c_double), intent(out) :: scrn

real(c_double)           :: fac,rcutsq
real(c_double),parameter :: frcut=0.9
real(c_double),parameter :: xncut=2.0
real(c_double),parameter :: xmcut=4.0

rcutsq=meatable%model_cutsq

fac=(rij-frcut*rcutsq)/((1.0-frcut)*rcutsq)
if(fac .le. 0.0_cd) then
    scrn=1.0_cd
elseif(fac .ge. 1.0_cd) then
    scrn=0.0_cd
else
    scrn=(1.0_cd-(fac**xmcut))**xncut
endif
return
end subroutine rscrn
!-------------------------------------------------------------------------------
!
!  Calculate screening parameter
!
!-------------------------------------------------------------------------------
real(c_double) function get_c (rij,ril,rjl)
implicit none

real(c_double), intent(in) :: rij,ril,rjl
real(c_double) :: dotil,dotjl,argb,argt,Xi,Xj

real(c_double),parameter :: sconst=1.125

get_c=3.0_cd
if ((ril.gt.sconst*rij).or.(rjl.gt.sconst*rij)) return
!
!   If l atom is outside of ellipse there is no screening
!   dotx = 2 r_ij dot r_x
!   this eliminates screeners with distances greater than 0.2*rij
!
dotil=ril+rij-rjl
if(dotil.le.0.027_cd*rij) return
dotjl=rij+rjl-ril
if(dotjl.le.0.027_cd*rij) return

Xi=(ril/rij)
Xj=(rjl/rij)
argb=1.0_cd-(Xi-Xj)**2
argt=2.0_cd*(Xi+Xj)-1.0_cd-(Xi-Xj)**2
get_c=argt/argb

return
end function get_c
!-------------------------------------------------------------------------------
!
!  Calculate pair potential phi(i) for atom i
!
!-------------------------------------------------------------------------------
real(c_double) function calc_phi(r,it,jt,meatable)
implicit none

!-- Transferred variables
real(c_double), intent(in)          :: r
type(model_parameters), intent(in)  :: meatable
integer(c_int), intent(in)          :: it,jt

!-- Local variables
integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double)                      :: p
real(c_double),dimension(1:Np)      :: phiarray
integer(c_int)                      :: k,ctr

p= r * ( meatable%ntab - 1.0_cd ) / meatable%rrtabmax + 1.0_cd
k = int(p)

if ( k .gt. meatable%ntab - 1 ) then
     calc_phi = 0.0_cd
     return
endif

k = min(k,meatable%ntab - Pr)
k = max(k,Pr)
p = p - dble(k)

!       make sure that p is less than 2.0
!       then if r is out of range, p = 2.0 and r = last value of phitab
!

p = min(p,dble(Pr))

do ctr=1,Np
   phiarray(ctr)= meatable%phitab( k - Pr + ctr , it, jt )
enddo

calc_phi = lagrangefxn(p,phiarray)

return
end function calc_phi
!-------------------------------------------------------------------------------
!
!  Calculate rhol+ and rhol- at a given r with atom j
!
!-------------------------------------------------------------------------------
real(c_double) function rhof(r,l,ispin,it,meatable)
implicit none

!-- Transferred variables
type(model_parameters), intent(in)  :: meatable
integer(c_int),intent(in)           :: l
integer(c_int),intent(in)           :: ispin
integer(c_int),intent(in)           :: it
real(c_double),intent(in)           :: r

!-- Local variables
integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double)                      :: p
real(c_double),dimension(1:Np)      :: rhoarray
integer(c_int)                      :: k,ctr

p= r * ( meatable%ntab - 1.0_cd ) / meatable%rrtabmax + 1.0_cd
k = int(p)

if ( k .gt. meatable%ntab - 1 ) then
     rhof = 0.0_cd
     return
endif

k = min ( k , meatable%ntab - Pr )
k = max ( k , Pr )
p = p - dble(k)

!       make sure that p is less than 2.0
!       then if r is out of range, p = 2.0 and r = last value of phitab
!

p = min ( p , dble(Pr))

do ctr=1,Np
   rhoarray(ctr)=meatable%rhotab( k - Pr + ctr , l , ispin , it )
enddo

rhof = lagrangefxn(p,rhoarray)

return
end function rhof
!-------------------------------------------------------------------------------
!
!  Calculate the derivatives of rhol+ and rhol- at a given r with atom j
!
!-------------------------------------------------------------------------------
real(c_double) function drhof(r,l,ispin,it,meatable)
implicit none

integer(c_int),intent(in)           :: l
integer(c_int),intent(in)           :: ispin
integer(c_int),intent(in)           :: it
real(c_double),intent(in)           :: r
type(model_parameters), intent(in)  :: meatable

integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double)                      :: p
real(c_double),dimension(1:Np)      :: rhoarray
integer(c_int)                      :: k,ctr

p = r * ( dble( meatable%ntab ) - 1.0_cd ) / meatable%rrtabmax + 1.0_cd
k = int(p)

if ( k .gt. meatable%ntab-1) then
     drhof = 0.0_cd
     return
endif

k = min( k , meatable%ntab - Pr )
k = max( k , Pr )
p = p - dble( k )

!       make sure that p is less than 2.0
!       then if r is out of range, p = 2.0 and r = last value of rhotab
!

p = min( p , dble(Pr) )

do ctr=1,Np
   rhoarray(ctr)=meatable%rhotab( k - Pr + ctr, l , ispin , it )
enddo

drhof = lagrangedrv(p,rhoarray) * ( dble( meatable%ntab ) - 1.0_cd ) / meatable%rrtabmax

return
end function drhof
!-------------------------------------------------------------------------------
!
!  Calculate the derivative of the pair potential dphi(r)
!
!-------------------------------------------------------------------------------
real(c_double) function dphif(r,it,jt,meatable)
implicit none

!-- Transferred variables
real(c_double), intent(in)          :: r
integer(c_int), intent(in)          :: it,jt
type(model_parameters), intent(in)  :: meatable

!-- Local variables
integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double)                      :: p
real(c_double),dimension(1:Np)      :: phiarray
integer(c_int)                      :: k,ctr

!
!  do four-point lagrange interpolation
!

p = r * ( dble( meatable%ntab ) - 1.0_cd ) / meatable%rrtabmax  + 1.0_cd
k = int(p)

if (k .gt. meatable%ntab-1 ) then
     dphif = 0.0_cd
     return
endif

k = min( k , meatable%ntab - Pr )
k = max( k , Pr )
p = p - dble( k )

!       make sure that p is less than 3.0
!       then if r is out of range, p = 3.0 and r = last value of phitab
!

p = min( p , dble(Pr) )

do ctr=1,Np
   phiarray(ctr)= meatable%phitab( k - Pr + ctr , it, jt )
enddo

dphif = lagrangedrv(p,phiarray) * ( dble( meatable%ntab ) - 1.0_cd ) / meatable%rrtabmax

return
end function dphif
!-------------------------------------------------------------------------------
!
!  Calculate embedding function for atom i
!
!-------------------------------------------------------------------------------
real(c_double) function frhoi(rho,it,asub,bsub,esub,meatable)
implicit none

!-- Transferred variables
real(c_double), intent(in)          :: rho,asub,bsub,esub
integer(c_int), intent(in)          :: it
type(model_parameters), intent(in)  :: meatable

!-- Local variables
integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double)                      :: xrho,p
integer(c_int)                      :: k,ctr
real(c_double),dimension(1:Np)      :: embarray

if ( rho .le. meatable%rhotabmin ) then
     frhoi = esub * ( asub * log( rho ) + bsub ) * rho
else
     xrho = log ( rho / meatable%rhotabmin )
     p= xrho * ( dble(meatable%ntab) - 1.0_cd ) / meatable%rhotabmax + 1.0_cd
     k = int(p)
     k = min ( k, meatable%ntab - Pr )
     k = max ( k , Pr )
     p = p - dble( k )

     p = min( p , dble(Pr) )

   do ctr=1,Np
      embarray(ctr)=meatable%embedtab( k - Pr + ctr , it ) 
   enddo

frhoi = lagrangefxn(p,embarray)

endif
return
end function frhoi

!-------------------------------------------------------------------------------
!
!  Calculate first derivative of embedding function for atom i
!
!-------------------------------------------------------------------------------
real(c_double) function dfrhoi(rho,it,asub,bsub,esub,meatable)
implicit none

!-- Transferred variables
real(c_double), intent(in)          :: rho,asub,bsub,esub
integer(c_int), intent(in)          :: it
type(model_parameters), intent(in)  :: meatable

!-- Local variables
integer(c_int), parameter           :: Np=4
integer(c_int), parameter           :: Pr=Np/2
real(c_double)                      :: xrhoi,p
integer(c_int)                      :: k,ctr
real(c_double),dimension(1:Np)      :: embarray

if( rho .lt. meatable%rhotabmin  ) then
     dfrhoi = esub * ( asub * log( rho ) + bsub + asub )
else
     xrhoi = log ( rho / meatable%rhotabmin  )
     p= xrhoi * ( dble( meatable%ntab ) - 1.0_cd ) / meatable%rhotabmax  + 1.0_cd
     k = int(p)
     k = min( k , meatable%ntab - Pr )
     k = max( k , Pr )
     p = p - dble( k )

!       make sure that p is less than 1.0
!       then if r is out of range, p = 1.0 and rho = last value of embedtab
!

     p = min( p , dble(Pr) )

     do ctr=1,Np
        embarray(ctr)=meatable%embedtab( k - Pr + ctr , it )
     enddo

     dfrhoi = lagrangedrv(p,embarray)*( dble( meatable%ntab ) - 1.0_cd ) / meatable%rhotabmax 
     dfrhoi = dfrhoi / max( rho , meatable%rhotabmin )

endif
return
end function dfrhoi
!-------------------------------------------------------------------------------
!
!  Compute derivative of vector components
!
!-------------------------------------------------------------------------------
subroutine dcmij(drcos,rr,rs)

implicit none

!-- Transferred variables
real(c_double), dimension(1:3,1:3), intent(out)   :: drcos
real(c_double), dimension(1:3), intent(in)        :: rr
real(c_double), intent(in)                        :: rs

!-- Local variables
real(c_double) :: r1,r13,tmp1
integer(c_int) :: i,j

!
!     Find direction cosines
!
        r1 = 1.0_cd / rs
        r13 = r1 * r1 * r1
        do i=1,3
           tmp1 = r13 * rr(i)
           do j=1,3
              drcos(i,j) = -tmp1 * rr(j)
           enddo
           drcos(i,i) = drcos(i,i) + r1
        enddo

end subroutine dcmij
!-------------------------------------------------------------------------------
!
! Compute energy and forces on particles from the positions.
!
!-------------------------------------------------------------------------------
integer(c_int) function Compute_Energy_Forces(pkim) bind(c)
implicit none

!-- Transferred variables
type(c_ptr), intent(in) :: pkim

!-- Local variables
integer(c_int)                                     :: comp_force,comp_energy
integer(c_int)                                     :: comp_enepot,comp_virial
integer(c_int)                                     :: comp_ptvirl
integer(c_int)                                     :: ispin,kk,l,m,k,s,w,ier
integer(c_int)                                     :: i,j,jj,ll,mm,mc,atom
integer(c_int)                                     :: numnei,part_ret
integer(c_int)                                     :: neighOnI

integer(c_int), parameter                          :: ity=1
integer(c_int), parameter                          :: jty=1
integer(c_int), parameter                          :: kty=1
integer(c_int), parameter                          :: neimax=2000
real(c_double), parameter                          :: h=1.0e-5_cd
real(c_double), parameter                          :: hh=h*h
real(c_double), parameter                          :: h2 = 2.0_cd * h
integer(c_int), parameter                          :: i2jl=1
integer(c_int), parameter                          :: j2il=2

real(c_double)          :: r_ij,r1,rej,rs,re
real(c_double)          :: negRpij,posRpij,negRpji,posRpji
real(c_double)          :: Rsqij,Rsqil,Rsqjl,Rsqim
real(c_double)          :: Rpjl,Rsqji,Rpji,Rpij
real(c_double)          :: Rsqmj,Rpim,Rpmj,Rpil
real(c_double)          :: dxij,dyij,dzij
real(c_double)          :: dxim,dyim,dzim
real(c_double)          :: dxmj,dymj,dzmj
real(c_double)          :: rho1p,rho2p,rho3p,rho54p
real(c_double)          :: rho1m,rho2m,rho3m,rho54m
real(c_double)          :: cp,cm,cgp,cgm
real(c_double)          :: AA,asub,bsub,esub
real(c_double)          :: cx,t1,t2,t3,t54
real(c_double)          :: dang1,dang2
real(c_double)          :: rho0,phi,plocal,rhobar,roz
real(c_double)          :: screenijp,screenij0,screenij1
real(c_double)          :: screenij2,screenij3,screenij54
real(c_double)          :: sp,s0,s1,s2,s3,s54
real(c_double)          :: dpp,dp0,dp1,dp2,dp3,dp54
real(c_double)          :: dmp,dm0,dm1,dm2,dm3,dm54
real(c_double)          :: screep,scree0,scree1
real(c_double)          :: scree2,scree3,scree54
real(c_double)          :: spp,sp0,sp1,sp2,sp3,sp54
real(c_double)          :: smp,sm0,sm1,sm2,sm3,sm54
real(c_double)          :: smax,smin,term1,term2,term3
real(c_double)          :: angular1,angular2
real(c_double)          :: dfrho,phiif,phiij
real(c_double)          :: rcosj,rcosm,rcosl,rcoss,rcosk,rcosw
real(c_double)          :: rml,rmls,rmlsk,rmlskw
real(c_double)          :: drmllm,drmlsslm,drmlskkslm,drmlskwwkslm
real(c_double)          :: drho0,drho0s
real(c_double)          :: drho1,drho1g,drho1gh
real(c_double)          :: drho1s,drho1sg,drho1sgh
real(c_double)          :: drho2,drho2s
real(c_double)          :: drho3,drho3s,drho3h,drho3sh
real(c_double)          :: drho5,drho5s
real(c_double)          :: dscrnabkp,dscrnabk0,dscrnabk1
real(c_double)          :: dscrnabk2,dscrnabk3,dscrnabk54
real(c_double)          :: tmp0p,tmp0m
real(c_double)          :: tmp1p,tmp1m,tmp1gp,tmp1gm,tmp1hp,tmp1hm
real(c_double)          :: tmp2p,tmp2m
real(c_double)          :: tmp3p,tmp3m,tmp3hp,tmp3hm
real(c_double)          :: tmp5p,tmp5m
real(c_double)          :: tmp4,tmp8

real(c_double),dimension(1:3)                      :: Rij,Ril,Rjl,Rim,rcos,rr,Rji
real(c_double),dimension(1:3)                      :: dscrpp,dscrp0,dscrp1
real(c_double),dimension(1:3)                      :: dscrp2,dscrp3,dscrp54
real(c_double),dimension(1:3)                      :: dscrmp,dscrm0,dscrm1
real(c_double),dimension(1:3)                      :: dscrm2,dscrm3,dscrm54
real(c_double),dimension(1:3)                      :: ap,am,agp,agm,ahp,ahm
real(c_double),dimension(1:3)                      :: dfij,dfim,dfji
real(c_double),dimension(1:3,1:3)                  :: drcos,bp,bm
real(c_double),dimension(1:3,1:3,1:3)              :: dp,dm,dhp,dhm
real(c_double),dimension(1:3,1:3,1:3,1:3,1:3)      :: gp,gm
real(c_double),dimension(0:4)                      :: rhjp,rhjm,drhjp,drhjm

integer(c_int),dimension(1:neimax)                 :: jIdxOnI
real(c_double),dimension(1:neimax)                 :: scrnabp,scrnab0,scrnab1
real(c_double),dimension(1:neimax)                 :: scrnab2,scrnab3,scrnab54
real(c_double),dimension(1:2,1:3,1:neimax)         :: dscrnabp,dscrnab0,dscrnab1
real(c_double),dimension(1:2,1:3,1:neimax)         :: dscrnab2,dscrnab3,dscrnab54
real(c_double),dimension(1:3,1:neimax)             :: Rij_xyzN

real(c_double), allocatable, dimension(:,:,:)      :: slocal

!-- KIM variables
integer(c_int), pointer :: N;                   type(c_ptr) :: pN
real(c_double), pointer :: energy;              type(c_ptr) :: penergy
real(c_double), pointer :: coor(:,:);           type(c_ptr) :: pcoor
real(c_double), pointer :: force(:,:);          type(c_ptr) :: pforce
real(c_double), pointer :: enepot(:);           type(c_ptr) :: penepot
real(c_double), pointer :: Rij_list(:,:);       type(c_ptr) :: pRij_list
integer(c_int), pointer :: numContrib;          type(c_ptr) :: pnumContrib
integer(c_int), pointer :: nei1atom(:);         type(c_ptr) :: pnei1atom
integer(c_int), pointer :: particleSpecies(:);  type(c_ptr) :: pparticleSpecies
real(c_double), pointer :: particleVirial(:,:); type(c_ptr) :: pparticleVirial
real(c_double), pointer :: virial(:);           type(c_ptr) :: pvirial
integer(c_int) numberContrib
integer(c_int) idum
integer(c_int) atom_ret

type(model_parameters), pointer :: meatable; type(c_ptr) :: pmeatable

numberContrib = 0 ! initialize

! get model buffer from KIM object
pmeatable = kim_api_get_model_buffer(pkim, ier)
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_get_model_buffer", ier)
   goto 42
endif
call c_f_pointer(pmeatable, meatable)

! Check to see if we have been asked to compute the forces, energyperpart,
! energy and virial
!
call kim_api_getm_compute(pkim, ier, &
     "energy",         comp_energy,     1, &
     "forces",         comp_force,      1, &
     "particleEnergy", comp_enepot,     1, &
     "particleVirial", comp_ptvirl,     1, &
     "virial",         comp_virial,     1  &
      )
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_getm_compute", &
                               ier)
   goto 42
endif

! Unpack data from KIM object
!
call kim_api_getm_data(pkim, ier, &
     "numberOfParticles",           pN,              1,                                    &
     "particleSpecies",             pparticleSpecies,1,                                    &
     "coordinates",                 pcoor,           1,                                    &
     "numberContributingParticles", pnumContrib,     TRUEFALSE(meatable%HalfOrFull.eq.1),  &
     "energy",                      penergy,         TRUEFALSE(comp_energy.eq.1),          &
     "forces",                      pforce,          TRUEFALSE(comp_force.eq.1),           &
     "particleEnergy",              penepot,         TRUEFALSE(comp_enepot.eq.1),          &
     "particleVirial",              pparticleVirial, TRUEFALSE(comp_ptvirl.eq.1),          &
     "virial",                      pvirial,         TRUEFALSE(comp_virial.eq.1)           &
     )
if (ier.lt.KIM_STATUS_OK) then
   idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                               "kim_api_getm_data", ier)
   goto 42
endif

call c_f_pointer(pN, N)
call c_f_pointer(pparticleSpecies, particleSpecies, [N])
call c_f_pointer(pcoor, coor, [DIM,N])
if (meatable%HalfOrFull.eq.1)   call c_f_pointer(pnumContrib, numContrib)
if (comp_energy.eq.1)      call c_f_pointer(penergy, energy)
if (comp_force.eq.1)       call c_f_pointer(pforce,force, [DIM,N])
if (comp_enepot.eq.1)      call c_f_pointer(penepot,enepot, [N])
if (comp_ptvirl.eq.1)      call c_f_pointer(pparticleVirial, particleVirial, [6,N])
if (comp_virial.eq.1)      call c_f_pointer(pvirial, virial, [6])
if (meatable%HalfOrFull.eq.1)   then
   if (meatable%NBC.ne.0) then ! non-CLUSTER cases
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

!  Initialize energies, forces, and virials
!
if (comp_enepot.eq.1)     enepot(1:N)                = 0.0_cd
if (comp_energy.eq.1)     energy                     = 0.0_cd
if (comp_force.eq.1)      force(1:3,1:N)             = 0.0_cd
if (comp_ptvirl.eq.1)     particleVirial(1:6,1:N)    = 0.0_cd
if (comp_virial.eq.1)     virial(1:6)                = 0.0_cd

! Reset iterator if one is being used
!
if (meatable%IterOrLoca.eq.1) then
   ier = kim_api_get_neigh(pkim,0,0,atom_ret,numnei,pnei1atom,pRij_list)
   if (ier.lt.KIM_STATUS_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "kim_api_get_neigh", ier)
      goto 42
   endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Loop over particles in the neighbor list
!  to get neighbor information, screening, and screening derivatives,
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(slocal(1:3,1:3,1:N))
slocal(1:3,1:3,1:N) = 0.0_cd

do i = 1,N

   !
   ! Initialize arrays
   !

   neighOnI                     = 0
   dfij(1:3)                    = 0.0_cd
   dfji(1:3)                    = 0.0_cd
   dfim(1:3)                    = 0.0_cd
   Rij_xyzN(1:3,1:neimax)       = 0.0_cd
   jIdxOnI(1:neimax)            = 0
   scrnabp(1:neimax)            = 0.0_cd
   scrnab0(1:neimax)            = 0.0_cd
   scrnab1(1:neimax)            = 0.0_cd
   scrnab2(1:neimax)            = 0.0_cd
   scrnab3(1:neimax)            = 0.0_cd
   scrnab54(1:neimax)           = 0.0_cd
   dscrnabp(1:2,1:3,1:neimax)   = 0.0_cd
   dscrnab0(1,2:1:3,1:neimax)   = 0.0_cd
   dscrnab1(1:2,1:3,1:neimax)   = 0.0_cd
   dscrnab2(1:2,1:3,1:neimax)   = 0.0_cd
   dscrnab3(1:2,1:3,1:neimax)   = 0.0_cd
   dscrnab54(1:2,1:3,1:neimax)  = 0.0_cd

   rho0                     = 0.0_cd
   dang1                    = 0.0_cd
   dang2                    = 0.0_cd
   phi                      = 0.0_cd
   plocal                   = 0.0_cd
   rho1p                    = 0.0_cd
   rho2p                    = 0.0_cd
   rho3p                    = 0.0_cd
   rho54p                   = 0.0_cd
   rho1m                    = 0.0_cd
   rho2m                    = 0.0_cd
   rho3m                    = 0.0_cd
   rho54m                   = 0.0_cd
   cp                       = 0.0_cd
   cm                       = 0.0_cd
   cgp                      = 0.0_cd
   cgm                      = 0.0_cd
   angular1                 = 0.0_cd
   angular2                 = 0.0_cd
   rhjp(0:4)                = 0.0_cd
   rhjm(0:4)                = 0.0_cd
   drhjp(0:4)               = 0.0_cd
   drhjm(0:4)               = 0.0_cd
   ap(1:3)                  = 0.0_cd
   am(1:3)                  = 0.0_cd
   agp(1:3)                 = 0.0_cd
   agm(1:3)                 = 0.0_cd
   ahp(1:3)                 = 0.0_cd
   ahm(1:3)                 = 0.0_cd
   bp(1:3,1:3)              = 0.0_cd
   bm(1:3,1:3)              = 0.0_cd
   dp(1:3,1:3,1:3)          = 0.0_cd
   dm(1:3,1:3,1:3)          = 0.0_cd
   dhp(1:3,1:3,1:3)         = 0.0_cd
   dhm(1:3,1:3,1:3)         = 0.0_cd
   gp(1:3,1:3,1:3,1:3,1:3)  = 0.0_cd
   gm(1:3,1:3,1:3,1:3,1:3)  = 0.0_cd

   ier = kim_api_get_neigh(pkim,1,i,part_ret,numnei,pnei1atom,pRij_list)
   if (ier.ne.KIM_STATUS_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  "kim_api_get_neigh", ier)
      goto 42
   endif

   if (meatable%NBC.ne.3) call c_f_pointer(pnei1atom, nei1atom, [numnei])
   if (meatable%NBC.eq.2) call c_f_pointer(pRij_list, Rij_list, [DIM,numnei])

   if (numnei .gt. 0) then

      do jj = 1,numnei

         j = nei1atom(jj)

         Rij(:) = -(coor(:,j) - coor(:,i))
         Rsqij = dot_product(Rij,Rij)

         Rji(:) = -(coor(:,i) - coor(:,j))
         Rsqji = dot_product(Rji,Rji)

         if (Rsqij .ge. meatable%model_cutsq) cycle

         neighOnI = neighOnI + 1
         Rij_xyzN(:,neighOnI) = -Rij(:)
         jIdxOnI(neighOnI) = j

         call rscrn(Rsqij,meatable,screenijp)
         call rscrn(Rsqij,meatable,screenij0)
         call rscrn(Rsqij,meatable,screenij1)
         call rscrn(Rsqij,meatable,screenij2)
         call rscrn(Rsqij,meatable,screenij3)
         call rscrn(Rsqij,meatable,screenij54)

         do ll = 1,numnei

            if (ll .eq. jj) cycle

            l = nei1atom(ll)

            Ril(:) = -(coor(:,l) - coor(:,i))
            Rsqil = dot_product(Ril,Ril)

            Rjl(:) = -(coor(:,l) - coor(:,j))
            Rsqjl = dot_product(Rjl,Rjl)

            cx = get_c(Rsqij,Rsqil,Rsqjl)

            call get_s(cx,ity,jty,kty,meatable,sp,s0,s1,s2,s3,s54)
            screenijp   = screenijp * sp
            screenij0   = screenij0 * s0
            screenij1   = screenij1 * s1
            screenij2   = screenij2 * s2
            screenij3   = screenij3 * s3
            screenij54  = screenij54 * s54

            smax=max(screenijp,screenij0,screenij1,screenij2,screenij3,screenij54)
            if (smax .lt. scnres) exit

         enddo ! end loop on ll

         scrnabp(neighOnI) = screenijp
            screep = scrnabp(neighOnI)
         scrnab0(neighOnI) = screenij0
            scree0 = scrnab0(neighOnI)
         scrnab1(neighOnI) = screenij1
            scree1 = scrnab1(neighOnI)
         scrnab2(neighOnI) = screenij2
            scree2 = scrnab2(neighOnI)
         scrnab3(neighOnI)= screenij3
            scree3 = scrnab3(neighOnI)
         scrnab54(neighOnI) = screenij54
            scree54 = scrnab54(neighOnI)

         smin=min(screep,scree0,scree1,scree2,scree3,scree54)
         smax=max(screep,scree0,scree1,scree2,scree3,scree54)

         if (smin .gt. 1.0_cd-scnres .or. smax .lt. scnres) cycle

         ! Compute dscrnab for i->j
         !

         dscrpp(1:3) = 0.0_cd
         dscrp0(1:3) = 0.0_cd
         dscrp1(1:3) = 0.0_cd
         dscrp2(1:3) = 0.0_cd
         dscrp3(1:3) = 0.0_cd
         dscrp54(1:3) = 0.0_cd

         dscrmp(1:3) = 0.0_cd
         dscrm0(1:3) = 0.0_cd
         dscrm1(1:3) = 0.0_cd
         dscrm2(1:3) = 0.0_cd
         dscrm3(1:3) = 0.0_cd
         dscrm54(1:3) = 0.0_cd

         do m = 1,3

            Rpij = Rsqij - h2 * Rij(m) + hh
            call rscrn(Rpij,meatable,posRpij)
            dscrpp(m) = posRpij
            dscrp0(m) = posRpij
            dscrp1(m) = posRpij
            dscrp2(m) = posRpij
            dscrp3(m) = posRpij
            dscrp54(m) = posRpij

            Rpij = Rsqij + h2 * Rij(m) + hh
            call rscrn(Rpij,meatable,negRpij)
            dscrmp(m) = negRpij
            dscrm0(m) = negRpij
            dscrm1(m) = negRpij
            dscrm2(m) = negRpij
            dscrm3(m) = negRpij
            dscrm54(m) = negRpij

         enddo

         do ll = 1,numnei
  
            if (ll .eq. jj) cycle

            l = nei1atom(ll)

            Ril(:) = -(coor(:,l) - coor(:,i))
            Rsqil = dot_product(Ril,Ril)

            Rjl(:) = -(coor(:,l) - coor(:,j))
            Rsqjl = dot_product(Rjl,Rjl)

            do m = 1,3
 
               Rpij = Rsqij - h2 * Rij(m) + hh
               Rpil = Rsqil - h2 * Ril(m) + hh
               cx = get_c(Rpij,Rpil,Rsqjl)
               call get_s (cx,ity,jty,kty,meatable,dpp,dp0,dp1,dp2,dp3,dp54)

               dscrpp(m) = dscrpp(m)*dpp
               dscrp0(m) = dscrp0(m)*dp0
               dscrp1(m) = dscrp1(m)*dp1
               dscrp2(m) = dscrp2(m)*dp2
               dscrp3(m) = dscrp3(m)*dp3
               dscrp54(m) = dscrp54(m)*dp54

               Rpij = Rsqij + h2 * Rij(m) + hh
               Rpil = Rsqil + h2 * Ril(m) + hh
               cx = get_c(Rpij,Rpil,Rsqjl)
               call get_s (cx,ity,jty,kty,meatable,dmp,dm0,dm1,dm2,dm3,dm54)

               dscrmp(m) = dscrmp(m)*dmp
               dscrm0(m) = dscrm0(m)*dm0
               dscrm1(m) = dscrm1(m)*dm1
               dscrm2(m) = dscrm2(m)*dm2
               dscrm3(m) = dscrm3(m)*dm3
               dscrm54(m) = dscrm54(m)*dm54

            enddo

         enddo ! end loop on ll

         do m = 1,3

            dscrnabp(i2jl,m,neighOnI) = 0.5_cd*(dscrpp(m)-dscrmp(m))/h
            dscrnab0(i2jl,m,neighOnI) = 0.5_cd*(dscrp0(m)-dscrm0(m))/h
            dscrnab1(i2jl,m,neighOnI) = 0.5_cd*(dscrp1(m)-dscrm1(m))/h
            dscrnab2(i2jl,m,neighOnI) = 0.5_cd*(dscrp2(m)-dscrm2(m))/h
            dscrnab3(i2jl,m,neighOnI) = 0.5_cd*(dscrp3(m)-dscrm3(m))/h
            dscrnab54(i2jl,m,neighOnI)= 0.5_cd*(dscrp54(m)-dscrm54(m))/h

         enddo

         ! Compute dscrnab for j->i
         !

         dscrpp(1:3) = 0.0_cd
         dscrp0(1:3) = 0.0_cd
         dscrp1(1:3) = 0.0_cd
         dscrp2(1:3) = 0.0_cd
         dscrp3(1:3) = 0.0_cd
         dscrp54(1:3) = 0.0_cd

         dscrmp(1:3) = 0.0_cd
         dscrm0(1:3) = 0.0_cd
         dscrm1(1:3) = 0.0_cd
         dscrm2(1:3) = 0.0_cd
         dscrm3(1:3) = 0.0_cd
         dscrm54(1:3) = 0.0_cd

         do m = 1,3

            Rpji = Rsqji - h2 * Rji(m) + hh
            call rscrn(Rpji,meatable,posRpji)
            dscrpp(m) = posRpji
            dscrp0(m) = posRpji
            dscrp1(m) = posRpji
            dscrp2(m) = posRpji
            dscrp3(m) = posRpji
            dscrp54(m) = posRpji

            Rpji = Rsqji + h2 * Rji(m) + hh
            call rscrn(Rpji,meatable,negRpji)
            dscrmp(m) = negRpji
            dscrm0(m) = negRpji
            dscrm1(m) = negRpji
            dscrm2(m) = negRpji
            dscrm3(m) = negRpji
            dscrm54(m) = negRpji

         enddo

         do ll = 1,numnei

            if (ll .eq. jj) cycle

            l = nei1atom(ll)

            Ril(:) = -(coor(:,l) - coor(:,i))
            Rsqil = dot_product(Ril,Ril)

            Rjl(:) = -(coor(:,l) - coor(:,j))
            Rsqjl = dot_product(Rjl,Rjl)

            do m = 1,3

               Rpji = Rsqji - h2 * Rji(m) + hh
               Rpjl = Rsqjl - h2 * Rjl(m) + hh
               cx = get_c(Rpji,Rsqil,Rpjl)
               call get_s (cx,ity,jty,kty,meatable,dpp,dp0,dp1,dp2,dp3,dp54)

               dscrpp(m) = dscrpp(m)*dpp
               dscrp0(m) = dscrp0(m)*dp0
               dscrp1(m) = dscrp1(m)*dp1
               dscrp2(m) = dscrp2(m)*dp2
               dscrp3(m) = dscrp3(m)*dp3
               dscrp54(m) = dscrp54(m)*dp54

               Rpji = Rsqji + h2 * Rji(m) + hh
               Rpjl = Rsqjl + h2 * Rjl(m) + hh
               cx = get_c(Rpji,Rsqil,Rpjl)
               call get_s (cx,ity,jty,kty,meatable,dmp,dm0,dm1,dm2,dm3,dm54)

               dscrmp(m) = dscrmp(m)*dmp
               dscrm0(m) = dscrm0(m)*dm0
               dscrm1(m) = dscrm1(m)*dm1
               dscrm2(m) = dscrm2(m)*dm2
               dscrm3(m) = dscrm3(m)*dm3
               dscrm54(m) = dscrm54(m)*dm54

            enddo

         enddo ! end loop on ll

         do m = 1,3

            dscrnabp(j2il,m,neighOnI) = 0.5_cd*(dscrpp(m)-dscrmp(m))/h
            dscrnab0(j2il,m,neighOnI) = 0.5_cd*(dscrp0(m)-dscrm0(m))/h
            dscrnab1(j2il,m,neighOnI) = 0.5_cd*(dscrp1(m)-dscrm1(m))/h
            dscrnab2(j2il,m,neighOnI) = 0.5_cd*(dscrp2(m)-dscrm2(m))/h
            dscrnab3(j2il,m,neighOnI) = 0.5_cd*(dscrp3(m)-dscrm3(m))/h
            dscrnab54(j2il,m,neighOnI)= 0.5_cd*(dscrp54(m)-dscrm54(m))/h

         enddo

      enddo  ! end loop on jj

   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Now compute the pair potential and embedding function for each atom i
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (neighOnI .gt. 0) then

      do jj = 1,neighOnI

         Rij(:) = -Rij_xyzN(:,jj)
         Rsqij = dot_product(Rij,Rij)
         r_ij = sqrt(Rsqij)
         screep = scrnabp(jj)
         scree0 = scrnab0(jj)
         scree1 = scrnab1(jj)
         scree2 = scrnab2(jj)
         scree3 = scrnab3(jj)
         scree54 = scrnab54(jj)

         smin=min(screep,scree0,scree1,scree2,scree3,scree54)
         smax=max(screep,scree0,scree1,scree2,scree3,scree54)

         if (smax .le. scnres) cycle

         !
         ! compute phi for atom i
         !

         esub=meatable%esubs(ity,jty)
         phi=phi+0.5_cd*calc_phi(r_ij,ity,jty,meatable)*screep

         !
         ! compute background electron density for atom i
         !

         rej=meatable%res(ity,jty)

         do ll=0,4
            ispin=1
            rhjp(ll)=rhof(r_ij,ll,ispin,jty,meatable)
            ispin=2
            rhjm(ll)=rhof(r_ij,ll,ispin,jty,meatable)
         enddo

         r1 = 1.0_cd / r_ij
         rcos(:)= r1 * Rij(:)
         do m=1,3
            t1=rcos(m)
            ap(m)=ap(m)+rhjp(1)*t1*scree1
            am(m)=am(m)+rhjm(1)*t1*scree1
            agp(m)=agp(m)+rhjp(3)*t1*scree3
            agm(m)=agm(m)+rhjm(3)*t1*scree3
            ahp(m)=ahp(m)+rhjp(4)*t1*scree54
            ahm(m)=ahm(m)+rhjm(4)*t1*scree54
            do l=1,3
               t2=rcos(m)*rcos(l)
               bp(l,m)=bp(l,m)+rhjp(2)*t2*scree2
               bm(l,m)=bm(l,m)+rhjm(2)*t2*scree2
               do s=1,3
                  t3=rcos(m)*rcos(l)*rcos(s)
                  dp(s,l,m)= dp(s,l,m)+rhjp(3)*t3*scree3
                  dm(s,l,m)= dm(s,l,m)+rhjm(3)*t3*scree3
                  dhp(s,l,m)=dhp(s,l,m)+rhjp(4)*t3*scree54
                  dhm(s,l,m)=dhm(s,l,m)+rhjm(4)*t3*scree54
                  do k=1,3
                     do w=1,3
                        t54=rcos(m)*rcos(l)*rcos(s)*rcos(k)*rcos(w)
                        gp(w,k,s,l,m)=gp(w,k,s,l,m)+rhjp(4)*t54*scree54
                        gm(w,k,s,l,m)=gm(w,k,s,l,m)+rhjm(4)*t54*scree54
                     enddo
                  enddo
               enddo
            enddo
         enddo

         cp=cp+rhjp(0)*scree0
         cm=cm+rhjm(0)*scree0

         cgp=cgp+rhjp(2)*scree2
         cgm=cgm+rhjm(2)*scree2

      enddo  ! loop on jj

   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Now that we know the + and - partial electron densities, calculate embedding
!  energy and the deriviative of embedding energy, the derivative of the
!  pair potential, to obtain the particle energy, total energy, partical forces,
!  and viral forces.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (neighOnI .gt. 0) then

      rho0=cp**2-cm**2
      if(rho0.gt.0.0_cd) then

         rho0=sqrt(rho0)
         asub=meatable%asubs(ity)
         bsub=meatable%bsubs(ity)
         esub=meatable%esubs(ity,jty)

         rho1p=0.0_cd
         rho1m=0.0_cd
         rho2p=-(cgp**2)/3.0_cd
         rho2m=-(cgm**2)/3.0_cd
         rho3p=0.0_cd
         rho3m=0.0_cd
         rho54p=0.0_cd
         rho54m=0.0_cd
         do m=1,3
            rho3p=rho3p-legend00*agp(m)**2
            rho3m=rho3m-legend00*agm(m)**2
            rho54p=rho54p+legend02*ahp(m)**2
            rho54m=rho54m+legend02*ahm(m)**2
            do l=1,3
               do s=1,3
                  rho54p=rho54p-legend01*dhp(s,l,m)**2
                  rho54m=rho54m-legend01*dhm(s,l,m)**2
               enddo
            enddo
         enddo
         do m=1,3
            rho1p=rho1p+ap(m)**2
            rho1m=rho1m+am(m)**2
            do l=1,3
               rho2p=rho2p+bp(l,m)**2
               rho2m=rho2m+bm(l,m)**2
               do s=1,3
                  rho3p=rho3p+dp(s,l,m)**2
                  rho3m=rho3m+dm(s,l,m)**2
                  do k=1,3
                     do w=1,3
                        rho54p=rho54p+gp(w,k,s,l,m)**2
                        rho54m=rho54m+gm(w,k,s,l,m)**2
                     enddo
                  enddo
               enddo
            enddo
          enddo

!        Eqn. (9a)
!
          AA=rho1p+rho2p+rho3p+rho54p-(rho1m+rho2m+rho3m+rho54m)

          call elecden(rho0,AA,angular1,angular2,rhobar)
          dang1=angular1
          dang2=angular2
          plocal=frhoi(rhobar,ity,asub,bsub,esub,meatable)

          if (rhobar.lt.0.0_cd) then
              print *,'**********warning**********'
              print *, 'for atom ',i,' rhobar= ',rhobar,plocal
              goto 42
          endif

      endif

      if(comp_enepot.eq.1) then
         enepot(i) = enepot(i) + plocal + phi
!!         print *, 'rho1p,rho2p,rho3p,rho5p',rho1p,rho2p,rho3p,rho54p
!!         print *, 'rho1m,rho2m,rho3m,rho5m',rho1m,rho2m,rho3m,rho54m
!!         print *, 'rho0,rho1,rho2,rho3,rho5',rho0,rho1p-rho1m,rho2p-rho2m,rho3p-rho3m,rho54p-rho54m
!!         print *, 'i,enepot,plocal,phi,rhobar ', i,enepot(i),plocal,phi,rhobar
      endif

      if(comp_energy.eq.1) then
         energy = energy + plocal + phi
!!         print *,'atom = ',i,'; energy = ',energy
      endif

   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Also, get the particle forces and the viral forces (if requested)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (neighOnI .gt. 0) then

      rho0=cp**2-cm**2
      if(rho0.lt.0.0_cd) goto 42
      rho0=sqrt(rho0)
      asub=meatable%asubs(ity)
      bsub=meatable%bsubs(ity)
      esub=meatable%esubs(ity,jty)
      roz=rhobar
      dfrho=dfrhoi(roz,ity,asub,bsub,esub,meatable)

      do jj = 1,neighOnI

         j = jIdxOnI(jj)

         Rij(:) = -Rij_xyzN(:,jj)
         Rsqij = dot_product(Rij,Rij)
         r_ij = sqrt(Rsqij)
         phiij = calc_phi(r_ij,ity,jty,meatable)
         screep = scrnabp(jj)
         scree0 = scrnab0(jj)
         scree1 = scrnab1(jj)
         scree2 = scrnab2(jj)
         scree3 = scrnab3(jj)
         scree54 = scrnab54(jj)

         smin=min(screep,scree0,scree1,scree2,scree3,scree54)
         smax=max(screep,scree0,scree1,scree2,scree3,scree54)

         if (smax .le. scnres) cycle

         !
         !     Find interatomic distances, normalized directional vectors, etc.
         !

         rs=r_ij
         dxij=Rij(1)
         dyij=Rij(2)
         dzij=Rij(3)
         rr(:)=Rij(:)
         rcos(:)=rr(:)/rs

         call dcmij(drcos,rr,rs)

         phiif=dphif(r_ij,ity,jty,meatable)

         !
         !     Compile force on i and jindx atoms
         !
         !     note the minus sign on j term as the direc. cosine is opposite
         !

         do m=1,3

            dfij(m)=phiif*rcos(m)*scrnabp(jj)-phiij*dscrnabp(i2jl,m,jj)

            if(comp_force.eq.1) then
               force(m,i) = force(m,i) - dfij(m)
            endif

         enddo

         slocal(1,1,i) = slocal(1,1,i) - dfij(1)*dxij
         slocal(2,2,i) = slocal(2,2,i) - dfij(2)*dyij
         slocal(3,3,i) = slocal(3,3,i) - dfij(3)*dzij
         slocal(2,3,i) = slocal(2,3,i) - dfij(2)*dzij
         slocal(1,3,i) = slocal(1,3,i) - dfij(1)*dzij
         slocal(1,2,i) = slocal(1,2,i) - dfij(1)*dyij

         !
         !     Parameters for calculating embedding potentials
         !

         re=meatable%res(ity,jty)

         do ll=0,4
            ispin=1
            drhjp(ll)=drhof(rs,ll,ispin,jty,meatable)
            rhjp(ll)=rhof(rs,ll,ispin,jty,meatable)
            ispin=2
            drhjm(ll)=drhof(rs,ll,ispin,jty,meatable)
            rhjm(ll)=rhof(rs,ll,ispin,jty,meatable)
         enddo

         tmp0p = (2.0_cd / 3.0_cd) * cgp
         tmp0m = (2.0_cd / 3.0_cd) * cgm

         do kk=1,3
            rcosj = rcos(kk)

            !
            !     Initialize accumulation variables
            !

            drho0 = (cp*drhjp(0) - cm*drhjm(0))/rho0 * scree0 * rcosj
            drho0s = (cp*rhjp(0) - cm*rhjm(0))/rho0
            drho1 = 0.0_cd
            drho1g = 0.0_cd
            drho1gh = 0.0_cd
            drho1s = 0.0_cd
            drho1sg = 0.0_cd
            drho1sgh = 0.0_cd
            drho2 = -(tmp0p * drhjp(2) - tmp0m * drhjm(2)) * scree2 * rcosj
            drho2s = -rhjp(2) * tmp0p + rhjm(2) * tmp0m
            drho3 = 0.0_cd
            drho3h = 0.0_cd
            drho3s = 0.0_cd
            drho3sh = 0.0_cd
            drho5 = 0.0_cd
            drho5s = 0.0_cd

            do m = 1 , 3
               rcosm = rcos(m)
               tmp1p = 2.0_cd * ap(m)
               tmp1m = 2.0_cd * am(m)
               tmp1gp = 2.0_cd * agp(m)
               tmp1gm = 2.0_cd * agm(m)
               tmp1hp = 2.0_cd * ahp(m)
               tmp1hm = 2.0_cd * ahm(m)

               drho1 = drho1 + tmp1p * (drhjp(1) * rcosj * rcosm + rhjp(1) * drcos(m,kk)) * scree1 &
                             - tmp1m * (drhjm(1) * rcosj * rcosm + rhjm(1) * drcos(m,kk)) * scree1

               drho1g = drho1g + tmp1gp * (drhjp(3) * rcosj * rcosm + rhjp(3) * drcos(m,kk)) * scree3 &
                               - tmp1gm * (drhjm(3) * rcosj * rcosm + rhjm(3) * drcos(m,kk)) * scree3

               drho1gh = drho1gh + tmp1hp * (drhjp(4) * rcosj * rcosm + rhjp(4) * drcos(m,kk)) * scree54 &
                                 - tmp1hm * (drhjm(4) * rcosj * rcosm + rhjm(4) * drcos(m,kk)) * scree54

               drho1s = drho1s + tmp1p * rhjp(1) * rcosm &
                               - tmp1m * rhjm(1) * rcosm

               drho1sg = drho1sg + tmp1gp * rhjp(3) * rcosm &
                                 - tmp1gm * rhjm(3) * rcosm

               drho1sgh = drho1sgh + tmp1hp * rhjp(4) * rcosm &
                                   - tmp1hm * rhjm(4) * rcosm

               do l = 1, 3
                  rcosl = rcos(l)
                  rml = rcosm * rcosl
                  drmllm = drcos(m,kk) * rcosl + rcosm * drcos(l,kk)

                  tmp2p = 2.0_cd * bp(l,m)
                  tmp2m = 2.0_cd * bm(l,m)

                  drho2 = drho2 + tmp2p * (drhjp(2) * rcosj * rml + rhjp(2) * drmllm) * scree2 &
                                - tmp2m * (drhjm(2) * rcosj * rml + rhjm(2) * drmllm) * scree2

                  drho2s = drho2s + tmp2p * rhjp(2) * rml &
                                  - tmp2m * rhjm(2) * rml

                  do s = 1, 3
                     rcoss = rcos(s)
                     rmls = rml * rcoss
                     tmp3p = 2.0_cd * dp(s,l,m)
                     tmp3m = 2.0_cd * dm(s,l,m)
                     tmp3hp = 2.0_cd * dhp(s,l,m)
                     tmp3hm = 2.0_cd * dhm(s,l,m)
                     drmlsslm = rml * drcos(s,kk) + rcoss * drmllm

                     drho3 = drho3 + tmp3p * (drhjp(3) * rcosj * rmls + rhjp(3) * drmlsslm) * scree3 &
                                   - tmp3m * (drhjm(3) * rcosj * rmls + rhjm(3) * drmlsslm) * scree3

                     drho3h = drho3h + tmp3hp * (drhjp(4) * rcosj * rmls + rhjp(4) * drmlsslm) * scree54 &
                                     - tmp3hm * (drhjm(4) * rcosj * rmls + rhjm(4) * drmlsslm) * scree54

                     drho3s = drho3s + tmp3p * rhjp(3) * rmls &
                                     - tmp3m * rhjm(3) * rmls

                     drho3sh = drho3sh + tmp3hp * rhjp(4) * rmls &
                                       - tmp3hm * rhjm(4) * rmls

                     do k = 1 , 3
                        rcosk = rcos(k)
                        rmlsk = rmls * rcosk
                        drmlskkslm = drmlsslm * rcosk + rmls * drcos(k,kk)

                        do w = 1, 3
                           rcosw = rcos(w)
                           rmlskw = rmlsk * rcosw
                           drmlskwwkslm = drmlskkslm * rcosw + rmlsk * drcos(w,kk)

                           tmp5p = 2.0_cd * gp(w,k,s,l,m)
                           tmp5m = 2.0_cd * gm(w,k,s,l,m)

                           drho5 = drho5 + tmp5p * (drhjp(4) * rcosj * rmlskw + rhjp(4) * drmlskwwkslm) * scree54 &
                                         - tmp5m * (drhjm(4) * rcosj * rmlskw + rhjm(4) * drmlskwwkslm) * scree54

                           drho5s = drho5s + tmp5p * rhjp(4) * rmlskw &
                                           - tmp5m * rhjm(4) * rmlskw

                        enddo
                     enddo
                  enddo
               enddo
            enddo

            !
            !     Compile forces for i -> j and j -> i
            !

            term1 =  drho1 &
                  +  drho2 &
                  + (drho3 - legend00*drho1g) &
                  + (drho5 - legend01*drho3h + legend02*drho1gh)

            term2 =  drho1s * dscrnab1(i2jl,kk,jj) &
                  +  drho2s * dscrnab2(i2jl,kk,jj) &
                  + (drho3s - legend00*drho1sg) * dscrnab3(i2jl,kk,jj)  &
                  + (drho5s - legend01*drho3sh  + legend02*drho1sgh) * dscrnab54(i2jl,kk,jj)

            tmp4 = dscrnab0(i2jl,kk,jj)

            term3 =  drho1s * dscrnab1(j2il,kk,jj) &
                  +  drho2s * dscrnab2(j2il,kk,jj) &
                  + (drho3s - legend00*drho1sg) * dscrnab3(j2il,kk,jj)  &
                  + (drho5s - legend01*drho3sh  + legend02*drho1sgh) * dscrnab54(j2il,kk,jj)

            tmp8 = dscrnab0(j2il,kk,jj)

            dfij(kk) = dfrho * ( (drho0 - drho0s * tmp4) * dang1 + (term1 - term2) * dang2 )
            dfji(kk) = dfrho * ( (drho0 + drho0s * tmp8) * dang1 + (term1 + term3) * dang2 )

            if(comp_force.eq.1) then
               force(kk,i) = force(kk,i) - dfij(kk)
               force(kk,j) = force(kk,j) + dfji(kk)
            endif

         enddo

         slocal(1,1,i) = slocal(1,1,i) - dfij(1)*dxij
         slocal(2,2,i) = slocal(2,2,i) - dfij(2)*dyij
         slocal(3,3,i) = slocal(3,3,i) - dfij(3)*dzij
         slocal(2,3,i) = slocal(2,3,i) - dfij(2)*dzij
         slocal(1,3,i) = slocal(1,3,i) - dfij(1)*dzij
         slocal(1,2,i) = slocal(1,2,i) - dfij(1)*dyij

         slocal(1,1,j) = slocal(1,1,j) - dfji(1)*dxij
         slocal(2,2,j) = slocal(2,2,j) - dfji(2)*dyij
         slocal(3,3,j) = slocal(3,3,j) - dfji(3)*dzij
         slocal(2,3,j) = slocal(2,3,j) - dfji(2)*dzij
         slocal(1,3,j) = slocal(1,3,j) - dfji(1)*dzij
         slocal(1,2,j) = slocal(1,2,j) - dfji(1)*dyij

         do mm = 1,neighOnI

            if (mm .eq. jj) cycle

            mc = jIdxOnI(mm)

            Rim(:) = -Rij_xyzN(:,mm)
            Rsqim = dot_product(Rim,Rim)

            dxim = Rim(1)
            dyim = Rim(2)
            dzim = Rim(3)

            Rsqmj = Rsqij + Rsqim - 2.0_cd * (dxij*dxim + dyij*dyim + dzij*dzim)

            dxmj = 2.0 * dxim - dxij
            dymj = 2.0 * dyim - dyij
            dzmj = 2.0 * dzim - dzij

            cx =  get_c(Rsqij,Rsqim,Rsqmj)
            call get_s (cx,ity,jty,kty,meatable,sp,s0,s1,s2,s3,s54)

            do kk = 1,3

               Rpim = Rsqim + h2 * Rim(kk) + hh
               Rpmj = Rsqmj + h2 * (Rim(kk) - Rij(kk)) + hh
               cx = get_c(Rsqij,Rpim,Rpmj)
               call get_s (cx,ity,jty,kty,meatable,spp,sp0,sp1,sp2,sp3,sp54)

               Rpim = Rsqim - h2 * Rim(kk) + hh
               Rpmj = Rsqmj - h2 * (Rim(kk) - Rij(kk)) + hh
               cx = get_c(Rsqij,Rpim,Rpmj)
               call get_s (cx,ity,jty,kty,meatable,smp,sm0,sm1,sm2,sm3,sm54)

               if(sp.gt.scnres .and. sp.lt.1.0_cd) then
                  dscrnabkp=0.5_cd*((spp-smp)*screep)/(sp*h)
               else
                  dscrnabkp=0.0_cd
               endif

               if(s0.gt.scnres .and. s0.lt.1.0_cd) then
                  dscrnabk0=0.5_cd*((sp0-sm0)*scree0)/(s0*h)
               else
                  dscrnabk0=0.0_cd
               endif

               if(s1.gt.scnres .and. s1.lt.1.0_cd) then
                  dscrnabk1=0.5_cd*((sp1-sm1)*scree1)/(s1*h)
               else
                  dscrnabk1=0.0_cd
               endif

               if(s2.gt.scnres .and. s2.lt.1.0_cd) then
                  dscrnabk2=0.5_cd*((sp2-sm2)*scree2)/(s2*h)
               else
                  dscrnabk2=0.0_cd
               endif

               if(s3.gt.scnres .and. s3.lt.1.0_cd) then
                  dscrnabk3=0.5_cd*((sp3-sm3)*scree3)/(s3*h)
               else
                  dscrnabk3=0.0_cd
               endif

               if(s54.gt.scnres .and. s54.lt.1.0_cd) then
                  dscrnabk54=0.5_cd*((sp54-sm54)*scree54)/(s54*h)
               else
                  dscrnabk54=0.0_cd
               endif

               rcosj = rcos(kk)
               drho0s = (cp*rhjp(0)-cm*rhjm(0))/rho0
               drho1s = 0.0_cd
               drho1sg = 0.0_cd
               drho1sgh = 0.0_cd
               drho2s = -(rhjp(2) * tmp0p - rhjm(2) * tmp0m)
               drho3s = 0.0_cd
               drho3sh = 0.0_cd
               drho5s = 0.0_cd

               do m = 1,3
                  rcosm = rcos(m)
                  tmp1p = 2.0_cd * ap(m)
                  tmp1gp = 2.0_cd * agp(m)
                  tmp1hp = 2.0_cd * ahp(m)
                  tmp1m = 2.0_cd * am(m)
                  tmp1gm = 2.0_cd * agm(m)
                  tmp1hm = 2.0_cd * ahm(m)
                  drho1s = drho1s + tmp1p * rhjp(1) * rcosm &
                                  - tmp1m * rhjm(1) * rcosm
                  drho1sg = drho1sg + tmp1gp * rhjp(3) * rcosm &
                                    - tmp1gm * rhjm(3) * rcosm
                  drho1sgh = drho1sgh + tmp1hp * rhjp(4) * rcosm &
                                         - tmp1hm * rhjm(4) * rcosm

                  do l = 1,3
                     rcosl = rcos(l)
                     rml = rcosm * rcosl
                     tmp2p = 2.0_cd * bp(l,m)
                     tmp2m = 2.0_cd * bm(l,m)
                     drho2s = drho2s + tmp2p * rhjp(2) * rml &
                                     - tmp2m * rhjm(2) * rml

                     do s = 1,3
                        rcoss = rcos(s)
                        rmls = rml * rcoss
                        tmp3p = 2.0_cd * dp(s,l,m)
                        tmp3m = 2.0_cd * dm(s,l,m)
                        tmp3hp = 2.0_cd * dhp(s,l,m)
                        tmp3hm = 2.0_cd * dhm(s,l,m)
                        drho3s = drho3s + tmp3p * rhjp(3) * rmls &
                                        - tmp3m * rhjm(3) * rmls
                        drho3sh = drho3sh + tmp3hp * rhjp(4) * rmls &
                                          - tmp3hm * rhjm(4) * rmls

                        do k = 1,3
                           rcosk = rcos(k)
                           rmlsk = rmls * rcosk

                           do w = 1,3
                              rcosw = rcos(w)
                              rmlskw = rmlsk * rcosw
                              tmp5p = 2.0_cd * gp(w,k,s,l,m)
                              tmp5m = 2.0_cd * gm(w,k,s,l,m)
                              drho5s = drho5s + tmp5p * rhjp(4) * rmlskw &
                                              - tmp5m * rhjm(4) * rmlskw

                           enddo

                        enddo

                     enddo

                  enddo

               enddo

               term2 =  drho1s * dscrnabk1 &
                     +  drho2s * dscrnabk2 &
                     + (drho3s - legend00 * drho1sg) * dscrnabk3 &
                     + (drho5s - legend01 * drho3sh + legend02 * drho1sgh) * dscrnabk54

               dfim(kk) = 0.5 * phiij * dscrnabkp + dfrho * (dang1 * drho0s * dscrnabk0 + dang2 * term2)

               if(comp_force.eq.1) then
                  force(kk,mc) = force(kk,mc) + dfim(kk)
               endif

            enddo

            slocal(1,1,mc) = slocal(1,1,mc) - dfim(1)*dxmj
            slocal(2,2,mc) = slocal(2,2,mc) - dfim(2)*dymj
            slocal(3,3,mc) = slocal(3,3,mc) - dfim(3)*dzmj
            slocal(2,3,mc) = slocal(2,3,mc) - dfim(2)*dzmj
            slocal(1,3,mc) = slocal(1,3,mc) - dfim(1)*dzmj
            slocal(1,2,mc) = slocal(1,2,mc) - dfim(1)*dymj

         enddo

      enddo

   endif

enddo

if(comp_ptvirl.eq.1) then
   do atom = 1,N 
      particleVirial(1,atom) = particleVirial(1,atom) - 0.5_cd*slocal(1,1,atom)
      particleVirial(2,atom) = particleVirial(2,atom) - 0.5_cd*slocal(2,2,atom)
      particleVirial(3,atom) = particleVirial(3,atom) - 0.5_cd*slocal(3,3,atom)
      particleVirial(4,atom) = particleVirial(4,atom) - 0.5_cd*slocal(2,3,atom)
      particleVirial(5,atom) = particleVirial(5,atom) - 0.5_cd*slocal(1,3,atom)
      particleVirial(6,atom) = particleVirial(6,atom) - 0.5_cd*slocal(1,2,atom)
   enddo
endif

if(comp_virial.eq.1) then
   do atom = 1,N
      virial(1) = virial(1) - 0.5_cd*slocal(1,1,atom)
      virial(2) = virial(2) - 0.5_cd*slocal(2,2,atom)
      virial(3) = virial(3) - 0.5_cd*slocal(3,3,atom)
      virial(4) = virial(4) - 0.5_cd*slocal(2,3,atom)
      virial(5) = virial(5) - 0.5_cd*slocal(1,3,atom)
      virial(6) = virial(6) - 0.5_cd*slocal(1,2,atom)
   enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

deallocate(slocal)

! Everything is great
!
ier = KIM_STATUS_OK
42 continue
Compute_Energy_Forces = ier
return

end function Compute_Energy_Forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Model driver destroy routine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer(c_int) function destroy(pkim) bind(c)
use KIM_API_F03
implicit none

!-- Transferred variables
type(c_ptr), intent(in) :: pkim

!-- Local variables
integer(c_int) :: ier, idum
type(model_parameters), pointer :: meatable; type(c_ptr) :: pmeatable

pmeatable = kim_api_get_model_buffer(pkim, ier)
if (ier.ne.KIM_STATUS_OK) then
  idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
    "kim_api_get_model_buffer", ier)
  destroy = ier
  return
endif
call c_f_pointer(pmeatable, meatable)

deallocate(meatable%asubs)
deallocate(meatable%bsubs)
deallocate(meatable%esubs)
deallocate(meatable%res)
deallocate(meatable%embedtab)
deallocate(meatable%phitab)
deallocate(meatable%rhotab)
deallocate(meatable%scrp)
deallocate(meatable%scr0)
deallocate(meatable%scr1)
deallocate(meatable%scr2)
deallocate(meatable%scr3)
deallocate(meatable%scr54)

deallocate( meatable )

destroy = KIM_STATUS_OK
return

end function destroy

end module MSMEAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Model driver initialization routine (REQUIRED)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer(c_int) function model_driver_init(pkim, pparamfile, nmstrlen, &
                                          numparamfiles) bind(c)
use, intrinsic :: iso_c_binding
use MSMEAM
use KIM_API_F03
implicit none
integer(c_int), parameter :: cd = c_double ! used for literal constants

!-- Transferred variables
type(c_ptr),        intent(in) :: pkim
type(c_ptr), value, intent(in) :: pparamfile
integer(c_int),     intent(in) :: nmstrlen
integer(c_int),     intent(in) :: numparamfiles
character(len=nmstrlen), pointer :: paramfile(:)
character (len=80) :: error_message

!-- Local variables
integer(c_int), parameter :: zero=0
integer(c_int), parameter :: one=1
integer(c_int), parameter :: two=2
integer(c_int), parameter :: three=3
integer(c_int), parameter :: four=4
integer(c_int), parameter :: five=5
integer(c_int)            :: i,j,ier,idum,itab,nelem,mil1
integer(c_int)            :: ispin,it,jt,kt,return_error
type(model_parameters), pointer :: meatable

integer(c_int), parameter :: nelmax=2
real(c_double)            :: rx,embedtemp,phitemp

!-- KIM variables
real(c_double), pointer :: cutoff; type(c_ptr) :: pcutoff
character(len=KIM_KEY_STRING_LENGTH) :: NBC_Method

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
allocate( meatable )
call kim_api_set_model_buffer(pkim, c_loc(meatable), ier)
if (ier.ne.KIM_STATUS_OK) then
  idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
    "kim_api_set_model_buffer", ier)
  return_error = ier
  goto 1000
endif

! Read in model parameters from parameter file
!
open(10,file=paramfile(1),status="old")
read(10,*,iostat=ier,err=100)
read(10,*,iostat=ier,err=100)
read(10,*,iostat=ier,err=100) &
      &  meatable%ntab,meatable%nel,meatable%rtabmax,meatable%rhotabmin,meatable%rhotabmax
read(10,*,iostat=ier,err=100) meatable%model_cutoff

nelem=meatable%nel
mil1=meatable%ntab

allocate(meatable%asubs(one:nelem))
allocate(meatable%bsubs(one:nelem))
allocate(meatable%esubs(one:nelem,one:nelem))
allocate(meatable%res(one:nelem,one:nelem))
allocate(meatable%embedtab(one:mil1,one:nelem))
allocate(meatable%phitab(one:mil1,one:nelem,one:nelem))
allocate(meatable%rhotab(one:mil1,zero:four,one:two,one:nelem))
allocate(meatable%scrp(one:mil1,one:nelem,one:nelem,one:nelem))
allocate(meatable%scr0(one:mil1,one:nelem,one:nelem,one:nelem))
allocate(meatable%scr1(one:mil1,one:nelem,one:nelem,one:nelem))
allocate(meatable%scr2(one:mil1,one:nelem,one:nelem,one:nelem))
allocate(meatable%scr3(one:mil1,one:nelem,one:nelem,one:nelem))
allocate(meatable%scr54(one:mil1,one:nelem,one:nelem,one:nelem))

do i=1,meatable%nel
    read(10,*,iostat=ier,err=100) &
      &  meatable%asubs(i),meatable%bsubs(i),meatable%res(i,i),meatable%esubs(i,i)
    meatable%rrtabmax=meatable%rtabmax*meatable%res(i,i)
         ispin=1
              do itab=1,mil1
                   read(10,*,iostat=ier,err=100) &
                   meatable%rhotab(itab,zero,ispin,i), &
                   meatable%rhotab(itab,one,ispin,i), &
                   meatable%rhotab(itab,two,ispin,i), &
                   meatable%rhotab(itab,three,ispin,i),&
                   meatable%rhotab(itab,four,ispin,i), &
                   embedtemp
                   meatable%embedtab(itab,i)=embedtemp*meatable%esubs(i,i)

                   !meatable%rhotab(itab,one,ispin,i)=0.0_cd
                   !meatable%rhotab(itab,two,ispin,i)=0.0_cd
                   !meatable%rhotab(itab,three,ispin,i)=0.0_cd
                   !meatable%rhotab(itab,four,ispin,i)=0.0_cd

              enddo
         ispin=2
              do itab=1,mil1
                   read(10,*,iostat=ier,err=100) &
                   meatable%rhotab(itab,zero,ispin,i), &
                   meatable%rhotab(itab,one,ispin,i), &
                   meatable%rhotab(itab,two,ispin,i), &
                   meatable%rhotab(itab,three,ispin,i), &
                   meatable%rhotab(itab,four,ispin,i)

                   !meatable%rhotab(itab,one,ispin,i)=0.0_cd
                   !meatable%rhotab(itab,two,ispin,i)=0.0_cd
                   !meatable%rhotab(itab,three,ispin,i)=0.0_cd
                   !meatable%rhotab(itab,four,ispin,i)=0.0_cd

              enddo
              do itab=1,mil1
                   read(10,*,iostat=ier,err=100) phitemp
                   meatable%phitab(itab,i,i)=phitemp*meatable%esubs(i,i)
                   rx=(dble(itab)/(dble(meatable%ntab)-1.0d0)*&
                     &(meatable%rrtabmax/meatable%res(i,i)))

              enddo
              it=i
              jt=i
              kt=i
              do itab=1,mil1
                  read(10,*,iostat=ier,err=100) &
                  meatable%scrp(itab,it,jt,kt), &
                  meatable%scr0(itab,it,jt,kt), &
                  meatable%scr1(itab,it,jt,kt), &
                  meatable%scr2(itab,it,jt,kt), &
                  meatable%scr3(itab,it,jt,kt), &
                  meatable%scr54(itab,it,jt,kt) 

                  !meatable%scrp(itab,it,jt,kt)=1.0_cd
                  !meatable%scr0(itab,it,jt,kt)=1.0_cd
                  !meatable%scr1(itab,it,jt,kt)=1.0_cd
                  !meatable%scr2(itab,it,jt,kt)=1.0_cd
                  !meatable%scr3(itab,it,jt,kt)=1.0_cd
                  !meatable%scr54(itab,it,jt,kt)=1.0_cd

              enddo
enddo
close(10)

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
   goto 1000
endif
if (index(NBC_Method,"NEIGH_PURE_F").eq.1) then
    meatable%NBC = 2
    meatable%HalfOrFull = 2
else
    ier = KIM_STATUS_FAIL
    idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                              "Unknown NBC method.", ier)
    goto 1000
endif


! Determine neighbor list handling mode
!
   !*****************************
   !* IterOrLoca = 1 -- Iterator
   !*            = 2 -- Locator
   !*****************************
   meatable%IterOrLoca = kim_api_get_neigh_mode(pkim, ier)
   if (ier.lt.KIM_STATUS_OK) then
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                    "kim_api_get_neigh_mode", ier)
      goto 1000
   endif
   if (meatable%IterOrLoca.ne.1 .and. meatable%IterOrLoca.ne.2) then
      ier = KIM_STATUS_FAIL
      write(error_message,'(a,i1)') &
         'Unsupported IterOrLoca mode = ',meatable%IterOrLoca
      idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                                  error_message, ier)
      goto 1000
   endif


goto 200
100 continue
! reading parameters failed
idum = kim_api_report_error(__LINE__, THIS_FILE_NAME, &
                            "Unable to read potential parameters, kimerror = ",KIM_STATUS_FAIL)
return_error = KIM_STATUS_FAIL
goto 1000

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
cutoff = meatable%model_cutoff
meatable%model_cutsq   = meatable%model_cutoff**2

1000 continue
model_driver_init = return_error
return

end function model_driver_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

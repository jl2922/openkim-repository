#include "KIM_API_status.h"
#define THIS_FILE_NAME __FILE__

      subroutine MEAM_potential(pkim,mp,N,coors,species,boxsize_in,
     &                          energy,ene_pot_kim,forces,virial,
     &                          comp_energy,comp_enepot,comp_force,
     &                          comp_virial,ier)

      use, intrinsic :: iso_c_binding
      use KIM_API_F03
      use meam_parameters
      implicit none

      type(c_ptr), intent(in) :: pkim
      type(meam_params) :: mp
      integer(c_int), intent(in)  :: N
      real(c_double), intent(in)  :: coors(3, N)
      integer(c_int), intent(in)  :: species(N)
      real(c_double), intent(in)  :: boxsize_in(3)
      real(c_double), intent(out) :: energy
      real(c_double), intent(out) :: ene_pot_kim(N)
      real(c_double), intent(out) :: forces(3, N)
      real(c_double), intent(out) :: virial(6)
      integer(c_int), intent(in)  :: comp_energy
      integer(c_int), intent(in)  :: comp_enepot
      integer(c_int), intent(in)  :: comp_force
      integer(c_int), intent(in)  :: comp_virial
      integer(c_int), intent(out) :: ier

c
c  Define Pairwise Potential and Make Tables
c
      common /abcI/ Natom,Idim,Iperform,Ishear
      common /atom/ IAtomType(200000)
      common /lat1/ pos(3,200000),vel(3,200000),acc(3,200000)
      common /lat2/ eps(3,3),eps_inv(3,3),BoxSize(3)
      common /nei4/ ListMC(1000)
      common /ener/ ene_pot(200000),ene_kin(200000)
      integer(c_int) Natom, Idim, Iperform, Ishear
      integer(c_int) IAtomType
      real(c_double) pos, vel, acc
      real(c_double) eps, eps_inv, BoxSize
      integer(c_int) ListMC
      real(c_double) ene_pot, ene_kin

      integer(c_int) Ipbc(3)
      real(c_double) f_stress(3,3)
      integer(c_int) i, NListMC
      real(c_double) ene_sum

      ! Initialize variables
      !
      eps(:,:) = 0.0_cd
      eps(1,1) = 1.0_cd
      eps(2,2) = 1.0_cd
      eps(3,3) = 1.0_cd
      BoxSize(:)=boxsize_in(:)
      Natom = N
      do i = 1,N
         pos(1,i) = coors(1,i)/BoxSize(1)
         pos(2,i) = coors(2,i)/BoxSize(2)
         pos(3,i) = coors(3,i)/BoxSize(3)
         iatomtype(i)=species(i)
      end do
      Ipbc(1)=1
      Ipbc(2)=1
      Ipbc(3)=1

      ! Compute energy and forces
      !
      ier = KIM_STATUS_OK
      NListMC=0

      if (comp_energy.eq.1) then
         call Compute_EAM_Energy(Ipbc,ene_sum,NListMC,pkim,mp,ier)
         if (ier.lt.KIM_STATUS_OK) return
      endif

      if (comp_force.eq.1) then
         call Compute_EAM_Forces(Ipbc,f_stress,pkim,mp,ier)
         if (ier.lt.KIM_STATUS_OK) return
      endif

      ! Place results in KIM variables
      !
      if (comp_energy.eq.1) energy = ene_sum
      if (comp_enepot.eq.1) ene_pot_kim(1:N) = ene_pot(1:N)
      if (comp_force.eq.1) then
         do i=1,N
            forces(:,i) = acc(:,i)*BoxSize(:)*mp%Amass(species(i))
         enddo
      endif
      if (comp_virial.eq.1) then
         virial(1) = f_stress(1,1)
         virial(2) = f_stress(2,2)
         virial(3) = f_stress(3,3)
         virial(4) = f_stress(2,3)
         virial(5) = f_stress(1,3)
         virial(6) = f_stress(1,2)
      endif

      return

      end subroutine MEAM_potential

c
c Read in MEAM parameters and generate necessary tables
c
      subroutine Define_EAM_Potential(mp,paramfile_name,nmstrlen,
     &           Ncomp,Icomp,Ierr)

      use, intrinsic :: iso_c_binding
      use KIM_API_F03
      use meam_parameters
c
c  Define Pairwise Potential and Make Tables
c
      implicit none

      type(meam_params) mp
      integer(c_int) nmstrlen
      character(len=nmstrlen) paramfile_name
      integer(c_int) Ncomp
      character(len=2) Icomp(Ncomp)
      integer(c_int) Ierr

      real(c_double) data_ele(14)
      character(len=2) element
      character(len=80) card,card2
      character(len=6) lattice,lattice_element(mp%max_comp)
      character(len=5) alloy,alloy_1_name,alloy_2_name
      character(len=8) alloy3,ternary_name
      real(c_double) DeltaRsq
      integer(c_int) i, j, k, n
      real(c_double) a2N, alat, alpha_ij, alphaN, astar, astar3, B_ij,
     &   C_iji, C_jij, C_iij, C_ijj, C_ijk, C_jik, C_ikj,
     &   Cx_iji, Cx_jij, Cx_iij, Cx_ijj, Cx_ikj, Cx_jik, Cx_ijk,
     &   d_ij, d_input, dEdR, delta_Ec, DeltaR, dFdR_i, dFdR_j, Dij,
     &   Phi, Phi_i, Phi_j, DPhi, DPhi_i, DPhi_j, dRho_i, dRho_j, Eu,
     &   F_i, F_j, GAMMA, R11factor, R12factor, R21factor, R22factor,
     &   R2factor, Re_ij, Re_input, Rfactor, Rho0, Rho1, Rho1P, Rho2,
     &   Rho2_1, Rho2_2, Rho2_3, Rho2P, Rho3, Rho3P, Rho_Bar_i,
     &   Rho_Bar_j, RhoFactor, rr, rr0_i, rr0_j, scr, scr2, scr_mix,
     &   scr_pure, t_1i, t_1j, t_2i, t_2j, term2, tmp_ij, Z2_i, Z2_j,
     &   Z_ij, ZR2, Rcutoff
      real(c_double), external :: screen

c
c Make sure there aren't too many components
c
      if (Ncomp .gt. mp%max_comp) then
         Ierr = KIM_STATUS_FAIL
         write(*,'(/,a,/)') ' **** ERROR : Too many components.'
         return
      endif

c
c Read in cutoff radius
c
      call Read_Cutoff(paramfile_name,nmstrlen,Ncomp,Icomp,Rcutoff,Ierr)
      if (Ierr.lt.KIM_STATUS_OK) return

c
c  ===========================================================================
c     Initialize
c  ===========================================================================
c
      mp%Rmin = 0.5_cd
      mp%Rsqmin = mp%Rmin*mp%Rmin
      mp%Rcutoff = Rcutoff
      DeltaRsq = ( mp%Rcutoff**2 - mp%RsqMin ) / ( mp%LineTable - 1 )
      mp%OverDeltaRsq = 1.0_cd / DeltaRsq
      mp%Dcutoff = 0.1_cd

c
c  ===========================================================================
c     Search for the elements and Read MEAM parameters
c  ===========================================================================
c
      open(unit=2,file=paramfile_name,status='old',err=4448)
      do n = 1, Ncomp
         read(2,'(a)') card
         read(2,'(a)') card2
    1    continue
         read(2,'(a)',end=4444) card
         read(2,'(a)') card2
         read(card,'(a2)') element
         if(element.ne.Icomp(n)) goto 1
         if(element.eq.'XX') goto 4444
         read(card,2000) lattice,mp%amass(n),(data_ele(k), k=1,6)
 2000    format(3x,a6,d12.3,d9.2,d10.3,d11.4,3d9.2)
         read(card2,2001) mp%Eu_d(n),(data_ele(k), k=7,14)
 2001    format(2x,f5.3,5x,2d9.2,3d10.2,d9.2,1x,2f4.2)
c
         mp%Ec(n,n) = Data_Ele(1)
         mp%Re(n)   = Data_Ele(2)
         mp%B(n,n)  = Data_Ele(3)
         mp%A(n)    = Data_Ele(4)
         mp%beta(n,4) = Data_Ele(5)
         mp%beta(n,1) = Data_Ele(6)
         mp%beta(n,2) = Data_Ele(7)
         mp%beta(n,3) = Data_Ele(8)
         mp%t(n,1) = Data_Ele(9)
         mp%t(n,2) = Data_Ele(10)
         mp%t(n,3) = Data_Ele(11)
         mp%RhoZ(n)   = Data_Ele(12)
         mp%Cmin(n,n,n) = Data_Ele(13)
         mp%Cmax(n,n,n) = Data_Ele(14)
         lattice_element(n) = lattice

         rewind(2)
c
c
c ---------------This Part should be further extended -------------------------
c
         if(lattice.eq.'FCC_A1') then
c
            mp%a2(n)   = dsqrt(2.0_cd)
            mp%Z(n)    = 12.0_cd
            scr = screen(0.5_cd,0.5_cd,mp%Cmax(n,n,n),mp%Cmin(n,n,n))
            mp%scr_2nn = scr * scr * scr * scr
            mp%Z2(n) = 6.0_cd * mp%scr_2nn
            ZR2 = mp%Z2(n) / mp%Z(n)
            alat = dsqrt(2.0_cd) * mp%Re(n)
            mp%omega(n,n) = alat * alat * alat / 4.0_cd
c
         elseif(lattice.eq.'BCC_A2') then
c
            mp%a2(n)   = dsqrt(4.0_cd/3.0_cd)
            mp%Z(n)    = 8.0_cd
            scr = screen(0.75_cd,0.75_cd,mp%Cmax(n,n,n),mp%Cmin(n,n,n))
            mp%scr_2nn = scr * scr * scr * scr
            mp%Z2(n) = 6.0_cd * mp%scr_2nn
            ZR2 = mp%Z2(n) / mp%Z(n)
            alat = 2.0_cd * mp%Re(n) / dsqrt(3.0_cd)
            mp%omega(n,n) = alat * alat * alat / 2.0_cd
c
         elseif(lattice.eq.'HCP_A3') then
c
            mp%a2(n)   = dsqrt(2.0_cd)
            mp%Z(n)    = 12.0_cd
            scr = screen(0.5_cd,0.5_cd,mp%Cmax(n,n,n),mp%Cmin(n,n,n))
            mp%scr_2nn = scr * scr * scr * scr
            mp%Z2(n) = 6.0_cd * mp%scr_2nn
            ZR2 = mp%Z2(n) / mp%Z(n)
            alat = mp%Re(n)
            mp%omega(n,n) = alat * alat * alat / mp%a2(n)
c
         elseif(lattice.eq.'DIA_A4') then
c
            mp%a2(n)   = 4.0_cd / dsqrt(6.0_cd)
            mp%Z(n)    = 4.0_cd
            mp%scr_2nn = screen(0.375_cd,0.375_cd,mp%Cmax(n,n,n),
     &                          mp%Cmin(n,n,n))
            mp%Z2(n) = 12.0_cd * mp%scr_2nn
            ZR2 = mp%Z2(n) / mp%Z(n)
            alat = 4.0_cd * mp%Re(n) / dsqrt(3.0_cd)
            mp%omega(n,n) = alat * alat * alat / 8.0_cd
c
         elseif(lattice.eq.'DIMER') then
c  This is an arbitrary treatment, especially those for lat and omega(n,n)
            mp%a2(n)   = dsqrt(2.0_cd)
            mp%Z(n)    = 1.0_cd
            mp%scr_2nn = 0.0_cd
            mp%Z2(n) = mp%scr_2nn
            ZR2 = mp%Z2(n) / mp%Z(n)
            alat = mp%Re(n)
            mp%omega(n,n) = alat * alat * alat
c
         else
            goto 4445
         endif
c
         mp%alpha(n,n) = dsqrt(9.0_cd*mp%B(n,n)*mp%omega(n,n)/
     &                   mp%Ec(n,n))
         if(lattice.eq.'DIMER') mp%alpha(n,n) = mp%B(n,n)
         alphaN = mp%alpha(n,n)
         a2N = mp%a2(n)
         R2factor = mp%a2(n) - 1.0_cd
         mp%Rho0_Bar(n) = mp%RhoZ(n) * mp%Z(n)
     &       + mp%Z2(n) * mp%RhoZ(n) * dexp ( - mp%beta(n,4) * R2factor)
c
c  add t(3) term in the case of HCP_A3 or DIA_A4 lattice
c
         if(lattice.eq.'HCP_A3' .or. lattice.eq.'DIA_A4') then
            if(lattice.eq.'HCP_A3') then
               Rho3 = 0.5_cd * mp%RhoZ(n)
     &    - mp%scr_2nn * dsqrt(2.0_cd) * mp%RhoZ(n) * dexp(-mp%beta(n,3)
     &                                            * R2factor)
               Rho3P = 4.0_cd * Rho3 * Rho3 / 3.0_cd
            elseif(lattice.eq.'DIA_A4') then
               Rho3 = mp%RhoZ(n)
               Rho3P = 32.0_cd * Rho3 * Rho3 / 9.0_cd
            endif
c
            Rho0 = mp%Rho0_Bar(n)
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.0_cd) RhoFactor = 1.e-20_cd
            GAMMA = mp%t(n,3) * Rho3P / RhoFactor
            mp%Rho0_Bar(n) = Rho0 * 2.0_cd / ( 1.0_cd + dexp(-GAMMA) )
c           Rho0_Bar(n) = Rho0 * dsqrt( 1.0_cd + GAMMA )
         endif
c
c  add t(1), t(2) and t(3) term in the case of DIMER elements (Gases)
c
         if(lattice.eq.'DIMER') then
            Rho1 = mp%RhoZ(n)
            Rho1P = Rho1 * Rho1
            Rho2 = mp%RhoZ(n)
            Rho2P = 2.0_cd / 3.0_cd * Rho2 * Rho2
            Rho3 = mp%RhoZ(n)
            Rho3P = 2.0_cd / 5.0_cd * Rho3 * Rho3
c
            Rho0 = mp%Rho0_Bar(n)
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.0_cd) RhoFactor = 1.e-20_cd
            GAMMA = (mp%t(n,1)*Rho1P+mp%t(n,2)*Rho2P+mp%t(n,3)*Rho3P)
     &            / RhoFactor
            mp%Rho0_Bar(n) = Rho0 * 2.0_cd / ( 1.0_cd + dexp(-GAMMA) )
c           Rho0_Bar(n) = Rho0 * dsqrt( 1.0_cd + GAMMA )
         endif
c
c -----------------------------------------------------------------------------
c
c     Make a table for pair potential and its derivative between each pair
c
         DeltaR = ( mp%Rcutoff - mp%RMin ) / dble( mp%LineTable - 1 )
         mp%OverDeltaR = 1.0_cd / DeltaR
c
         do k=1,mp%LineTable
            Dij = mp%RMin + (k-1) * DeltaR
            call Define_Pair_Potential
     &          (mp,N,lattice,alphaN,ZR2,a2N,Dij,Phi,DPhi)
            mp%PhiTab(k,n,n)  = Phi
            mp%DPhiTab(k,n,n) = DPhi
         enddo
      enddo
c
c  ===========================================================================
c     Search for the alloy parameters
c  ===========================================================================
c
      if(Ncomp.gt.1) then
      do i = 1, Ncomp-1
      do j = i+1, Ncomp
         alloy_1_name(1:2) = Icomp(i)
         alloy_1_name(3:3) = '-'
         alloy_1_name(4:5) = Icomp(j)
         alloy_2_name(1:2) = Icomp(j)
         alloy_2_name(3:3) = '-'
         alloy_2_name(4:5) = Icomp(i)
c
         read(2,'(a)') card
         read(2,'(a)') card2
    2    continue
         read(2,'(a)',end=4446) card
         read(card,'(a5)') alloy
         if(alloy.ne.alloy_1_name .and. alloy.ne.alloy_2_name) goto 2
         read(card,2010) lattice,delta_Ec,Re_input,B_ij,d_input
     &                 , C_iji, C_jij, C_iij, C_ijj
         read(2,2012) Cx_iji, Cx_jij, Cx_iij, Cx_ijj
c2010    format(7x,a6,5d12.4)
 2010    format(7x,a6,3d12.4,1X,5f5.2)
 2012    format(55x,4f5.2)
c
         rewind(2)
c
         if(C_iji.eq.0.0) then
            mp%Cmin(i,j,i) = mp%Cmin(i,i,i)
         else
            mp%Cmin(i,j,i) = C_iji
         endif
         if(C_jij.eq.0.0) then
            mp%Cmin(j,i,j) = mp%Cmin(j,j,j)
         else
            mp%Cmin(j,i,j) = C_jij
         endif
         if(Cx_iji.ne.0.0) mp%Cmax(i,j,i) = Cx_iji
         if(Cx_jij.ne.0.0) mp%Cmax(j,i,j) = Cx_jij
         if(Cx_iij.ne.0.0) then
            mp%Cmax(i,i,j) = Cx_iij
            mp%Cmax(j,i,i) = Cx_iij
         endif
         if(Cx_ijj.ne.0.0) then
            mp%Cmax(i,j,j) = Cx_ijj
            mp%Cmax(j,j,i) = Cx_ijj
         endif
c
c ---------------This Part should be further extended -------------------------
c
         if(lattice.eq.'FCC_B1' .or. lattice.eq.'NaCl') then
c
c          for fcc_B1/NaCl
c
            a2N      = dsqrt(2.0_cd)
            mp%Ec(i,j) = 0.5_cd * (mp%Ec(i,i) + mp%Ec(j,j)) - delta_Ec
            if(Re_input.eq.0.0_cd) then
               mp%omega(i,j) = 0.5_cd * (mp%omega(i,i) + mp%omega(j,j))
                   tmp_ij = 8.0_cd * mp%omega(i,j)
               Re_ij      = tmp_ij**(1.0_cd/3.0_cd) / 2.0_cd
            else
               Re_ij      =  Re_input
                   tmp_ij =  2.0_cd * Re_ij
               mp%omega(i,j) = 0.125_cd * tmp_ij * tmp_ij * tmp_ij
            endif
            if(B_ij.eq.0.0_cd) then
               mp%B(i,j) = 0.5_cd * mp%B(i,i) + 0.5_cd * mp%B(j,j)
            else
               mp%B(i,j) =  B_ij
            endif
            alpha_ij  = dsqrt(9.0_cd*mp%B(i,j)*mp%omega(i,j)/mp%Ec(i,j))

            if(d_input.eq.0.0_cd) then
               d_ij   = 0.5_cd * mp%Eu_d(i) + 0.5_cd * mp%Eu_d(j)
            elseif(d_input.lt.0.0_cd) then
               d_ij   =  0.0_cd
            else
               d_ij   =  d_input
            endif

            Z_ij      = 6.0_cd
              scr = screen(0.5_cd,0.5_cd,mp%Cmax(i,j,i),mp%Cmin(i,j,i))
              mp%scr_2nn = scr * scr
            Z2_i      = 12.0_cd * mp%scr_2nn
              scr = screen(0.5_cd,0.5_cd,mp%Cmax(j,i,j),mp%Cmin(j,i,j))
              mp%scr_2nn = scr * scr
            Z2_j      = 12.0_cd * mp%scr_2nn

            if(C_iij.eq.0.0) then
               tmp_ij = 0.5_cd*(dsqrt(mp%Cmin(i,i,i))+
     &                         dsqrt(mp%Cmin(j,j,j)))
               mp%Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,i,j) = C_iij
            endif
            mp%Cmin(j,i,i) = mp%Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij = 0.5_cd*(dsqrt(mp%Cmin(i,i,i))+
     &                         dsqrt(mp%Cmin(j,j,j)))
               mp%Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,j,j) = C_ijj
            endif
            mp%Cmin(j,j,i) = mp%Cmin(i,j,j)
c
         elseif(lattice.eq.'BCC_B2') then
c
c        for bcc_B2
c
            a2N     = dsqrt(4.0_cd/3.0_cd)
            mp%Ec(i,j) = 0.5_cd * (mp%Ec(i,i) + mp%Ec(j,j)) - delta_Ec
            if(Re_input.eq.0.0_cd) then
               mp%omega(i,j) = 0.5_cd * (mp%omega(i,i) + mp%omega(j,j))
                   tmp_ij = 2.0_cd * mp%omega(i,j)
               Re_ij      = tmp_ij**(1.0_cd/3.0_cd) / a2N
            else
               Re_ij      =  Re_input
                   tmp_ij =  a2N * Re_ij
               mp%omega(i,j) = 0.5_cd * tmp_ij * tmp_ij * tmp_ij
            endif
            if(B_ij.eq.0.0_cd) then
               mp%B(i,j) = 0.5_cd * mp%B(i,i) + 0.5_cd * mp%B(j,j)
            else
               mp%B(i,j) =  B_ij
            endif
            alpha_ij = dsqrt(9.0_cd*mp%B(i,j)*mp%omega(i,j)/mp%Ec(i,j))

            if(d_input.eq.0.0_cd) then
               d_ij   = 0.5_cd * mp%Eu_d(i) + 0.5_cd * mp%Eu_d(j)
            elseif(d_input.lt.0.0_cd) then
               d_ij   =  0.0_cd
            else
               d_ij   =  d_input
            endif

            Z_ij      = 8.0_cd
             scr = screen(0.75_cd,0.75_cd,mp%Cmax(i,j,i),mp%Cmin(i,j,i))
              mp%scr_2nn = scr * scr * scr * scr
            Z2_i      = 6.0_cd * mp%scr_2nn
             scr = screen(0.75_cd,0.75_cd,mp%Cmax(j,i,j),mp%Cmin(j,i,j))
              mp%scr_2nn = scr * scr * scr * scr
            Z2_j      = 6.0_cd * mp%scr_2nn

            if(C_iij.eq.0.0) then
               tmp_ij = 0.5_cd*(dsqrt(mp%Cmin(i,i,i))+
     &                         dsqrt(mp%Cmin(j,j,j)))
               mp%Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,i,j) = C_iij
            endif
            mp%Cmin(j,i,i) = mp%Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij = 0.5_cd*(dsqrt(mp%Cmin(i,i,i))+
     &                         dsqrt(mp%Cmin(j,j,j)))
               mp%Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,j,j) = C_ijj
            endif
            mp%Cmin(j,j,i) = mp%Cmin(i,j,j)
c
         elseif(lattice.eq.'L12AB3') then
c
c        for fcc_L12 (AB3 type)
c
            a2N     = dsqrt(2.0_cd)
            mp%Ec(i,j) = 0.25_cd * mp%Ec(i,i) + 0.75_cd * mp%Ec(j,j)
     &                 - delta_Ec
            if(Re_input.eq.0.0_cd) then
               mp%omega(i,j) = 0.25_cd * mp%omega(i,i)
     &                       + 0.75_cd * mp%omega(j,j)
                   tmp_ij = 4.0_cd * mp%omega(i,j)
               Re_ij      = tmp_ij**(1.0_cd/3.0_cd) / a2N
            else
               Re_ij      =  Re_input
                   tmp_ij =  a2N * Re_ij
               mp%omega(i,j) = 0.25_cd * tmp_ij * tmp_ij * tmp_ij
            endif
            if(B_ij.eq.0.0_cd) then
               mp%B(i,j) = 0.25_cd * mp%B(i,i) + 0.75_cd * mp%B(j,j)
            else
               mp%B(i,j) =  B_ij
            endif
            alpha_ij = dsqrt(9.0_cd*mp%B(i,j)*mp%omega(i,j)/mp%Ec(i,j))

            if(d_input.eq.0.0_cd) then
               d_ij   = 0.25_cd * mp%Eu_d(i) + 0.75_cd * mp%Eu_d(j)
            elseif(d_input.lt.0.0_cd) then
               d_ij   =  0.0_cd
            else
               d_ij   =  d_input
            endif

            Z_ij      = 12.0_cd
              scr = screen(0.5_cd,0.5_cd,mp%Cmax(i,j,i),mp%Cmin(i,j,i))
              mp%scr_2nn = scr * scr * scr * scr
            Z2_i      = 6.0_cd * mp%scr_2nn
c  Here, Cmin(j,i,j) should be equal to Cmin(j,j,j) .. 2002. 6. 14, B.-J. Lee
c  This limitation was removed .. 2011. 4. 5, B.-J. Lee
              scr = screen(0.5_cd,0.5_cd,mp%Cmax(j,j,j),mp%Cmin(j,j,j))
              scr_pure = scr * scr * scr * scr
              scr2 = screen(0.5_cd,0.5_cd,mp%Cmax(j,i,j),mp%Cmin(j,i,j))
              scr_mix = scr * scr * scr2 * scr2
            Z2_j      = 2.0_cd * scr_pure + 4.0_cd * scr_mix

            if(C_iij.eq.0.0) then
               tmp_ij =0.5_cd*dsqrt(mp%Cmin(i,i,i))
     &                +0.5_cd*dsqrt(mp%Cmin(j,j,j))
               mp%Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,i,j) = C_iij
            endif
            mp%Cmin(j,i,i) = mp%Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij =0.5_cd*dsqrt(mp%Cmin(i,i,i))
     &                +0.5_cd*dsqrt(mp%Cmin(j,j,j))
               mp%Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,j,j) = C_ijj
            endif
            mp%Cmin(j,j,i) = mp%Cmin(i,j,j)
c
         elseif(lattice.eq.'L12A3B') then
c
c        for fcc_L12 (A3B type)
c
            a2N     = dsqrt(2.0_cd)
            mp%Ec(i,j) = 0.75_cd * mp%Ec(i,i) + 0.25_cd * mp%Ec(j,j)
     &                 - delta_Ec
            if(Re_input.eq.0.0_cd) then
               mp%omega(i,j) = 0.75_cd * mp%omega(i,i)
     &                       + 0.25_cd * mp%omega(j,j)
                   tmp_ij = 4.0_cd * mp%omega(i,j)
               Re_ij      = tmp_ij**(1.0_cd/3.0_cd) / a2N
            else
               Re_ij      =  Re_input
                   tmp_ij =  a2N * Re_ij
               mp%omega(i,j) = 0.25_cd * tmp_ij * tmp_ij * tmp_ij
            endif
            if(B_ij.eq.0.0_cd) then
               mp%B(i,j) = 0.75_cd * mp%B(i,i) + 0.25_cd * mp%B(j,j)
            else
               mp%B(i,j) =  B_ij
            endif
            alpha_ij = dsqrt(9.0_cd*mp%B(i,j)*mp%omega(i,j)/mp%Ec(i,j))

            if(d_input.eq.0.0_cd) then
               d_ij   = 0.75_cd * mp%Eu_d(i) + 0.25_cd * mp%Eu_d(j)
            elseif(d_input.lt.0.0_cd) then
               d_ij   =  0.0_cd
            else
               d_ij   =  d_input
            endif

            Z_ij      = 12.0_cd
              scr = screen(0.5_cd,0.5_cd,mp%Cmax(j,i,j),mp%Cmin(j,i,j))
              mp%scr_2nn = scr * scr * scr * scr
            Z2_j      = 6.0_cd * mp%scr_2nn
c  Here, Cmin(i,j,i) should be equal to Cmin(i,i,i) .. 2002. 6. 14, B.-J. Lee
c  This limitation was removed .. 2011. 4. 5, B.-J. Lee
              scr = screen(0.5_cd,0.5_cd,mp%Cmax(i,i,i),mp%Cmin(i,i,i))
              scr_pure = scr * scr * scr * scr
              scr2 = screen(0.5_cd,0.5_cd,mp%Cmax(i,j,i),mp%Cmin(i,j,i))
              scr_mix = scr * scr * scr2 * scr2
            Z2_i      = 2.0_cd * scr_pure + 4.0_cd * scr_mix

            if(C_iij.eq.0.0) then
               tmp_ij =0.5_cd*dsqrt(mp%Cmin(i,i,i))
     &                +0.5_cd*dsqrt(mp%Cmin(j,j,j))
               mp%Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,i,j) = C_iij
            endif
            mp%Cmin(j,i,i) = mp%Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij =0.5_cd*dsqrt(mp%Cmin(i,i,i))
     &                +0.5_cd*dsqrt(mp%Cmin(j,j,j))
               mp%Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,j,j) = C_ijj
            endif
            mp%Cmin(j,j,i) = mp%Cmin(i,j,j)
c
         elseif(lattice.eq.'ZnS_B3') then
c
c          for ZnS_B3
c
            a2N     = 4.0_cd / dsqrt(6.0_cd)
            mp%Ec(i,j) = 0.5_cd * (mp%Ec(i,i) + mp%Ec(j,j)) - delta_Ec
            if(Re_input.eq.0.0_cd) then
               mp%omega(i,j) = 0.5_cd * (mp%omega(i,i) + mp%omega(j,j))
                   tmp_ij = 8.0_cd * mp%omega(i,j)
               Re_ij      = tmp_ij**(1.0_cd/3.0_cd) * 
     &                      dsqrt(3.0_cd) / 4.0_cd
            else
               Re_ij      =  Re_input
                   tmp_ij =  4.0_cd * Re_ij / dsqrt(3.0_cd)
               mp%omega(i,j) = 0.125_cd * tmp_ij * tmp_ij * tmp_ij
            endif
            if(B_ij.eq.0.0_cd) then
               mp%B(i,j) = 0.5_cd * mp%B(i,i) + 0.5_cd * mp%B(j,j)
            else
               mp%B(i,j) =  B_ij
            endif
            alpha_ij  = dsqrt(9.0_cd*mp%B(i,j)*mp%omega(i,j)/
     &                  mp%Ec(i,j))

            if(d_input.eq.0.0_cd) then
               d_ij   = 0.5_cd * mp%Eu_d(i) + 0.5_cd * mp%Eu_d(j)
            elseif(d_input.lt.0.0_cd) then
               d_ij   =  0.0_cd
            else
               d_ij   =  d_input
            endif

            Z_ij      = 4.0_cd
              mp%scr_2nn = screen(0.375_cd,0.375_cd,mp%Cmax(i,j,i),
     &                         mp%Cmin(i,j,i))
            Z2_i      = 12.0_cd * mp%scr_2nn
              mp%scr_2nn = screen(0.375_cd,0.375_cd,mp%Cmax(j,i,j),
     &                         mp%Cmin(j,i,j))
            Z2_j      = 12.0_cd * mp%scr_2nn

            if(C_iij.eq.0.0) then
               tmp_ij = 0.5_cd*(dsqrt(mp%Cmin(i,i,i))+
     &                         dsqrt(mp%Cmin(j,j,j)))
               mp%Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,i,j) = C_iij
            endif
            mp%Cmin(j,i,i) = mp%Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij = 0.5_cd*(dsqrt(mp%Cmin(i,i,i))+
     &                         dsqrt(mp%Cmin(j,j,j)))
               mp%Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,j,j) = C_ijj
            endif
            mp%Cmin(j,j,i) = mp%Cmin(i,j,j)
c
         elseif(lattice.eq.'DIMER') then
c
c          for Diatomic Gases
c
            a2N     = dsqrt(2.0_cd)
            mp%Ec(i,j) = 0.5_cd * (mp%Ec(i,i) + mp%Ec(j,j)) - delta_Ec
            if(Re_input.eq.0.0_cd) then
               mp%omega(i,j) = 0.5_cd * (mp%omega(i,i) + mp%omega(j,j))
               Re_ij      = mp%omega(i,j)**(1.0_cd/3.0_cd)
            else
               Re_ij      =  Re_input
               mp%omega(i,j) = Re_ij * Re_ij * Re_ij
            endif
            if(B_ij.eq.0.0_cd) then
               mp%B(i,j) = 0.5_cd * mp%B(i,i) + 0.5_cd * mp%B(j,j)
            else
               mp%B(i,j) =  B_ij
            endif
            alpha_ij  = dsqrt(9.0_cd*mp%B(i,j)*mp%omega(i,j)/mp%Ec(i,j))

            if(d_input.eq.0.0_cd) then
               d_ij   = 0.5_cd * mp%Eu_d(i) + 0.5_cd * mp%Eu_d(j)
            elseif(d_input.lt.0.0_cd) then
               d_ij   =  0.0_cd
            else
               d_ij   =  d_input
            endif

            Z_ij      = 1.0_cd
              mp%scr_2nn = 0.0_cd
            Z2_i      = mp%scr_2nn
            Z2_j      = mp%scr_2nn

            if(C_iij.eq.0.0) then
               tmp_ij = 0.5_cd*(dsqrt(mp%Cmin(i,i,i))+
     &                         dsqrt(mp%Cmin(j,j,j)))
               mp%Cmin(i,i,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,i,j) = C_iij
            endif
            mp%Cmin(j,i,i) = mp%Cmin(i,i,j)
            if(C_ijj.eq.0.0) then
               tmp_ij = 0.5_cd*(dsqrt(mp%Cmin(i,i,i))+
     &                         dsqrt(mp%Cmin(j,j,j)))
               mp%Cmin(i,j,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,j,j) = C_ijj
            endif
            mp%Cmin(j,j,i) = mp%Cmin(i,j,j)
c
         else
            goto 4445
         endif
c
         mp%Re_Alloy(i,j) = Re_ij
         mp%alpha(i,j) = alpha_ij
         write(*,'(1x,a5,a7,f7.4,a10,f7.4)')
     &   alloy_1_name,':  Re =',Re_ij,',  alpha =',alpha_ij
c
c -----------------------------------------------------------------------------
c
c  Make a table for pair potential and its derivative between different element
c
         do k=1,mp%LineTable+1
            Dij = mp%RMin + (k-1) * DeltaR
c
c        Calculation of Eu(R) and dEdR for reference structure at r = Dij
c
            Rfactor = Dij / Re_ij - 1.0_cd
            astar   = alpha_ij * Rfactor
            astar3  = astar * astar * astar
            Eu = - mp%Ec(i,j) * (1.0_cd + astar + d_ij*astar3)
     &           * dexp(-astar)
            dEdR = mp%Ec(i,j) * astar * alpha_ij / Re_ij * dexp(-astar)
     &         * (1.0_cd + d_ij * astar * astar - 3.0_cd * d_ij * astar)
c
c        Calculation of Phi and DPhi at r = a2N * Dij for each element
c
            if(lattice(1:3).ne.'DIM') then
               rr = a2N * Dij
            ZR2 = mp%Z2(i) / mp%Z(i)
            alphaN = mp%alpha(i,i)
            if(lattice_element(i).eq.'FCC_A1') then
               scr = screen(0.5_cd,0.5_cd,mp%Cmax(i,i,i),mp%Cmin(i,i,i))
               mp%scr_2nn = scr * scr * scr * scr
            elseif(lattice_element(i).eq.'BCC_A2') then
               scr = screen(0.75_cd,0.75_cd,mp%Cmax(i,i,i),
     &                      mp%Cmin(i,i,i))
               mp%scr_2nn = scr * scr * scr * scr
            elseif(lattice_element(i).eq.'HCP_A3') then
               scr = screen(0.5_cd,0.5_cd,mp%Cmax(i,i,i),mp%Cmin(i,i,i))
               mp%scr_2nn = scr * scr * scr * scr
            elseif(lattice_element(i).eq.'DIA_A4') then
               mp%scr_2nn = screen(0.375_cd,0.375_cd,mp%Cmax(i,i,i),
     &                          mp%Cmin(i,i,i))
            endif
            call Define_Pair_Potential(mp,i,lattice_element(i),alphaN,
     &                                 ZR2,mp%a2(i),rr,Phi_i,DPhi_i)
            ZR2 = mp%Z2(j) / mp%Z(j)
            alphaN = mp%alpha(j,j)
            if(lattice_element(j).eq.'FCC_A1') then
               scr = screen(0.5_cd,0.5_cd,mp%Cmax(j,j,j),mp%Cmin(j,j,j))
               mp%scr_2nn = scr * scr * scr * scr
            elseif(lattice_element(j).eq.'BCC_A2') then
               scr = screen(0.75_cd,0.75_cd,mp%Cmax(j,j,j),
     &                      mp%Cmin(j,j,j))
               mp%scr_2nn = scr * scr * scr * scr
            elseif(lattice_element(j).eq.'HCP_A3') then
               scr = screen(0.5_cd,0.5_cd,mp%Cmax(j,j,j),mp%Cmin(j,j,j))
               mp%scr_2nn = scr * scr * scr * scr
            elseif(lattice_element(j).eq.'DIA_A4') then
               mp%scr_2nn = screen(0.375_cd,0.375_cd,mp%Cmax(j,j,j),
     &                             mp%Cmin(j,j,j))
            endif
            call Define_Pair_Potential(mp,j,lattice_element(j),alphaN,
     &                                 ZR2,mp%a2(j),rr,Phi_j,DPhi_j)
            endif
c
            R12factor = dexp(- mp%beta(j,4) * (Dij / mp%Re(j) - 1.0_cd))
            R21factor = dexp(- mp%beta(i,4) * (a2N*Dij/mp%Re(i)
     &                                                        - 1.0_cd))
            R11factor = dexp(- mp%beta(i,4) * (Dij / mp%Re(i) - 1.0_cd))
            R22factor = dexp(- mp%beta(j,4) * (a2N*Dij/mp%Re(j)
     &                                                        - 1.0_cd))
c
         if(lattice.eq.'FCC_B1' .or. lattice.eq.'NaCl' .or.
     &      lattice.eq.'BCC_B2') then
c
c           calculation of Rho_Bar, F, dFdR for reference structure at r = Dij
c
            Rho_Bar_i = Z_ij * mp%RhoZ(j) * R12factor
     &                + Z2_i * mp%RhoZ(i) * R21factor
            Rho_Bar_j = Z_ij * mp%RhoZ(i) * R11factor
     &                + Z2_j * mp%RhoZ(j) * R22factor
c
            if(rho_bar_i.eq.0.0_cd) rho_bar_i = 1.0e-20_cd
            if(rho_bar_j.eq.0.0_cd) rho_bar_j = 1.0e-20_cd
            rr0_i = Rho_Bar_i / mp%Rho0_Bar(i)
            rr0_j = Rho_Bar_j / mp%Rho0_Bar(j)
            F_i = mp%A(i) * mp%Ec(i,i) * rr0_i * dlog(rr0_i)
            F_j = mp%A(j) * mp%Ec(j,j) * rr0_j * dlog(rr0_j)
c
c           dFdR = - beta(4)/Re * A * Ec * rr0 * ( dlog(rr0) + 1.0_cd )
c
            dRho_i = mp%beta(j,4) / mp%Re(j) * Z_ij * mp%RhoZ(j)
     &                                              * R12factor
     &       + a2N * mp%beta(i,4) / mp%Re(i) * Z2_i * mp%RhoZ(i)
     &                                              * R21factor
            dRho_j = mp%beta(i,4) / mp%Re(i) * Z_ij * mp%RhoZ(i)
     &                                              * R11factor
     &       + a2N * mp%beta(j,4) / mp%Re(j) * Z2_j * mp%RhoZ(j)
     &                                              * R22factor
            dFdR_i = - mp%A(i) * mp%Ec(i,i) / mp%Rho0_Bar(i)
     &             * ( dlog(rr0_i) + 1.0_cd ) * dRho_i
            dFdR_j = - mp%A(j) * mp%Ec(j,j) / mp%Rho0_Bar(j)
     &             * ( dlog(rr0_j) + 1.0_cd ) * dRho_j
c
c           Pair potential and its derivative between i and j atoms at r = Dij
c
            mp%PhiTab(k,i,j) = 1.0_cd / Z_ij *
     &        (2.0_cd*Eu - F_i - F_j - 0.5_cd*(Z2_i*Phi_i + Z2_j*Phi_j))
            mp%PhiTab(k,j,i) = mp%PhiTab(k,i,j)
            mp%DPhiTab(k,i,j) = (2.0_cd*dEdR -dFdR_i -dFdR_j) /Z_ij/ Dij
     &          - 0.5_cd * a2N * a2N * (Z2_i*DPhi_i+Z2_j*DPhi_j) / Z_ij
            mp%DPhiTab(k,j,i) = mp%DPhiTab(k,i,j)
c
         elseif(lattice.eq.'L12AB3') then
c
c           calculation of Rho_Bar and F for reference structure at r = Dij
c
            Rho_Bar_i = Z_ij * mp%RhoZ(j) * R12factor
     &                + Z2_i * mp%RhoZ(i) * R21factor

            Rho0 = 4.0_cd * mp%RhoZ(i) * R11factor
     &           + 8.0_cd * mp%RhoZ(j) * R12factor
     &           + Z2_j * mp%RhoZ(j) * R22factor
            t_1i = mp%RhoZ(i) * dexp( - mp%beta(i,2)
     &                     * (Dij/mp%Re(i) - 1.0_cd) )
            t_1j = mp%RhoZ(j) * dexp( - mp%beta(j,2)
     &                     * (Dij/mp%Re(j) - 1.0_cd) )
            term2 = mp%RhoZ(j) * dexp( - mp%beta(j,2)
     &                      * (a2N*Dij/mp%Re(j)-1.0_cd) )
            t_2i = scr_mix * term2
            t_2j = scr_pure * term2
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.0_cd) RhoFactor = 1.e-20_cd
            Rho2_1 = 24.0_cd * (t_1i + t_1j + t_2i) * 
     &               (t_1i + t_1j + t_2i)
            Rho2_2 = 12.0_cd * (2.0_cd * t_1j + t_2j) * 
     &               (2.0_cd * t_1j + t_2j)
            Rho2_3 = 8.0_cd*t_1j + 4.0_cd*t_1i + 2.0_cd*t_2j + 
     &               4.0_cd*t_2i
            Rho2P = ( Rho2_1 + Rho2_2 - Rho2_3 * Rho2_3 ) / 3.0_cd
            GAMMA = mp%t(j,2) * Rho2P / RhoFactor
            Rho_Bar_j = Rho0 * 2.0_cd / ( 1.0_cd + dexp(-GAMMA) )
c
            if(rho_bar_i.eq.0.0_cd) rho_bar_i = 1.0e-20_cd
            if(rho_bar_j.eq.0.0_cd) rho_bar_j = 1.0e-20_cd
            rr0_i = Rho_Bar_i / mp%Rho0_Bar(i)
            rr0_j = Rho_Bar_j / mp%Rho0_Bar(j)
            F_i = mp%A(i) * mp%Ec(i,i) * rr0_i * dlog(rr0_i)
            F_j = mp%A(j) * mp%Ec(j,j) * rr0_j * dlog(rr0_j)
c
c           Pair potential and its derivative between i and j atoms at r = Dij
c
            mp%PhiTab(k,i,j) = 1.0_cd / Z_ij *
     &        ( 4.0_cd * Eu - F_i - 3.0_cd*F_j - Z_ij * mp%PhiTab(k,j,j)
     &          - 0.5_cd * (Z2_i * Phi_i + 3.0_cd * Z2_j * Phi_j) )
            mp%PhiTab(k,j,i) = mp%PhiTab(k,i,j)
            if(k.gt.2) then
               mp%DPhiTab(k-1,i,j) = 0.5_cd*(mp%PhiTab(k,i,j)
     &                             - mp%PhiTab(k-2,i,j))
     &                          / DeltaR / (Dij-DeltaR)
               mp%DPhiTab(k-1,j,i) = mp%DPhiTab(k-1,i,j)
            endif
c
         elseif(lattice.eq.'L12A3B') then
c
c           calculation of Rho_Bar and F for reference structure at r = Dij
c
            Rho_Bar_j = Z_ij * mp%RhoZ(i) * R11factor
     &                + Z2_j * mp%RhoZ(j) * R22factor

            Rho0 = 4.0_cd * mp%RhoZ(j) * R12factor
     &           + 8.0_cd * mp%RhoZ(i) * R11factor
     &           + Z2_i * mp%RhoZ(i) * R21factor
            t_1i = mp%RhoZ(i) * dexp( - mp%beta(i,2)
     &                     * (Dij/mp%Re(i) - 1.0_cd) )
            t_1j = mp%RhoZ(j) * dexp( - mp%beta(j,2)
     &                     * (Dij/mp%Re(j) - 1.0_cd) )
            term2 = mp%RhoZ(i) * dexp( - mp%beta(i,2)
     &                      * (a2N*Dij/mp%Re(i)-1.0_cd) )
            t_2i = scr_pure * term2
            t_2j = scr_mix * term2
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.0_cd) RhoFactor = 1.e-20_cd
            Rho2_1 = 24.0_cd * (t_1i + t_1j + t_2j) * 
     &               (t_1i + t_1j + t_2j)
            Rho2_2 = 12.0_cd * (2.0_cd * t_1i + t_2i) * 
     &               (2.0_cd * t_1i + t_2i)
            Rho2_3 = 8.0_cd*t_1i + 4.0_cd*t_1j + 2.0_cd*t_2i + 
     &               4.0_cd*t_2j
            Rho2P = ( Rho2_1 + Rho2_2 - Rho2_3 * Rho2_3 ) / 3.0_cd
            GAMMA = mp%t(i,2) * Rho2P / RhoFactor
            Rho_Bar_i = Rho0 * 2.0_cd / ( 1.0_cd + dexp(-GAMMA) )
c
            if(rho_bar_i.eq.0.0_cd) rho_bar_i = 1.0e-20_cd
            if(rho_bar_j.eq.0.0_cd) rho_bar_j = 1.0e-20_cd
            rr0_i = Rho_Bar_i / mp%Rho0_Bar(i)
            rr0_j = Rho_Bar_j / mp%Rho0_Bar(j)
            F_i = mp%A(i) * mp%Ec(i,i) * rr0_i * dlog(rr0_i)
            F_j = mp%A(j) * mp%Ec(j,j) * rr0_j * dlog(rr0_j)
c
c           Pair potential and its derivative between i and j atoms at r = Dij
c
            mp%PhiTab(k,i,j) = 1.0_cd / Z_ij *
     &        ( 4.0_cd * Eu - 3.0_cd*F_i - F_j - Z_ij * mp%PhiTab(k,i,i)
     &          - 0.5_cd * (Z2_j * Phi_j + 3.0_cd * Z2_i * Phi_i) )
            mp%PhiTab(k,j,i) = mp%PhiTab(k,i,j)
            if(k.gt.2) then
               mp%DPhiTab(k-1,i,j) = 0.5_cd*(mp%PhiTab(k,i,j)
     &                          - mp%PhiTab(k-2,i,j))
     &                          / DeltaR / (Dij-DeltaR)
               mp%DPhiTab(k-1,j,i) = mp%DPhiTab(k-1,i,j)
            endif
c
         elseif(lattice.eq.'ZnS_B3') then
c
c           calculation of Rho_Bar and F for reference structure at r = Dij
c
            Rho0 = Z_ij * mp%RhoZ(j) * R12factor
     &           + Z2_i * mp%RhoZ(i) * R21factor
            Rho3 = mp%RhoZ(j) * dexp( - mp%beta(j,3)
     &                     * (Dij/mp%Re(j) - 1.0_cd) )
            Rho3P = 32.0_cd * Rho3 * Rho3 / 9.0_cd
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.0_cd) RhoFactor = 1.e-20_cd
            GAMMA = mp%t(i,3) * Rho3P / RhoFactor
            Rho_Bar_i = Rho0 * 2.0_cd / ( 1.0_cd + dexp(-GAMMA) )

            Rho0 = Z_ij * mp%RhoZ(i) * R11factor
     &           + Z2_j * mp%RhoZ(j) * R22factor
            Rho3 = mp%RhoZ(i) * dexp( - mp%beta(i,3)
     &                     * (Dij/mp%Re(i) - 1.0_cd) )
            Rho3P = 32.0_cd * Rho3 * Rho3 / 9.0_cd
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.0_cd) RhoFactor = 1.e-20_cd
            GAMMA = mp%t(j,3) * Rho3P / RhoFactor
            Rho_Bar_j = Rho0 * 2.0_cd / ( 1.0_cd + dexp(-GAMMA) )
c
            if(rho_bar_i.eq.0.0_cd) rho_bar_i = 1.0e-20_cd
            if(rho_bar_j.eq.0.0_cd) rho_bar_j = 1.0e-20_cd
            rr0_i = Rho_Bar_i / mp%Rho0_Bar(i)
            rr0_j = Rho_Bar_j / mp%Rho0_Bar(j)
            F_i = mp%A(i) * mp%Ec(i,i) * rr0_i * dlog(rr0_i)
            F_j = mp%A(j) * mp%Ec(j,j) * rr0_j * dlog(rr0_j)
c
c           Pair potential and its derivative between i and j atoms at r = Dij
c
            mp%PhiTab(k,i,j) = 1.0_cd / Z_ij *
     &        (2.0_cd*Eu - F_i - F_j - 0.5_cd*(Z2_i*Phi_i + Z2_j*Phi_j))
            mp%PhiTab(k,j,i) = mp%PhiTab(k,i,j)
            if(k.gt.2) then
               mp%DPhiTab(k-1,i,j) = 0.5_cd*(mp%PhiTab(k,i,j)
     &                          - mp%PhiTab(k-2,i,j))
     &                          / DeltaR / (Dij-DeltaR)
               mp%DPhiTab(k-1,j,i) = mp%DPhiTab(k-1,i,j)
            endif
c
         elseif(lattice.eq.'DIMER') then
c
c           calculation of Rho_Bar and F for reference structure at r = Dij
c
            Rho0 = Z_ij * mp%RhoZ(j) * R12factor
            Rho1 = mp%RhoZ(j) * dexp( - mp%beta(j,1)
     &                     * (Dij/mp%Re(j) - 1.0_cd) )
            Rho2 = mp%RhoZ(j) * dexp( - mp%beta(j,2)
     &                     * (Dij/mp%Re(j) - 1.0_cd) )
            Rho3 = mp%RhoZ(j) * dexp( - mp%beta(j,3)
     &                     * (Dij/mp%Re(j) - 1.0_cd) )
            Rho1P = Rho1 * Rho1
            Rho2P = 2.0_cd * Rho2 * Rho2 / 3.0_cd
            Rho3P = 2.0_cd * Rho3 * Rho3 / 5.0_cd
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.0_cd) RhoFactor = 1.e-20_cd
            GAMMA = (mp%t(i,1)*Rho1P+mp%t(i,2)*Rho2P+mp%t(i,3)*Rho3P)
     &            / RhoFactor
            Rho_Bar_i = Rho0 * 2.0_cd / ( 1.0_cd + dexp(-GAMMA) )
c           Rho_Bar_i = Rho0 * dsqrt ( 1.0_cd + GAMMA )

            Rho0 = Z_ij * mp%RhoZ(i) * R11factor
            Rho1 = mp%RhoZ(i) * dexp( - mp%beta(i,1)
     &                     * (Dij/mp%Re(i) - 1.0_cd) )
            Rho2 = mp%RhoZ(i) * dexp( - mp%beta(i,2)
     &                     * (Dij/mp%Re(i) - 1.0_cd) )
            Rho3 = mp%RhoZ(i) * dexp( - mp%beta(i,3)
     &                     * (Dij/mp%Re(i) - 1.0_cd) )
            Rho1P = Rho1 * Rho1
            Rho2P = 2.0_cd * Rho2 * Rho2 / 3.0_cd
            Rho3P = 2.0_cd * Rho3 * Rho3 / 5.0_cd
            RhoFactor = Rho0 * Rho0
            if(RhoFactor.eq.0.0_cd) RhoFactor = 1.e-20_cd
            GAMMA = (mp%t(j,1)*Rho1P+mp%t(j,2)*Rho2P+mp%t(j,3)*Rho3P)
     &            / RhoFactor
            Rho_Bar_j = Rho0 * 2.0_cd / ( 1.0_cd + dexp(-GAMMA) )
c           Rho_Bar_j = Rho0 * dsqrt ( 1.0_cd + GAMMA )
c
            if(rho_bar_i.eq.0.0_cd) rho_bar_i = 1.0e-20_cd
            if(rho_bar_j.eq.0.0_cd) rho_bar_j = 1.0e-20_cd
            rr0_i = Rho_Bar_i / mp%Rho0_Bar(i)
            rr0_j = Rho_Bar_j / mp%Rho0_Bar(j)
            F_i = mp%A(i) * mp%Ec(i,i) * rr0_i * dlog(rr0_i)
            F_j = mp%A(j) * mp%Ec(j,j) * rr0_j * dlog(rr0_j)
c
c           Pair potential and its derivative between i and j atoms at r = Dij
c
            mp%PhiTab(k,i,j) = (2.0_cd*Eu - F_i - F_j) / Z_ij
            mp%PhiTab(k,j,i) = mp%PhiTab(k,i,j)
            if(k.gt.2) then
               mp%DPhiTab(k-1,i,j) = 0.5_cd*(mp%PhiTab(k,i,j)
     &                          - mp%PhiTab(k-2,i,j))
     &                          / DeltaR / (Dij-DeltaR)
               mp%DPhiTab(k-1,j,i) = mp%DPhiTab(k-1,i,j)
            endif
c
c        elseif(lattice.eq.'HCP_A3' .or. lattice.eq.'PtCo') then
c     Description for PtCo type ordered HCP structure should be given
c     2000. 5. 4.  B.-J. Lee
c
         endif
         enddo
      enddo
      enddo
      endif
c
      if(Ncomp.gt.2) then
c
         do i = 1, Ncomp-2
         do j = i+1, Ncomp-1
         do k = j+1, Ncomp
c
            ternary_name(1:2) = Icomp(i)
            ternary_name(3:3) = '-'
            ternary_name(4:5) = Icomp(j)
            ternary_name(6:6) = '-'
            ternary_name(7:8) = Icomp(k)
c
            read(2,'(a)') card
            read(2,'(a)') card2
    3       continue
            read(2,'(a)',end=4447) card
            read(card,'(a8)') alloy3
            if(alloy3.ne.ternary_name) goto 3

            read(card,2011) Cx_ikj, Cx_ijk, Cx_jik, C_ikj, C_ijk, C_jik
 2011       format(40x,6f5.2)

            if(C_ikj.eq.0.0) then
               tmp_ij = 0.5_cd*(dsqrt(mp%Cmin(i,i,i))+
     &                         dsqrt(mp%Cmin(j,j,j)))
               mp%Cmin(i,k,j) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,k,j) = C_ikj
            endif
            mp%Cmin(j,k,i) = mp%Cmin(i,k,j)

            if(C_ijk.eq.0.0) then
               tmp_ij = 0.5_cd*(dsqrt(mp%Cmin(i,i,i))+
     &                         dsqrt(mp%Cmin(k,k,k)))
               mp%Cmin(i,j,k) = tmp_ij * tmp_ij
            else
               mp%Cmin(i,j,k) = C_ijk
            endif
            mp%Cmin(k,j,i) = mp%Cmin(i,j,k)

            if(C_jik.eq.0.0) then
               tmp_ij = 0.5_cd*(dsqrt(mp%Cmin(j,j,j))+
     &                         dsqrt(mp%Cmin(k,k,k)))
               mp%Cmin(j,i,k) = tmp_ij * tmp_ij
            else
               mp%Cmin(j,i,k) = C_jik
            endif
            mp%Cmin(k,i,j) = mp%Cmin(j,i,k)

            if(Cx_ikj.ne.0.0) then
               mp%Cmax(i,k,j) = Cx_ikj
               mp%Cmax(j,k,i) = Cx_ikj
            endif
            if(Cx_ijk.ne.0.0) then
               mp%Cmax(i,j,k) = Cx_ijk
               mp%Cmax(k,j,i) = Cx_ijk
            endif
            if(Cx_jik.ne.0.0) then
               mp%Cmax(j,i,k) = Cx_jik
               mp%Cmax(k,i,j) = Cx_jik
            endif
c
         enddo
         enddo
         enddo
      endif
c
      close(2)
      Ierr = KIM_STATUS_OK
      return
c
 4444 continue
      close(2)
      Ierr = KIM_STATUS_FAIL
      write(*,'(/,a,a,/)') ' **** ERROR : Non-Registered Element ',
     &      Icomp(n)
      return
 4445 continue
      close(2)
      Ierr = KIM_STATUS_FAIL
      write(*,'(/,a,a,/)') ' **** ERROR : Non-Registered Structure ',
     &      lattice
      return
 4446 continue
      close(2)
      Ierr = KIM_STATUS_FAIL
      write(*,'(/,a,a,/)') ' **** ERROR : Non-Registered Alloy ',
     &      alloy_1_name
      return
 4447 continue
      close(2)
      Ierr = KIM_STATUS_FAIL
      write(*,'(/,a,a,/)') ' **** ERROR : Non-Registered Alloy ',
     &      ternary_name
      return
 4448 continue
      close(2)
      Ierr = KIM_STATUS_FAIL
      write(*,'(/,a,/)') ' **** ERROR : Unable to open file '
      return
      end subroutine Define_EAM_Potential


      subroutine Define_Pair_Potential
     &          (mp,N,lattice,alphaN,ZR2,a2N,Dij,Phi,DPhi)
c
c  Compute phi(R) and dphi(R) for element according to 2nd NN approximation
c
      use, intrinsic :: iso_c_binding
      use meam_parameters
      implicit none

      type(meam_params) :: mp
      integer(c_int) N
      character lattice*6
      real(c_double) alphaN
      real(c_double) ZR2
      real(c_double) a2N
      real(c_double) Dij
      real(c_double) Phi
      real(c_double) DPhi

      real(c_double) P_Phi, P_DPhi
c
      call Compute_Phi(mp,N,lattice,alphaN,1.0_cd,Dij,P_Phi,P_DPhi)
           Phi  = P_Phi
           DPhi = P_DPhi
c
      if(ZR2.eq.0.0_cd) return
c
      call Compute_Phi(mp,N,lattice,alphaN,a2N,Dij,P_Phi,P_DPhi)
           Phi  = Phi  - ZR2 * P_Phi
           DPhi = DPhi - ZR2 * P_DPhi
      call Compute_Phi(mp,N,lattice,alphaN,a2N*a2N,Dij,P_Phi,P_DPhi)
           Phi  = Phi  + ZR2*ZR2 * P_Phi
           DPhi = DPhi + ZR2*ZR2 * P_DPhi
      call Compute_Phi(mp,N,lattice,alphaN,a2N**3,Dij,P_Phi,P_DPhi)
           Phi  = Phi  - ZR2*ZR2*ZR2 * P_Phi
           DPhi = DPhi - ZR2*ZR2*ZR2 * P_DPhi
      call Compute_Phi(mp,N,lattice,alphaN,a2N**4,Dij,P_Phi,P_DPhi)
           Phi  = Phi  + ZR2*ZR2*ZR2*ZR2 * P_Phi
           DPhi = DPhi + ZR2*ZR2*ZR2*ZR2 * P_DPhi
      call Compute_Phi(mp,N,lattice,alphaN,a2N**5,Dij,P_Phi,P_DPhi)
           Phi  = Phi  - ZR2*ZR2*ZR2*ZR2*ZR2 * P_Phi
           DPhi = DPhi - ZR2*ZR2*ZR2*ZR2*ZR2 * P_DPhi
      call Compute_Phi(mp,N,lattice,alphaN,a2N**6,Dij,P_Phi,P_DPhi)
           Phi  = Phi  + ZR2*ZR2*ZR2*ZR2*ZR2*ZR2 * P_Phi
           DPhi = DPhi + ZR2*ZR2*ZR2*ZR2*ZR2*ZR2 * P_DPhi
      call Compute_Phi(mp,N,lattice,alphaN,a2N**7,Dij,P_Phi,P_DPhi)
           Phi  = Phi  - ZR2**7 * P_Phi
           DPhi = DPhi - ZR2**7 * P_DPhi
      call Compute_Phi(mp,N,lattice,alphaN,a2N**8,Dij,P_Phi,P_DPhi)
           Phi  = Phi  + ZR2**8 * P_Phi
           DPhi = DPhi + ZR2**8 * P_DPhi
      call Compute_Phi(mp,N,lattice,alphaN,a2N**9,Dij,P_Phi,P_DPhi)
           Phi  = Phi  - ZR2**9 * P_Phi
           DPhi = DPhi - ZR2**9 * P_DPhi
      call Compute_Phi(mp,N,lattice,alphaN,a2N**10,Dij,P_Phi,P_DPhi)
           Phi  = Phi  + ZR2**10 * P_Phi
           DPhi = DPhi + ZR2**10 * P_DPhi
      call Compute_Phi(mp,N,lattice,alphaN,a2N**11,Dij,P_Phi,P_DPhi)
           Phi  = Phi  - ZR2**11 * P_Phi
           DPhi = DPhi - ZR2**11 * P_DPhi
c
      return
      end


      subroutine Compute_Phi(mp,i,lattice,alphaN,aa,Dij,Phi,DPhi)
c
c  Compute Eu(r) and F(rho(r)) for calculation of phi-potential
c
      use, intrinsic :: iso_c_binding
      use meam_parameters
      implicit none

      type(meam_params) :: mp
      integer(c_int) i
      character lattice*6
      real(c_double) alphaN
      real(c_double) aa
      real(c_double) Dij
      real(c_double) Phi
      real(c_double) DPhi

      real(c_double) Rfactor, astar, astar3, Eu, dEdR, R2factor,
     &   Rho_Bar, dRho0, dFdR, dRho1, dRho2, dRho3, dRho_App, dRho_dR,
     &   EXP_GAMMA, F, G_GAMMA, G_prime, GAMMA, Rho0, Rho1, Rho1P, Rho2,
     &   Rho2P, Rho3, Rho3P, RhoFactor, rr0
c
c  calculation of Eu(R) and dEdR for reference structure at r = aa*Dij
c
      Rfactor = aa * Dij / mp%Re(i) - 1.0_cd
      astar  = alphaN * Rfactor
      astar3 = astar * astar * astar
      Eu = - mp%Ec(i,i) * (1.0_cd + astar + mp%Eu_d(i) * astar3)
     &     * dexp(-astar)
      dEdR = mp%Ec(i,i) * astar * alphaN / mp%Re(i) * dexp(-astar) * aa
     &     * (1.0_cd + mp%Eu_d(i) * astar * astar - 3.0_cd * mp%Eu_d(i)
     &                                                 * astar)
c
c  for calculation of Rho_Bar, F and dFdR for reference structure at r = aa*Dij
c
      R2factor = aa * mp%a2(i) * Dij / mp%Re(i) - 1.0_cd
      Rho_Bar =  mp%Z(i) * mp%RhoZ(i) * dexp(-mp%beta(i,4)*Rfactor)
     &        + mp%Z2(i) * mp%RhoZ(i) * dexp(-mp%beta(i,4)*R2factor)
      dRho0 = - mp%beta(i,4) / mp%Re(i) * aa
     &      * ( mp%Z(i) * mp%RhoZ(i) * dexp(-mp%beta(i,4)*Rfactor)
     &      +  mp%a2(i) * mp%Z2(i) * mp%RhoZ(i)
     &                       * dexp(-mp%beta(i,4)*R2factor) )
      if(Rho_Bar.lt.1.0e-20_cd) then
         Phi   = 0.0_cd
         DPhi  = 0.0_cd
         return
      endif
c
      G_GAMMA  = 1.0_cd
      dRho_App = 0.0_cd
c
c  add t(3) term in the case of HCP_A3 or DIA_A4 lattice
c
      if(lattice.eq.'HCP_A3' .or. lattice.eq.'DIA_A4') then
         if(lattice.eq.'HCP_A3') then
            Rho3  = 0.5_cd * mp%RhoZ(i) * dexp(-mp%beta(i,3)*Rfactor)
     &  - mp%scr_2nn * dsqrt(2.0_cd) * mp%RhoZ(i)
     &                          * dexp(-mp%beta(i,3)*R2factor)
            Rho3P = 4.0_cd * Rho3 * Rho3 / 3.0_cd
            dRho3 = - mp%beta(i,3) / mp%Re(i) * aa * mp%RhoZ(i)
     &            * ( 0.5_cd * dexp(-mp%beta(i,3)*Rfactor)
     &  - mp%scr_2nn * mp%a2(i) * dsqrt(2.0_cd)
     &                       * dexp(-mp%beta(i,3)*R2factor) )
         elseif(lattice.eq.'DIA_A4') then
            Rho3  = mp%RhoZ(i) * dexp(-mp%beta(i,3)*Rfactor)
            Rho3P = 32.0_cd * Rho3 * Rho3 / 9.0_cd
            dRho3 = - mp%beta(i,3) / mp%Re(i) * aa * Rho3
         endif
c
         Rho0 = Rho_Bar
         RhoFactor = Rho0 * Rho0
         if(RhoFactor.eq.0.0_cd) RhoFactor = 1.e-20_cd
         GAMMA = mp%t(i,3) * Rho3P / RhoFactor
            if(GAMMA.lt.-500.0_cd) then
               Phi   = 0.0_cd
               DPhi  = 0.0_cd
               return
            endif
         EXP_GAMMA = dexp(-GAMMA)
         G_GAMMA = 2.0_cd / ( 1.0_cd + EXP_GAMMA )
         Rho_Bar = Rho0 * G_GAMMA
c        Rho_Bar = Rho0 * dsqrt ( 1.0_cd + GAMMA )
c
         G_prime = 2.0_cd*EXP_GAMMA / (1.0_cd+EXP_GAMMA) / 
     &             (1.0_cd+EXP_GAMMA)
         dRho_App = Rho0 * G_prime * 2.0_cd * GAMMA * Rho0 / Rho3
     &            * ( dRho3 / Rho0 - Rho3 / Rho0 / Rho0 * dRho0 )
c        dRho_App = Rho0 / dsqrt( 1.0_cd + GAMMA ) * GAMMA
c    &            * ( dRho3 / Rho3 - dRho0 / Rho0 )
      endif
c
c  add t(1), t(2) and t(3) terms in the case of DIMER species
c
      if(lattice.eq.'DIMER') then
         Rho1  = mp%RhoZ(i) * dexp(-mp%beta(i,1)*Rfactor)
         Rho2  = mp%RhoZ(i) * dexp(-mp%beta(i,2)*Rfactor)
         Rho3  = mp%RhoZ(i) * dexp(-mp%beta(i,3)*Rfactor)
         Rho1P = Rho1 * Rho1
         Rho2P = 2.0_cd * Rho2 * Rho2 / 3.0_cd
         Rho3P = 2.0_cd * Rho3 * Rho3 / 5.0_cd
         dRho1 = - mp%beta(i,1) / mp%Re(i) * aa * Rho1
         dRho2 = - mp%beta(i,2) / mp%Re(i) * aa * Rho2
         dRho3 = - mp%beta(i,3) / mp%Re(i) * aa * Rho3
c
         Rho0 = Rho_Bar
         RhoFactor = Rho0 * Rho0
         if(RhoFactor.eq.0.0_cd) RhoFactor = 1.e-20_cd
         GAMMA = (mp%t(i,1)*Rho1P + mp%t(i,2)*Rho2P + mp%t(i,3)*Rho3P)
     &         / RhoFactor
            if(GAMMA.lt.-500.0_cd) then
               Phi   = 0.0_cd
               DPhi  = 0.0_cd
               return
            endif
         EXP_GAMMA = dexp(-GAMMA)
         G_GAMMA = 2.0_cd / ( 1.0_cd + EXP_GAMMA )
         Rho_Bar = Rho0 * G_GAMMA
c        Rho_Bar = Rho0 * dsqrt ( 1.0_cd + GAMMA )
c
         G_prime = 2.0_cd*EXP_GAMMA / (1.0_cd+EXP_GAMMA) / 
     &             (1.0_cd+EXP_GAMMA)
         dRho_App = Rho0 * G_prime * 2.0_cd *
     1   ( mp%t(i,1)*Rho1/Rho0 * (dRho1*Rho0-Rho1*dRho0) / RhoFactor
     2   + mp%t(i,2)*2.0_cd/3.0_cd*Rho2/Rho0*(dRho2*Rho0-Rho2*dRho0)
     3   / RhoFactor + mp%t(i,3)*0.4_cd*Rho3/Rho0
     4   * (dRho3*Rho0-Rho3*dRho0) / RhoFactor )
c        dRho_App = Rho0 / dsqrt( 1.0_cd + GAMMA ) *
c    1   ( t(i,1)*Rho1/Rho0 * (dRho1*Rho0-Rho1*dRho0) / RhoFactor
c    2   + t(i,2)*2.0_cd/3.0_cd*Rho2/Rho0*(dRho2*Rho0-Rho2*dRho0) /RhoFactor
c    3   + t(i,3)*0.4_cd*Rho3/Rho0 * (dRho3*Rho0-Rho3*dRho0) / RhoFactor )
      endif
c
      if(rho_bar.eq.0.0_cd) rho_bar = 1.0e-20_cd
      rr0 = Rho_Bar / mp%Rho0_Bar(i)
      F = mp%A(i) * mp%Ec(i,i) * rr0 * dlog(rr0)

c     dFdR = A * Ec * ( 1 + dlog(rr0) ) / Rho0_Bar * dRho_dR
      dRho_dR = dRho0 * G_GAMMA + dRho_App
c     dRho_dR = dRho0 * dsqrt(1.0_cd + GAMMA) + dRho_App
      dFdR = mp%A(i) * mp%Ec(i,i) * (1.0_cd+dlog(rr0)) / mp%Rho0_Bar(i)
     &     * dRho_dR
c
c  Calculation of Phi and DPhi at r = aa * Dij
c
      Phi  = 2.0_cd * ( Eu - F ) / mp%Z(i)
      DPhi = 2.0_cd * ( dEdR - dFdR ) / mp%Z(i) / Dij
c
      return
      end

      subroutine Compute_EAM_Energy(Ipbc,ene_sum,NListMC,pkim,mp,ier)
c
c  Compute potential energy using MEAM for all atoms (NListMC=0) or
c  for selected atoms (NListMC.gt.0) for MC simulation
c
      use, intrinsic :: iso_c_binding
      use KIM_API_F03
      use meam_parameters
      implicit none

      integer(c_int) Ipbc(3)
      real(c_double) ene_sum
      integer(c_int) NListMC
      type(c_ptr) pkim
      type(meam_params) :: mp
      integer(c_int) ier

      common /abcI/ Natom,Idim,Iperform,Ishear
      common /atom/ IAtomType(200000)
      common /ener/ ene_pot(200000),ene_kin(200000)
      common /lat1/ pos(3,200000),vel(3,200000),acc(3,200000)
      common /lat2/ eps(3,3),eps_inv(3,3),BoxSize(3)
      common /nei4/ ListMC(1000)
      common /rho1/ RhoAtom(200000,4),Rho_0(200000),Rho_2(200000)
      common /rho2/ Rho1Sub(200000,3),Rho2Sub(200000,3,3)
     &             ,Rho3Sub(200000,3,3,3),Rho3Add(200000,3)
      integer(c_int) Natom, Idim, Iperform, Ishear
      integer(c_int) IAtomType
      real(c_double) ene_pot, ene_kin
      real(c_double) pos, vel, acc
      real(c_double) eps, eps_inv, BoxSize
      integer(c_int) ListMC
      real(c_double) RhoAtom, Rho_0, Rho_2
      real(c_double) Rho1Sub, Rho2Sub, Rho3Sub, Rho3Add

      real(c_double) Sij(3),Rij(3),Sik(3),Rik(3),Sjk(3),Rjk(3)
      real(c_double) Rho1S(3),Rho2S(3,3),Rho3S(3,3,3),Rho3A(3)
      integer(c_int) ListJ(100)
      integer(c_int) i, id, ip, ix, iy, iz, j, jp, k, kp, L
      integer(c_int) M, N_atoms, NListJ
      real(c_double) D, Dij, F, phi, Rfactor, Rho0, Rho2, Rho_Bar, rk,
     &   Rsqij, Rsqik, Rsqjk, scnres, screenij, Sikj, wt
      real(c_double), external :: cutoff
      real(c_double), external :: rhoi
      real(c_double), external :: screen

c
c  KIM variables
c
      integer(c_int) idum,atom_ret,numnei
      integer(c_int), pointer :: nei1atom(:); type(c_ptr) :: pnei1atom
      type(c_ptr) :: pRij_list

c
c  Reset to zero potential energies, forces and electron densities
c
      scnres = 1.0e-6_cd
      ene_sum = 0.0_cd
      ene_pot = 0.0_cd
      Rho_0   = 0.0_cd
      Rho_2   = 0.0_cd
      Rho1Sub = 0.0_cd
      Rho2Sub = 0.0_cd
      Rho3Sub = 0.0_cd
      Rho3Add = 0.0_cd
c
c  Loop over particles in the neighbor list for all or selected atoms
c
      if(NListMC.eq.0) then
         N_Atoms = Natom
      else
         N_Atoms = NListMC
      endif
c
      do id = 1, N_Atoms
         if(NListMC.eq.0) then
            i = id
         else
            i  = ListMC(id)
         endif
c
         ip = IAtomType(i)
         if(ip.lt.0) ip=-ip
         NListJ = 0
c
c Get neighbors of atom i via KIM API
c
         ier = kim_api_get_neigh(pkim,1,i,atom_ret,numnei,pnei1atom,
     &                             pRij_list)
         if (ier.ne.KIM_STATUS_OK) then ! some sort of problem, exit
            idum = kim_api_report_error(__LINE__, THIS_FILE_NAME,
     &                                  "kim_api_get_neigh", ier)
            ier = KIM_STATUS_FAIL
            return
         endif
         call c_f_pointer(pnei1atom, nei1atom, [numnei])

      do L = 1, numnei
         j = nei1atom(L)
         jp = IAtomType(j)
         if(jp.lt.0) jp=-jp
         Sij = pos(:,i) - pos(:,j)
         where ( abs(Sij) > 0.5_cd .and. Ipbc.eq.1 )
            Sij = Sij - sign(1.0_cd,Sij)
         end where
c
c  ============================================================================
c              Computation of pair potential and force due to that
c  ============================================================================
c
c        Rij = BoxSize*Sij
         call Rij_Real(Rij,Sij)
         Rsqij = dot_product(Rij,Rij)
         Dij = dsqrt(Rsqij)
         if( Dij < mp%Rcutoff ) then
c
c  ----------------------------------------------------------------------------
c        Compute radial cutoff function
c
            screenij = cutoff(mp%Rcutoff,mp%Dcutoff,Rsqij)
            if(screenij.le.scnres) goto 1000
c
c  ----------------------------------------------------------------------------
c  --  Here is the only part which is different from BJLEE version  -----------
c  ----------------------------------------------------------------------------
c        Compute screening function
c
            do M = 1,numnei
               k = nei1atom(M)
               if(k.ne.j) then
                  kp = IAtomType(k)
                  if(kp.lt.0) kp=-kp
c
c              evaluate Rik
c
                  Sik = pos(:,i) - pos(:,k)
                  where ( abs(Sik) > 0.5_cd .and. Ipbc.eq.1 )
                     Sik = Sik - sign(1.0_cd,Sik)
                  end where
c                 Rik = BoxSize*Sik
                  call Rij_Real(Rik,Sik)
                  Rsqik = dot_product(Rik,Rik) / Rsqij
                  if(Rsqik.gt.1.0_cd) goto 900
c
c              evaluate Rjk
c
                  Sjk = pos(:,j) - pos(:,k)
                  where ( abs(Sjk) > 0.5_cd .and. Ipbc.eq.1 )
                     Sjk = Sjk - sign(1.0_cd,Sjk)
                  end where
c                 Rjk = BoxSize*Sjk
                  call Rij_Real(Rjk,Sjk)
                  Rsqjk = dot_product(Rjk,Rjk) / Rsqij
                  if(Rsqjk.gt.1.0_cd) goto 900
c
c              compute Sikj: screening between i and j due to k
c
               Sikj = screen(Rsqik,Rsqjk,mp%Cmax(ip,kp,jp),
     &                       mp%Cmin(ip,kp,jp))
                  screenij = screenij * Sikj
                  if(screenij.le.scnres) goto 1000
               endif
  900          continue
            enddo
c
c  ----------------------------------------------------------------------------
c        "continuous" and discretized index in table
c
            rk = (Dij - mp%Rmin)*mp%OverDeltaR + 1.0_cd
            k = idint(rk)
            if (k < 1) k = 1
            wt = rk - k
c
c         linear interpolation of phi and 1/R * dphi/dR
c
            phi  = wt* mp%PhiTab(k+1,ip,jp)
     &           + (1.0_cd-wt)* mp%PhiTab(k,ip,jp)
c
c         accumulate energy
c
            ene_pot(i) = ene_pot(i) + 0.5_cd * phi * screenij
c           ene_pot(j) = ene_pot(j) + 0.5_cd * phi * screenij
c
c  ============================================================================
c              Computation of partial electron density from atom j
c  ============================================================================
c
            NListJ = NListJ + 1
            ListJ(NListJ) = j
            Rfactor = Dij / mp%Re(jp) - 1.0_cd
c
            do m = 1, 4
c              RhoAtom(i,m)=screenij*RhoZ(jp)*dexp(-beta(jp,m)*Rfactor)
               RhoAtom(j,m)=screenij*mp%RhoZ(jp)
     &                              *dexp(-mp%beta(jp,m)*Rfactor)
            enddo
c
            Rho_0(i) = Rho_0(i) + RhoAtom(j,4)
c           Rho_0(j) = Rho_0(j) + RhoAtom(i,4)
c
            Rho_2(i) = Rho_2(i) + RhoAtom(j,2)
c           Rho_2(j) = Rho_2(j) + RhoAtom(i,2)
c
            do ix = 1, 3
               Rho1Sub(i,ix) = Rho1Sub(i,ix) - Rij(ix)/Dij*RhoAtom(j,1)
c              Rho1Sub(j,ix) = Rho1Sub(j,ix) + Rij(ix)/Dij*RhoAtom(i,1)
c
c    Additive term after M.I. Baskes, Mater. Sci. Engineer. A261, 165 (1999)
c
               Rho3Add(i,ix) = Rho3Add(i,ix) - Rij(ix)/Dij*RhoAtom(j,3)
            do iy = 1, 3
               Rho2Sub(i,ix,iy) = Rho2Sub(i,ix,iy)
     &                          + Rij(ix)*Rij(iy)/Rsqij * RhoAtom(j,2)
c              Rho2Sub(j,ix,iy) = Rho2Sub(j,ix,iy)
c    &                          + Rij(ix)*Rij(iy)/Rsqij * RhoAtom(i,2)
            do iz = 1, 3
               Rho3Sub(i,ix,iy,iz) = Rho3Sub(i,ix,iy,iz)
     &            - Rij(ix)*Rij(iy)*Rij(iz)/Dij/Rsqij * RhoAtom(j,3)
c              Rho3Sub(j,ix,iy,iz) = Rho3Sub(j,ix,iy,iz)
c    &         + Rij(ix)*Rij(iy)*Rij(iz)/Dij/Rsqij * RhoAtom(i,3)
            enddo
            enddo
            enddo
 1000       continue
         endif
      enddo
c
c  ============================================================================
c      Computation of Rho_Bar, embedding potential and force due to that
c      in order to calculate local stress, force from atom j will be computed
c  ============================================================================
c
      Rho0 = Rho_0(i)
      Rho2 = Rho_2(i)
      Rho1S = Rho1Sub(i,:)
      Rho2S = Rho2Sub(i,:,:)
      Rho3S = Rho3Sub(i,:,:,:)
      Rho3A = Rho3Add(i,:)
c
      Rho_Bar = rhoi(mp%t(ip,:),rho0,rho2,rho1s,rho2s,rho3s,rho3a)
      D = Rho_Bar / mp%Rho0_Bar(ip)
      if(D.lt.1.e-20_cd) then
         D = 1.e-20_cd
      endif
      F = mp%A(ip) * mp%Ec(ip,ip) * D * dlog(D)
c      print *,"F",A(ip),Ec(ip,ip),D,dlog(D)
c
c     accumulate energy
c
      ene_pot(i) = ene_pot(i) + F
      ene_sum = ene_sum + ene_pot(i)

!      print *,"energy:",i,ene_pot(i),ene_sum

      enddo
      return
      end
c
c
      Subroutine Rij_Real(Rij,Sij)
c
c  Calculation of Real Distance between two atoms is made here from eps(3,3),
c  Sij(3) and BoxSize(3), considering shear,   2004. 1. 12.  B.-J. Lee
c
      use, intrinsic :: iso_c_binding
      implicit none

      real(c_double) Rij(3)
      real(c_double) Sij(3)

      common /lat2/ eps(3,3),eps_inv(3,3),BoxSize(3)
      real(c_double) eps, eps_inv, BoxSize
c
      Rij(1) = eps(1,1)*BoxSize(1)*Sij(1)
     &       + eps(1,2)*BoxSize(2)*Sij(2) + eps(1,3)*BoxSize(3)*Sij(3)
      Rij(2) = eps(2,1)*BoxSize(1)*Sij(1)
     &       + eps(2,2)*BoxSize(2)*Sij(2) + eps(2,3)*BoxSize(3)*Sij(3)
      Rij(3) = eps(3,1)*BoxSize(1)*Sij(1)
     &       + eps(3,2)*BoxSize(2)*Sij(2) + eps(3,3)*BoxSize(3)*Sij(3)
c
      return
      end

      Subroutine Rij_Inv(Rij,Sij)
c
c  Calculation of Real Distance between two atoms is made here from eps(3,3),
c  Sij(3) and BoxSize(3), considering shear,   2004. 1. 12.  B.-J. Lee
c
      use, intrinsic :: iso_c_binding
      implicit none

      real(c_double) Rij(3)
      real(c_double) Sij(3)

      common /lat2/ eps(3,3),eps_inv(3,3),BoxSize(3)
      real(c_double) eps, eps_inv, BoxSize
c
      Sij(1) = eps_inv(1,1)*Rij(1)
     &       + eps_inv(1,2)*Rij(2) + eps_inv(1,3)*Rij(3)
      Sij(2) = eps_inv(2,1)*Rij(1)
     &       + eps_inv(2,2)*Rij(2) + eps_inv(2,3)*Rij(3)
      Sij(3) = eps_inv(3,1)*Rij(1)
     &       + eps_inv(3,2)*Rij(2) + eps_inv(3,3)*Rij(3)
c
      return
      end

      subroutine Compute_EAM_Forces(Ipbc,f_stress,pkim,mp,ier)
c
c  Compute potential energy and forces on atoms, using MEAM potential.
c
      use, intrinsic :: iso_c_binding
      use KIM_API_F03
      use meam_parameters
      implicit none

      integer(c_int) Ipbc(3)
      real(c_double) f_stress(3,3)
      type(c_ptr) pkim
      type(meam_params) :: mp
      integer(c_int) ier

      common /abcI/ Natom,Idim,Iperform,Ishear
      common /atom/ IAtomType(200000)
      common /dbtt/ N_fix,N_rotate,N_dbtt
      common /ener/ ene_pot(200000),ene_kin(200000)
      common /lat1/ pos(3,200000),vel(3,200000),acc(3,200000)
      common /lat2/ eps(3,3),eps_inv(3,3),BoxSize(3)
      common /rho1/ RhoAtom(200000,4),Rho_0(200000),Rho_2(200000)
      common /rho2/ Rho1Sub(200000,3),Rho2Sub(200000,3,3)
     &             ,Rho3Sub(200000,3,3,3),Rho3Add(200000,3)
      integer(c_int) Natom, Idim, Iperform, Ishear
      integer(c_int) IAtomType
      integer(c_int) N_fix, N_rotate, N_dbtt
      real(c_double) ene_pot, ene_kin
      real(c_double) pos, vel, acc
      real(c_double) eps, eps_inv, BoxSize
      real(c_double) RhoAtom, Rho_0, Rho_2
      real(c_double) Rho1Sub, Rho2Sub, Rho3Sub, Rho3Add

      real(c_double) Sij(3),Rij(3),force(3,200000)
      real(c_double) Sik(3),Rik(3),Sjk(3),Rjk(3)
      real(c_double) dscreenij(3),DRij(3),scni(6),scnj(6)
      real(c_double) DRik(3,100),DRjk(3,100),DRsqik(100),DRsqjk(100)
      real(c_double) DDij(6),DDRij(6,3),DRfactor(6),DF(6),FF(3)
      real(c_double) RRij(3,100)
      real(c_double) DiRhoAtom(6,100,4),DiRho_0(6,100),DiRho_2(6,100)
      real(c_double) DiRho1Sub(6,100,3),DiRho2Sub(6,100,3,3)
      real(c_double) DiRho3Sub(6,100,3,3,3),DiRho3Add(6,100,3)
      real(c_double) DjRhoAtom(6,100,4),DjRho_0(6,100),DjRho_2(6,100)
      real(c_double) DjRho1Sub(6,100,3),DjRho2Sub(6,100,3,3)
      real(c_double) DjRho3Sub(6,100,3,3,3),DjRho3Add(6,100,3)
      real(c_double) Rho1S(3),Rho2S(3,3),Rho3S(3,3,3),Rho3A(3)
      real(c_double) ScreenK(6,100,100),RRik(3,100,100)
      real(c_double) DkRho_0(6,100,100),DkRho_2(6,100,100)
      real(c_double) DkRho1Sub(6,100,100,3),DkRho2Sub(6,100,100,3,3)
      real(c_double) DkRho3Sub(6,100,100,3,3,3),DkRhoAtom(6,100,100,4)
      real(c_double) DkRho3Add(6,100,100,3)
      integer(c_int) ListJ(100),NListM(100),ListK(100,100)
      integer(c_int) KAtomType(100)
      integer(c_int) i, id, ikj, ip, ix, iy, iz, j, jp, k, kp, L, M,
     &    NListJ, NListK, NLM
      real(c_double) D, DDsqij, denergy, Dij, phi, dphi, DRsqij,
     &   dscnmi, dscnmj, dscnpi, dscnpj, F, h, Rfactor, Rho0, Rho2,
     &   Rho_Bar, rk, Rsqij, Rsqik, Rsqjk, scnres, screenij,
     &   Sikj, Sikj_M, Sikj_P, wt
      real(c_double), external :: cutoff
      real(c_double), external :: screen
      real(c_double), external :: rhoi

c
c KIM variables
c
      integer(c_int) idum,atom_ret,numnei
      integer(c_int), pointer :: nei1atom(:); type(c_ptr) :: pnei1atom
      type(c_ptr) :: pRij_list

c
c  Initial values
c
      scnres = 1.0e-6_cd
      h      = 1.0e-5_cd
      f_stress = 0.0_cd
c
c  Reset to zero potential energies, forces and electron densities
c
      ene_pot = 0.0_cd
      force   = 0.0_cd
      Rho_0   = 0.0_cd
      Rho_2   = 0.0_cd
      Rho1Sub = 0.0_cd
      Rho2Sub = 0.0_cd
      Rho3Sub = 0.0_cd
      Rho3Add = 0.0_cd
c
c  Loop over particles in the neighbor list
c
      do i = 1,Natom
         ip = IAtomType(i)
         if(ip.lt.0) ip=-ip
         NListJ = 0

c
c Get neighbors of atom i via KIM API
c
         ier = kim_api_get_neigh(pkim,1,i,atom_ret,numnei,pnei1atom,
     &                             pRij_list)
         if (ier.ne.KIM_STATUS_OK) then ! some sort of problem, exit
            idum = kim_api_report_error(__LINE__, THIS_FILE_NAME,
     &                                  "kim_api_get_neigh", ier)
            ier = KIM_STATUS_FAIL
            return
         endif
         call c_f_pointer(pnei1atom, nei1atom, [numnei])
c
      do L = 1,numnei
         j = nei1atom(L)
         jp = IAtomType(j)
         if(jp.lt.0) jp=-jp
         Sij = pos(:,i) - pos(:,j)
         where ( abs(Sij) > 0.5_cd .and. Ipbc.eq.1 )
            Sij = Sij - sign(1.0_cd,Sij)
         end where
c
c  ============================================================================
c              Computation of pair potential and force due to that
c  ============================================================================
c
c        Rij = BoxSize*Sij
         call Rij_Real(Rij,Sij)
         Rsqij = dot_product(Rij,Rij)
         Dij = dsqrt(Rsqij)
         if( Dij < mp%Rcutoff ) then
c
c  ----------------------------------------------------------------------------
c        Compute radial cutoff function
c
            screenij = cutoff(mp%Rcutoff,mp%Dcutoff,Rsqij)
            if(screenij.le.scnres) goto 1000
c
c  ----------------------------------------------------------------------------
c        Compute screening function
c
            NListK = 0
            do M = 1,numnei
               k = nei1atom(M)
               if(k.ne.j) then
                  kp = IAtomType(k)
                  if(kp.lt.0) kp=-kp
c
c              evaluate Rik
c
                  Sik = pos(:,i) - pos(:,k)
                  where ( abs(Sik) > 0.5_cd .and. Ipbc.eq.1 )
                     Sik = Sik - sign(1.0_cd,Sik)
                  end where
c                 Rik = BoxSize*Sik
                  call Rij_Real(Rik,Sik)
                  Rsqik = dot_product(Rik,Rik) / Rsqij
                  if(Rsqik.gt.1.0_cd) goto 900
c
c              evaluate Rjk
c
                  Sjk = pos(:,j) - pos(:,k)
                  where ( abs(Sjk) > 0.5_cd .and. Ipbc.eq.1 )
                     Sjk = Sjk - sign(1.0_cd,Sjk)
                  end where
c                 Rjk = BoxSize*Sjk
                  call Rij_Real(Rjk,Sjk)
                  Rsqjk = dot_product(Rjk,Rjk) / Rsqij
                  if(Rsqjk.gt.1.0_cd) goto 900
c
c              compute Sikj: screening between i and j due to k
c
               Sikj = screen(Rsqik,Rsqjk,mp%Cmax(ip,kp,jp),
     &                       mp%Cmin(ip,kp,jp))
                  if(Sikj.ge.scnres.and. Sikj.le.1.0_cd-scnres) then
                     NListK = NListK + 1
                     DRik(:,NListK) = Rik
                     DRjk(:,NListK) = Rjk
                     DRsqik(NListK) = Rsqik * Rsqij
                     DRsqjk(NListK) = Rsqjk * Rsqij
                     KAtomType(NListK) = kp
                  endif
                  screenij = screenij * Sikj
                  if(screenij.le.scnres) goto 1000
               endif
  900          continue
            enddo
c
c  ----------------------------------------------------------------------------
c        Compute derivatives of screening function between atom i and atom j
c
            if(screenij.gt.1.0_cd-scnres .or. screenij.lt.scnres) then
               dscreenij = 0.0_cd
               scni = screenij
               scnj = screenij
            else
               DRij = Rij
               do id = 1, 3
c
c              Compute radial cutoff function at xi = xi + h or at xj = xj - h
c
                  DRij(id) = Rij(id) + h
                  DRsqij   = dot_product(DRij,DRij)
                  dscnpi  = cutoff(mp%Rcutoff,mp%Dcutoff,DRsqij)
                  dscnpj  = dscnpi
c
c              Compute screen function at xi = xi + h and at xj = xj - h
c              We will consider only the k atoms that don't give 1 or 0
c
                  do M = 1, NListK
                     kp = KAtomType(M)
                     Rik     = DRik(:,M)
                     Rik(id) = Rik(id) + h
                     Rsqik   = dot_product(Rik,Rik) / DRsqij
                     Rsqjk   = DRsqjk(M) / DRsqij
                  Sikj=screen(Rsqik,Rsqjk,mp%Cmax(ip,kp,jp),
     &                        mp%Cmin(ip,kp,jp))
                     dscnpi = dscnpi * Sikj
c
                     Rjk     = DRjk(:,M)
                     Rjk(id) = Rjk(id) - h
                     Rsqjk   = dot_product(Rjk,Rjk) / DRsqij
                     Rsqik   = DRsqik(M) / DRsqij
                  Sikj=screen(Rsqik,Rsqjk,mp%Cmax(ip,kp,jp),
     &                        mp%Cmin(ip,kp,jp))
                     dscnpj = dscnpj * Sikj
                  enddo
c
c              Compute radial cutoff function at xi = xi - h or at xj = xj + h
c
                  DRij(id) = Rij(id) - h
                  DRsqij   = dot_product(DRij,DRij)
                  dscnmi  = cutoff(mp%Rcutoff,mp%Dcutoff,DRsqij)
                  dscnmj  = dscnmi
c
c              Compute screen function at xi = xi - h and at xj = xj + h
c
                  do M = 1, NListK
                     kp = KAtomType(M)
                     Rik     = DRik(:,M)
                     Rik(id) = Rik(id) - h
                     Rsqik   = dot_product(Rik,Rik) / DRsqij
                     Rsqjk   = DRsqjk(M) / DRsqij
                  Sikj=screen(Rsqik,Rsqjk,mp%Cmax(ip,kp,jp),
     &                        mp%Cmin(ip,kp,jp))
                     dscnmi = dscnmi * Sikj
c
                     Rjk     = DRjk(:,M)
                     Rjk(id) = Rjk(id) + h
                     Rsqjk   = dot_product(Rjk,Rjk) / DRsqij
                     Rsqik   = DRsqik(M) / DRsqij
                  Sikj=screen(Rsqik,Rsqjk,mp%Cmax(ip,kp,jp),
     &                        mp%Cmin(ip,kp,jp))
                     dscnmj = dscnmj * Sikj
                  enddo
c
                  DRij(id) = Rij(id)
c
c              Compute derivatives of screening function
c
                  dscreenij(id) = 0.5_cd * (dscnpi - dscnmi) / h
                  scni(2*id-1) = dscnpi
                  scni(2*id)   = dscnmi
                  scnj(2*id-1) = dscnpj
                  scnj(2*id)   = dscnmj
               enddo
            endif
c
c  ----------------------------------------------------------------------------
c
c        "continuous" and discretized index in table
c
            rk = (Dij - mp%Rmin)*mp%OverDeltaR + 1.0_cd
            k = idint(rk)
            if (k < 1) k = 1
            wt = rk - k
c
c         linear interpolation of phi and 1/R * dphi/dR
c
            phi  = wt* mp%PhiTab(k+1,ip,jp)
     &           + (1.0_cd-wt)* mp%PhiTab(k,ip,jp)
            dphi = wt*mp%DPhiTab(k+1,ip,jp)
     &           + (1.0_cd-wt)*mp%DPhiTab(k,ip,jp)
c
c         accumulate energy
c
            ene_pot(i) = ene_pot(i) + 0.5_cd * phi * screenij
c           if(i.eq.1) write(*,*) j, dij, screenij
c           ene_pot(j) = ene_pot(j) + 0.5_cd * phi * screenij
c
c         accumulate force and local stress for i atom
c
            FF = - dscreenij * phi - screenij * Rij * dphi
            force(:,i) = force(:,i) + FF
c
            f_stress(1,1) = f_stress(1,1) + 0.5_cd * FF(1) * Rij(1)
            f_stress(1,2) = f_stress(1,2) + 0.5_cd * FF(1) * Rij(2)
            f_stress(1,3) = f_stress(1,3) + 0.5_cd * FF(1) * Rij(3)
            f_stress(2,2) = f_stress(2,2) + 0.5_cd * FF(2) * Rij(2)
            f_stress(2,3) = f_stress(2,3) + 0.5_cd * FF(2) * Rij(3)
            f_stress(3,3) = f_stress(3,3) + 0.5_cd * FF(3) * Rij(3)
c
c  ============================================================================
c   Compute energy changes and thus created forces due to movement of k atoms
c  ============================================================================
c
            NListJ = NListJ + 1
            ListJ(NListJ) = j
c
            NListM(NListJ) = 0
            NLM            = 0
            do M = 1,numnei
               k = nei1atom(M)
               if(k.eq.j) goto 950
               kp = IAtomType(k)
               if(kp.lt.0) kp=-kp
c
c           evaluate Rik
c
               Sik = pos(:,i) - pos(:,k)
               where ( abs(Sik) > 0.5_cd .and. Ipbc.eq.1 )
                  Sik = Sik - sign(1.0_cd,Sik)
               end where
c              Rik = BoxSize*Sik
               call Rij_Real(Rik,Sik)
               Rsqik = dot_product(Rik,Rik) / Rsqij
               if(Rsqik.ge.1.0_cd) goto 950
c
c           evaluate Rjk
c
               Sjk = pos(:,j) - pos(:,k)
               where ( abs(Sjk) > 0.5_cd .and. Ipbc.eq.1 )
                  Sjk = Sjk - sign(1.0_cd,Sjk)
               end where
c              Rjk = BoxSize*Sjk
               call Rij_Real(Rjk,Sjk)
               Rsqjk = dot_product(Rjk,Rjk) / Rsqij
               if(Rsqjk.ge.1.0_cd) goto 950
c
c           compute Sikj: screening between i and j due to k
c
               Sikj = screen(Rsqik,Rsqjk,mp%Cmax(ip,kp,jp),
     &                       mp%Cmin(ip,kp,jp))
               if(Sikj.gt.1.0_cd-scnres .or. Sikj.lt.scnres) goto 950
               NListM(NListJ) = NListM(NListJ) + 1
               NLM            = NListM(NListJ)
               ListK(NListJ,NLM) = k
c
               do id = 1, 3
c
c              Compute screen function at xk = xk + h
c
                  Rik(id) = Rik(id) - h
                  Rsqik   = dot_product(Rik,Rik) / Rsqij
                  Rjk(id) = Rjk(id) - h
                  Rsqjk   = dot_product(Rjk,Rjk) / Rsqij
                  Sikj_P =
     &            screen(Rsqik,Rsqjk,mp%Cmax(ip,kp,jp),
     &                   mp%Cmin(ip,kp,jp))/Sikj
c
                  Rik(id) = Rik(id) + 2.0_cd * h
                  Rsqik   = dot_product(Rik,Rik) / Rsqij
                  Rjk(id) = Rjk(id) + 2.0_cd * h
                  Rsqjk   = dot_product(Rjk,Rjk) / Rsqij
                  Sikj_M =
     &            screen(Rsqik,Rsqjk,mp%Cmax(ip,kp,jp),
     &                   mp%Cmin(ip,kp,jp))/Sikj
c
                  Rik(id) = Rik(id) - h
                  Rjk(id) = Rjk(id) - h
c
c              Compute energy change and force due to that
c
                  denergy = 0.5_cd * phi * screenij * (Sikj_P - Sikj_M)
                  FF(id) = - 0.5_cd * denergy / h

                  ScreenK(2*id-1,NListJ,NListM(NListJ)) = Sikj_P
                  ScreenK(2*id,NListJ,NListM(NListJ))   = Sikj_M
               enddo
               force(:,k) = force(:,k) + FF
c
               f_stress(1,1) = f_stress(1,1) - FF(1) * Rik(1)
               f_stress(1,2) = f_stress(1,2) - FF(1) * Rik(2)
               f_stress(1,3) = f_stress(1,3) - FF(1) * Rik(3)
               f_stress(2,2) = f_stress(2,2) - FF(2) * Rik(2)
               f_stress(2,3) = f_stress(2,3) - FF(2) * Rik(3)
               f_stress(3,3) = f_stress(3,3) - FF(3) * Rik(3)
               RRik(:,NListJ,NLM) = Rik
  950          continue
            enddo
c
c  ============================================================================
c              Computation of partial electron density from atom j
c  ============================================================================
c
            RRij(:,NListJ) = Rij
            Rfactor = Dij / mp%Re(jp) - 1.0_cd
c
            do id = 1, 6
               DDRij(id,:) = Rij
               if(id.eq.1) DDRij(id,1) = Rij(1) + h
               if(id.eq.2) DDRij(id,1) = Rij(1) - h
               if(id.eq.3) DDRij(id,2) = Rij(2) + h
               if(id.eq.4) DDRij(id,2) = Rij(2) - h
               if(id.eq.5) DDRij(id,3) = Rij(3) + h
               if(id.eq.6) DDRij(id,3) = Rij(3) - h
               DDsqij = DDRij(id,1)*DDRij(id,1)
     &                + DDRij(id,2)*DDRij(id,2)
     &                + DDRij(id,3)*DDRij(id,3)
               DDij(id) = dsqrt(DDsqij)
               DRfactor(id) = DDij(id) / mp%Re(jp) - 1.0_cd
            enddo
c
            do m = 1, 4
c              RhoAtom(i,m)=screenij*RhoZ(jp)*dexp(-beta(jp,m)*Rfactor)
               RhoAtom(j,m)=screenij*mp%RhoZ(jp)
     &                              *dexp(-mp%beta(jp,m)*Rfactor)
               do id = 1, 6
                  DiRhoAtom(id,NListJ,m) = scni(id) * mp%RhoZ(jp)
     &                               * dexp(-mp%beta(jp,m)
     &                               * DRfactor(id))
                  DjRhoAtom(id,NListJ,m) = scnj(id) * mp%RhoZ(jp)
     &                               * dexp(-mp%beta(jp,m)
     &                               * DRfactor(id))
                  do ikj = 1, NLM
                     DkRhoAtom(id,NListJ,ikj,m)
     &                 = ScreenK(id,NListJ,ikj) * RhoAtom(j,m)
                  enddo
               enddo
            enddo
c
            Rho_0(i) = Rho_0(i) + RhoAtom(j,4)
            do id = 1, 6
               DiRho_0(id,NListJ) = DiRhoAtom(id,NListJ,4)
     &                            - RhoAtom(j,4)
               DjRho_0(id,NListJ) = DjRhoAtom(id,NListJ,4)
     &                            - RhoAtom(j,4)
               do ikj = 1, NLM
                  DkRho_0(id,NListJ,ikj) = DkRhoAtom(id,NListJ,ikj,4)
     &                                   - RhoAtom(j,4)
               enddo
            enddo
c           Rho_0(j) = Rho_0(j) + RhoAtom(i,4)
c
            Rho_2(i) = Rho_2(i) + RhoAtom(j,2)
            do id = 1, 6
               DiRho_2(id,NListJ) = DiRhoAtom(id,NListJ,2)
     &                            - RhoAtom(j,2)
               DjRho_2(id,NListJ) = DjRhoAtom(id,NListJ,2)
     &                            - RhoAtom(j,2)
               do ikj = 1, NLM
                  DkRho_2(id,NListJ,ikj) = DkRhoAtom(id,NListJ,ikj,2)
     &                                   - RhoAtom(j,2)
               enddo
            enddo
c           Rho_2(j) = Rho_2(j) + RhoAtom(i,2)
c
            do ix = 1, 3
               Rho1Sub(i,ix) = Rho1Sub(i,ix) - Rij(ix)/Dij*RhoAtom(j,1)
               do id = 1, 6
                  DiRho1Sub(id,NListJ,ix) =
     &              - DDRij(id,ix)/DDij(id)*DiRhoAtom(id,NListJ,1)
     &              + Rij(ix)/Dij*RhoAtom(j,1)
                  DjRho1Sub(id,NListJ,ix) =
     &              - DDRij(id,ix)/DDij(id)*DjRhoAtom(id,NListJ,1)
     &              + Rij(ix)/Dij*RhoAtom(j,1)
                  do ikj = 1, NLM
                     DkRho1Sub(id,NListJ,ikj,ix) =
     &                 - Rij(ix)/Dij*DkRhoAtom(id,NListJ,ikj,1)
     &                 + Rij(ix)/Dij*RhoAtom(j,1)
                  enddo
               enddo
c
c    Additive term after M.I. Baskes, Mater. Sci. Engineer. A261, 165 (1999)
c
               Rho3Add(i,ix) = Rho3Add(i,ix) - Rij(ix)/Dij*RhoAtom(j,3)
               do id = 1, 6
                  DiRho3Add(id,NListJ,ix) =
     &              - DDRij(id,ix)/DDij(id)*DiRhoAtom(id,NListJ,3)
     &              + Rij(ix)/Dij*RhoAtom(j,3)
                  DjRho3Add(id,NListJ,ix) =
     &              - DDRij(id,ix)/DDij(id)*DjRhoAtom(id,NListJ,3)
     &              + Rij(ix)/Dij*RhoAtom(j,3)
                  do ikj = 1, NLM
                     DkRho3Add(id,NListJ,ikj,ix) =
     &                 - Rij(ix)/Dij*DkRhoAtom(id,NListJ,ikj,3)
     &                 + Rij(ix)/Dij*RhoAtom(j,3)
                  enddo
               enddo
c              Rho1Sub(j,ix) = Rho1Sub(j,ix) + Rij(ix)/Dij*RhoAtom(i,1)
            do iy = 1, 3
               Rho2Sub(i,ix,iy) = Rho2Sub(i,ix,iy)
     &                          + Rij(ix)*Rij(iy)/Rsqij * RhoAtom(j,2)
               do id = 1, 6
                  DiRho2Sub(id,NListJ,ix,iy) =
     &              + DDRij(id,ix)*DDRij(id,iy)/DDij(id)/DDij(id)
     &              * DiRhoAtom(id,NListJ,2)
     &              - Rij(ix)*Rij(iy)/Rsqij * RhoAtom(j,2)
                  DjRho2Sub(id,NListJ,ix,iy) =
     &              + DDRij(id,ix)*DDRij(id,iy)/DDij(id)/DDij(id)
     &              * DjRhoAtom(id,NListJ,2)
     &              - Rij(ix)*Rij(iy)/Rsqij * RhoAtom(j,2)
                  do ikj = 1, NLM
                     DkRho2Sub(id,NListJ,ikj,ix,iy) =
     &               + Rij(ix)*Rij(iy)/Rsqij*DkRhoAtom(id,NListJ,ikj,2)
     &               - Rij(ix)*Rij(iy)/Rsqij * RhoAtom(j,2)
                  enddo
               enddo
c              Rho2Sub(j,ix,iy) = Rho2Sub(j,ix,iy)
c    &                          + Rij(ix)*Rij(iy)/Rsqij * RhoAtom(i,2)
            do iz = 1, 3
               Rho3Sub(i,ix,iy,iz) = Rho3Sub(i,ix,iy,iz)
     &            - Rij(ix)*Rij(iy)*Rij(iz)/Dij/Rsqij * RhoAtom(j,3)
               do id = 1, 6
                  DiRho3Sub(id,NListJ,ix,iy,iz) =
     &              - DDRij(id,ix)*DDRij(id,iy)*DDRij(id,iz)
     &              / DDij(id) / DDij(id) / DDij(id)
     &              * DiRhoAtom(id,NListJ,3)
     &              + Rij(ix)*Rij(iy)*Rij(iz)/Dij/Rsqij * RhoAtom(j,3)
                  DjRho3Sub(id,NListJ,ix,iy,iz) =
     &              - DDRij(id,ix)*DDRij(id,iy)*DDRij(id,iz)
     &              / DDij(id) / DDij(id) / DDij(id)
     &              * DjRhoAtom(id,NListJ,3)
     &              + Rij(ix)*Rij(iy)*Rij(iz)/Dij/Rsqij * RhoAtom(j,3)
                  do ikj = 1, NLM
                     DkRho3Sub(id,NListJ,ikj,ix,iy,iz) =
     &               - Rij(ix)*Rij(iy)*Rij(iz)/Dij/Rsqij
     &               * DkRhoAtom(id,NListJ,ikj,3)
     &               + Rij(ix)*Rij(iy)*Rij(iz)/Dij/Rsqij * RhoAtom(j,3)
                  enddo
               enddo
c              Rho3Sub(j,ix,iy,iz) = Rho3Sub(j,ix,iy,iz)
c    &         + Rij(ix)*Rij(iy)*Rij(iz)/Dij/Rsqij * RhoAtom(i,3)
            enddo
            enddo
            enddo
 1000       continue
         endif
      enddo
c
c  ============================================================================
c      Computation of Rho_Bar, embedding potential and force due to that
c      in order to calculate local stress, force from atom j will be computed
c  ============================================================================
c
      Rho0 = Rho_0(i)
      Rho2 = Rho_2(i)
      Rho1S = Rho1Sub(i,:)
      Rho2S = Rho2Sub(i,:,:)
      Rho3S = Rho3Sub(i,:,:,:)
      Rho3A = Rho3Add(i,:)
c
      Rho_Bar = rhoi(mp%t(ip,:),rho0,rho2,rho1s,rho2s,rho3s,rho3a)
      D = Rho_Bar / mp%Rho0_Bar(ip)
      if(D.lt.1.e-20_cd) then
         D = 1.e-20_cd
      endif
      F = mp%A(ip) * mp%Ec(ip,ip) * D * dlog(D)
c
c     accumulate energy
c
      ene_pot(i) = ene_pot(i) + F
c
      do L = 1, NListJ
         j = ListJ(L)
c
c     compute force on atom i from neighbor atom j
c
         do id = 1, 6
            Rho0 = Rho_0(i) + DiRho_0(id,L)
            Rho2 = Rho_2(i) + DiRho_2(id,L)
            Rho1S = Rho1Sub(i,:) + DiRho1Sub(id,L,:)
            Rho2S = Rho2Sub(i,:,:) + DiRho2Sub(id,L,:,:)
            Rho3S = Rho3Sub(i,:,:,:) + DiRho3Sub(id,L,:,:,:)
            Rho3A = Rho3Add(i,:) + DiRho3Add(id,L,:)
            Rho_Bar = rhoi(mp%t(ip,:),rho0,rho2,rho1s,rho2s,rho3s,rho3a)
            D = Rho_Bar / mp%Rho0_Bar(ip)
            DF(id) = mp%A(ip) * mp%Ec(ip,ip) * D * dlog(D)
         enddo
         FF(1) = 0.5_cd * (DF(2) - DF(1)) / h
         FF(2) = 0.5_cd * (DF(4) - DF(3)) / h
         FF(3) = 0.5_cd * (DF(6) - DF(5)) / h
c
c     accumulate force and local stress
c
         force(:,i) = force(:,i) + FF
c
         f_stress(1,1) = f_stress(1,1) + 0.5_cd * FF(1) * RRij(1,L)
         f_stress(1,2) = f_stress(1,2) + 0.5_cd * FF(1) * RRij(2,L)
         f_stress(1,3) = f_stress(1,3) + 0.5_cd * FF(1) * RRij(3,L)
         f_stress(2,2) = f_stress(2,2) + 0.5_cd * FF(2) * RRij(2,L)
         f_stress(2,3) = f_stress(2,3) + 0.5_cd * FF(2) * RRij(3,L)
         f_stress(3,3) = f_stress(3,3) + 0.5_cd * FF(3) * RRij(3,L)
c
c     compute force on neighbor atom j from atom i
c
         do id = 1, 6
            Rho0 = Rho_0(i) + DjRho_0(id,L)
            Rho2 = Rho_2(i) + DjRho_2(id,L)
            Rho1S = Rho1Sub(i,:) + DjRho1Sub(id,L,:)
            Rho2S = Rho2Sub(i,:,:) + DjRho2Sub(id,L,:,:)
            Rho3S = Rho3Sub(i,:,:,:) + DjRho3Sub(id,L,:,:,:)
            Rho3A = Rho3Add(i,:) + DjRho3Add(id,L,:)
            Rho_Bar = rhoi(mp%t(ip,:),rho0,rho2,rho1s,rho2s,rho3s,rho3a)
            D = Rho_Bar / mp%Rho0_Bar(ip)
            DF(id) = mp%A(ip) * mp%Ec(ip,ip) * D * dlog(D)
         enddo
         FF(1) = 0.5_cd * (DF(2) - DF(1)) / h
         FF(2) = 0.5_cd * (DF(4) - DF(3)) / h
         FF(3) = 0.5_cd * (DF(6) - DF(5)) / h
c
c     accumulate force and local stress
c
         force(:,j) = force(:,j) - FF
c
         f_stress(1,1) = f_stress(1,1) + 0.5_cd * FF(1) * RRij(1,L)
         f_stress(1,2) = f_stress(1,2) + 0.5_cd * FF(1) * RRij(2,L)
         f_stress(1,3) = f_stress(1,3) + 0.5_cd * FF(1) * RRij(3,L)
         f_stress(2,2) = f_stress(2,2) + 0.5_cd * FF(2) * RRij(2,L)
         f_stress(2,3) = f_stress(2,3) + 0.5_cd * FF(2) * RRij(3,L)
         f_stress(3,3) = f_stress(3,3) + 0.5_cd * FF(3) * RRij(3,L)
c
c     Compute forces on screening atom k due to movement of the k atom
c
         do ikj = 1, NListM(L)
            k = ListK(L,ikj)
            do id = 1, 6
               Rho0 = Rho_0(i) + DkRho_0(id,L,ikj)
               Rho2 = Rho_2(i) + DkRho_2(id,L,ikj)
               Rho1S = Rho1Sub(i,:) + DkRho1Sub(id,L,ikj,:)
               Rho2S = Rho2Sub(i,:,:) + DkRho2Sub(id,L,ikj,:,:)
               Rho3S = Rho3Sub(i,:,:,:) + DkRho3Sub(id,L,ikj,:,:,:)
               Rho3A = Rho3Add(i,:) + DkRho3Add(id,L,ikj,:)
               Rho_Bar = rhoi(mp%t(ip,:),rho0,rho2,rho1s,rho2s,rho3s,
     &                        rho3a)
               D = Rho_Bar / mp%Rho0_Bar(ip)
               DF(id) = mp%A(ip) * mp%Ec(ip,ip) * D * dlog(D)
            enddo
            FF(1) = 0.5_cd * (DF(2) - DF(1)) / h
            FF(2) = 0.5_cd * (DF(4) - DF(3)) / h
            FF(3) = 0.5_cd * (DF(6) - DF(5)) / h
c
c        accumulate force and local stress
c
            force(:,k) = force(:,k) + FF
c
            f_stress(1,1) = f_stress(1,1) - FF(1)*RRik(1,L,ikj)
            f_stress(1,2) = f_stress(1,2) - FF(1)*RRik(2,L,ikj)
            f_stress(1,3) = f_stress(1,3) - FF(1)*RRik(3,L,ikj)
            f_stress(2,2) = f_stress(2,2) - FF(2)*RRik(2,L,ikj)
            f_stress(2,3) = f_stress(2,3) - FF(2)*RRik(3,L,ikj)
            f_stress(3,3) = f_stress(3,3) - FF(3)*RRik(3,L,ikj)
         enddo
      enddo
      enddo
c
c   compute acceleration on atom i
c
      do i = 1,Natom
         ip = IAtomType(i)
         if(ip.lt.0) then
            force(:,i) = 0.0_cd
            ip=-ip
         endif

         if(Ishear .eq. 0) then
            acc(:,i) = force(:,i) / mp%Amass(ip) / BoxSize
         else
            call Rij_Inv(force(:,i),Sij)
            acc(:,i) = Sij / mp%Amass(ip) / BoxSize
         endif
      enddo
c
      return
      end

      real(c_double) function cutoff(Rcutoff,Dcutoff,Rsq)
      use, intrinsic :: iso_c_binding
      implicit none

      integer(c_int), parameter :: cd = c_double
      real(c_double) Rcutoff
      real(c_double) Dcutoff
      real(c_double) Rsq

      real(c_double) RR, fac
c
      RR = dsqrt(Rsq)
      fac = ( Rcutoff - RR ) / Dcutoff
      if(fac.ge.1.0_cd) then
         cutoff = 1.0_cd
      else
         cutoff = ( 1.0_cd - (1.0_cd - fac)**4 )**2
      endif
      return
      end

      real(c_double) function screen(Rsqik,Rsqjk,Cmax,Cmin)
      use, intrinsic :: iso_c_binding
      implicit none

      integer(c_int), parameter :: cd = c_double
      real(c_double) Rsqik
      real(c_double) Rsqjk
      real(c_double) Cmax
      real(c_double) Cmin

      real(c_double) C, fac
c
      C = (2.0_cd*(Rsqik+Rsqjk) - (Rsqik-Rsqjk)**2 - 1.0_cd)
     &  / (1.0_cd - (Rsqik-Rsqjk)**2)
      fac = (C - Cmin) / (Cmax - Cmin)
      if(fac.ge.1.0_cd) then
         screen = 1.0_cd
      elseif(fac.le.0.0_cd) then
         screen = 0.0_cd
      else
         screen = ( 1.0_cd - (1.0_cd - fac)**4 )**2
      endif
      return
      end

      real(c_double) function rhoi
     &                 (t,rho0,rho2,rho1s,rho2s,rho3s,rho3a)
c
c     Computation of Rho_Bar
c
      use, intrinsic :: iso_c_binding
      implicit none

      integer(c_int), parameter :: cd = c_double
      real(c_double) t(3)
      real(c_double) rho0
      real(c_double) rho2
      real(c_double) rho1s(3)
      real(c_double) rho2s(3,3)
      real(c_double) rho3s(3,3,3)
      real(c_double) rho3a(3)

      integer(c_int) ix, iy, iz
      real(c_double) Rho1P, Rho2P, Rho3P, Rho33, RhoFactor, GAMMA
c
      Rho1P = 0.0_cd
      Rho3P = 0.0_cd
      Rho33 = 0.0_cd
c
      RhoFactor = Rho0 * Rho0
      Rho2P     = - Rho2 * Rho2 / 3.0_cd
c
      do ix = 1, 3
         Rho1P = Rho1P + Rho1S(ix) * Rho1S(ix)
         Rho33 = Rho33 + Rho3A(ix) * Rho3A(ix)
      do iy = 1, 3
         Rho2P = Rho2P + Rho2S(ix,iy) * Rho2S(ix,iy)
      do iz = 1, 3
         Rho3P = Rho3P + Rho3S(ix,iy,iz) * Rho3S(ix,iy,iz)
      enddo
      enddo
      enddo
      Rho3P = Rho3P - 0.6_cd * Rho33
c
      if(RhoFactor.eq.0.0_cd) then
         RhoFactor = 1.e-20_cd
      endif
      GAMMA = t(1) * Rho1P / RhoFactor
     &      + t(2) * Rho2P / RhoFactor
     &      + t(3) * Rho3P / RhoFactor
c
      rhoi = Rho0 * 2.0_cd / ( 1.0_cd + dexp(-GAMMA) )
c     rhoi = Rho0 * dsqrt ( 1.0_cd + GAMMA )
      return
      end


      subroutine Read_Cutoff(paramfile_name,nmstrlen,
     &           Ncomp,Icomp,Rcutoff,Ierr)
      use, intrinsic :: iso_c_binding
      implicit none

      integer(c_int),                 intent(in)  :: nmstrlen
      character(len=nmstrlen), intent(in)  :: paramfile_name
      integer(c_int),                 intent(in)  :: Ncomp
      character(len=2),        intent(in)  :: Icomp(Ncomp)
      real(c_double),        intent(out) :: Rcutoff
      integer(c_int),                 intent(out) :: Ierr

      character(len=80) line1,line2
c
c reset input file to read from beginning
c
      open(unit=2,file=paramfile_name,status='old',err=4448)

      line1 = " "
      line2 = " "
      if ( Ncomp.EQ.1 ) then

         do while ( Icomp(1).NE.line1(1:2) )
            read(2,'(a)',end=4444) line1
            read(2,'(a)',end=4444) line2
            if ( line1(1:2).EQ.'XX' ) goto 4444
         end do
         read(line2(8:12),'(f5.2)') Rcutoff

      else if ( Ncomp.EQ.2 ) then

         do while ( (Icomp(1)//'-'//Icomp(2).NE.line1(1:5)).AND.
     &              (Icomp(2)//'-'//Icomp(1).NE.line1(1:5)) )
            read(2,'(a)',end=4445) line1
            read(2,'(a)',end=4445) line2
            if ( line2(1:7).EQ.'Ternary' ) goto 4445
         end do
         read(line2(51:55),'(f5.2)') Rcutoff

      else if ( Ncomp.EQ.3 ) then

         do while( line2(1:7).NE.'Ternary' )
            read(2,'(a)') line1
            read(2,'(a)') line2
         end do
         do while( Icomp(1)//'-'//Icomp(2)//'-'//Icomp(3).NE.
     &             line1(1:8) )
            read(2,'(a)',end=4446) line1
         end do
         read(line1(36:40),'(f5.2)') Rcutoff

      else if ( Ncomp.GT.3 ) then

         do while( line1(1:14).NE.'Multicomponent' )
            read(2,'(a)',end=4447) line1
         end do
         read(2,'(10x,f5.2)') Rcutoff

      end if

      close(2)
      Ierr = KIM_STATUS_OK
      return

 4444 continue
      close(2)
      Ierr = KIM_STATUS_FAIL
      write(*,'(/,a,a,/)') ' **** ERROR : Non-Registered Element ',
     &      Icomp(1)
      return
 4445 continue
      close(2)
      Ierr = KIM_STATUS_FAIL
      write(*,'(/,a,a,/)') ' **** ERROR : Non-Registered Alloy ',
     &      Icomp(1)//'-'//Icomp(2)
      return
 4446 continue
      close(2)
      Ierr = KIM_STATUS_FAIL
      write(*,'(/,a,a,/)') ' **** ERROR : Non-Registered Ternary ',
     &      Icomp(1)//'-'//Icomp(2)//'-'//Icomp(3)
      return
 4447 continue
      close(2)
      Ierr = KIM_STATUS_FAIL
      write(*,'(/,a,/)') ' **** ERROR : No Rcutoff for Multicomponent '
      return
 4448 continue
      close(2)
      Ierr = KIM_STATUS_FAIL
      write(*,'(/,a,/)') ' **** ERROR : Unable to open file '
      return


      end subroutine Read_Cutoff

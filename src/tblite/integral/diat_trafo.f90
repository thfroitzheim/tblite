! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/integral/diat_trafo.f90
!> Evaluation of the diatomic scaled overlap
module tblite_integral_diat_trafo
   use iso_fortran_env, only: output_unit


   use mctc_env, only : wp
   use tblite_blas, only: gemm

   implicit none
   private

   public :: relvec, diat_trafo, diat_trafo_grad

   !> Dimension of trafo matrix for highest angular momentum.
   integer, parameter :: ndim(7) = [1, 4, 9, 16, 25, 36, 49]

contains

   !> Transformation to the diatomic frame and back: 
   !pure 
   subroutine diat_trafo(block_overlap, vec_diat_trafo, ksig, kpi, kdel, maxl)
      !> Diatomic block of CGTOs to be transformed (+ scaled)
      real(wp),intent(inout)    :: block_overlap(:,:)
      !> Transformation vector for the diatomic frame
      real(wp),intent(in)       :: vec_diat_trafo(3)
      !> Scaling parameters for different bonding contributions
      real(wp),intent(in)       :: ksig, kpi, kdel
      !> Highest angular momentum between the two shells
      integer,intent(in)        :: maxl

      integer  :: trafo_dim 
      real(wp) :: trafomat(ndim(maxl+1),ndim(maxl+1))
      real(wp), allocatable :: tmp(:,:), transformed_s(:,:)

      trafo_dim = ndim(maxl+1)

      ! Allocate arrays in the needed size
      allocate(transformed_s(trafo_dim,trafo_dim), tmp(trafo_dim,trafo_dim), source=0.0_wp)

      ! 1. Setup the transformation matrix
      call harmtr(maxl, vec_diat_trafo, trafomat)

      ! 2. Transform the overlap submatrix to the diatomic frame
      ! trans_block_s = O^T * S * O
      if (maxl > 0) then
         call gemm(amat=trafomat,bmat=block_overlap,cmat=tmp,transa='T',transb='N')
         call gemm(amat=tmp,bmat=trafomat,cmat=transformed_s,transa='N',transb='N')
      else
         transformed_s(1,1) = block_overlap(1,1)
      endif

      ! 3. Scale elements in the diatomic frame
      call scale_diatomic_frame(transformed_s, ksig, kpi, kdel, maxl) 
      
      block_overlap = 0.0_wp
      ! 4. Transform the overlap submatrix back to original frame
      ! block_overlap = O * S * O^T
      if (maxl > 0) then
         call gemm(amat=trafomat,bmat=transformed_s,cmat=tmp,transa='N',transb='N')
         call gemm(amat=tmp,bmat=trafomat,cmat=block_overlap,transa='N',transb='T')
      else
         block_overlap(1,1) = transformed_s(1,1)
      endif

   end subroutine diat_trafo


   !> Gradient of Transformation to the diatomic frame and back: 
   !pure 
   subroutine diat_trafo_grad(block_overlap, block_doverlap, vec_diat_trafo, ksig, kpi, kdel, maxl)
      !> Diatomic block of CGTO overlap to be transformed (+ scaled)
      real(wp),intent(inout)    :: block_overlap(:,:)
      !> Derivative of diatomic block of CGTO overlap to be transformed (+ scaled)
      real(wp),intent(inout)    :: block_doverlap(:,:,:)
      !> Transformation vector for the diatomic frame
      real(wp),intent(in)       :: vec_diat_trafo(3)
      !> Scaling parameters for different bonding contributions
      real(wp),intent(in)       :: ksig, kpi, kdel
      !> Highest angular momentum between the two shells
      integer,intent(in)        :: maxl
      
      integer :: ic, trafo_dim
      real(wp) :: trafomat(ndim(maxl+1),ndim(maxl+1)), dtrafomat(3,ndim(maxl+1),ndim(maxl+1))
      real(wp), allocatable :: eff_tra_mat(:,:), eff_dtra_mat(:,:,:)
      real(wp), allocatable :: eff_block_overlap(:,:), eff_block_doverlap(:,:,:)
      real(wp), allocatable :: tmp(:,:), tmp2(:,:), interm_oso(:,:), &
      & interm_doso(:,:,:), interm_odso(:,:,:), interm_osdo(:,:,:)

      trafo_dim = ndim(maxl+1)

      trafomat = 0.0_wp
      dtrafomat = 0.0_wp

      ! Allocate arrays in the needed size
      allocate(interm_oso(trafo_dim,trafo_dim), interm_doso(3,trafo_dim,trafo_dim), &
      & interm_osdo(3,trafo_dim,trafo_dim), interm_odso(3,trafo_dim,trafo_dim), &
      & tmp(trafo_dim,trafo_dim), tmp2(trafo_dim,trafo_dim), source=0.0_wp)

      ! 1. Setup the transformation matrix and its derivative
      call d_harmtr(maxl, vec_diat_trafo, trafomat, dtrafomat)

      ! 2. Transform the overlap submatrix to the diatomic frame
      ! interm_s = matmul(matmul(transpose(trafomat), block_overlap),trafomat)
      if (maxl > 0) then
         ! interm_oso = O^T * S * O
         call gemm(amat=trafomat,bmat=block_overlap,cmat=tmp,transa='T',transb='N')
         call gemm(amat=tmp,bmat=trafomat,cmat=interm_oso,transa='N',transb='N')
         do ic = 1, 3
            ! interm_doso = dO^T * S * O
            call gemm(amat=dtrafomat(ic,:,:),bmat=block_overlap,cmat=tmp,transa='T',transb='N')
            call gemm(amat=tmp,bmat=trafomat,cmat=interm_doso(ic,:,:),transa='N',transb='N')
            ! interm_osdo = O^T * S * dO
            call gemm(amat=trafomat,bmat=block_overlap,cmat=tmp,transa='T',transb='N')
            call gemm(amat=tmp,bmat=dtrafomat(ic,:,:),cmat=interm_osdo(ic,:,:),transa='N',transb='N')
            ! interm_odso = O^T * dS * O
            call gemm(amat=trafomat,bmat=block_doverlap(ic,:,:),cmat=tmp,transa='T',transb='N')
            call gemm(amat=tmp,bmat=trafomat,cmat=interm_odso(ic,:,:),transa='N',transb='N')
         end do
      else
         interm_oso(1,1) = block_overlap(1,1)
         interm_doso(:,1,1) = 0.0_wp
         interm_odso(:,1,1) = block_doverlap(:,1,1)
         interm_osdo(:,1,1) = 0.0_wp
      endif

      ! 3. Scale overlap and each dimension of the derivative in the diatomic frame
      call scale_diatomic_frame(interm_oso, ksig, kpi, kdel, maxl) 
      do ic = 1, 3
         call scale_diatomic_frame(interm_doso(ic,:,:), ksig, kpi, kdel, maxl)
         call scale_diatomic_frame(interm_osdo(ic,:,:), ksig, kpi, kdel, maxl) 
         call scale_diatomic_frame(interm_odso(ic,:,:), ksig, kpi, kdel, maxl) 
      end do

      ! 4. Transform diatomic frame quantities (S', (dOSO)', and (OdSO)') back to original frame
      if (maxl > 0) then
         ! block_overlap = O * S' * O^T
         call gemm(amat=trafomat,bmat=interm_oso,cmat=tmp,transa='N',transb='N')
         call gemm(amat=tmp,bmat=trafomat,cmat=block_overlap,transa='N',transb='T')

         do ic = 1, 3
            ! block_doverlap = dO * S' * O^T + O * S' * dO^T
            call gemm(amat=dtrafomat(ic,:,:),bmat=interm_oso,cmat=tmp,transa='N',transb='N')
            call gemm(alpha=1.0_wp,amat=tmp,bmat=trafomat,cmat=block_doverlap(ic,:,:),transa='N',transb='T')
            ! block_doverlap = dO * S' * O^T + O * S' * dO^T
            call gemm(amat=trafomat,bmat=interm_oso,cmat=tmp,transa='N',transb='N')
            call gemm(alpha=1.0_wp,amat=tmp,bmat=dtrafomat(ic,:,:),beta=1.0_wp,&
            &cmat=block_doverlap(ic,:,:),transa='N',transb='T')
            
            !eff_block_doverlap(ic,:,:) = tmp2 !+ transpose(tmp2)

            ! block_doverlap += O * ((dOSO)' + (OdSO)' + (OSdO)') * O^T 

            !tmp2 = interm_doso(ic,:,:) + interm_odso(ic,:,:) + interm_osdo(ic,:,:)
            call gemm(amat=trafomat,bmat=interm_doso(ic,:,:),cmat=tmp,transa='N',transb='N')
            call gemm(alpha=1.0_wp,amat=tmp,bmat=trafomat,beta=1.0_wp,cmat=block_doverlap(ic,:,:),transa='N',transb='T')

            call gemm(amat=trafomat,bmat=interm_odso(ic,:,:),cmat=tmp,transa='N',transb='N')
            call gemm(amat=tmp,bmat=trafomat,beta=1.0_wp,cmat=block_doverlap(ic,:,:),transa='N',transb='T')

            call gemm(amat=trafomat,bmat=interm_osdo(ic,:,:),cmat=tmp,transa='N',transb='N')
            call gemm(alpha=1.0_wp,amat=tmp,bmat=trafomat,beta=1.0_wp,cmat=block_doverlap(ic,:,:),transa='N',transb='T')
         end do
      else
         block_overlap(1,1) = interm_oso(1,1)
         block_doverlap(:,1,1) = interm_odso(:,1,1)
      endif

   end subroutine diat_trafo_grad


   subroutine relvec(vec, rkl, veckl)
      !> Original vector between atoms A and B
      real(wp), intent(in)             :: vec(3)
      !> Distance between the two atoms
      real(wp), intent(in)             :: rkl
      !> Normalized vector from atom k to atom l
      real(wp), intent(out)            :: veckl(3)

      real(wp), parameter              :: eps = 4.0e-08_wp
      real(wp)                         :: sq

      veckl(1:3) = vec(1:3) !/ rkl
      !if ( abs(1.0_wp-abs(veckl(1))) .lt. eps ) then
      !   veckl(1) = sign(1.0_wp,veckl(1))
      !   veckl(2) = 0.0_wp
      !   veckl(3) = 0.0_wp
      !else if ( abs(1.0_wp-abs(veckl(2))) .lt. eps ) then
      !   veckl(1) = 0.0_wp
      !   veckl(2) = sign(1.0_wp,veckl(2))
      !   veckl(3) = 0.0_wp
      !else if ( abs(1.0_wp-abs(veckl(3))) .lt. eps ) then
      !   veckl(1) = 0.0_wp
      !   veckl(2) = 0.0_wp
      !   veckl(3) = sign(1.0_wp,veckl(3))
      !else if ( (abs(veckl(1)) .lt. eps) .and. .not. eff_equality(veckl(1),0.0_wp) ) then
      !   veckl(1) = 0.0_wp
      !   sq = sqrt( veckl(2)**2 + veckl(3)**2 )
      !   veckl(2) = veckl(2)/sq
      !   veckl(3) = veckl(3)/sq
      !else if ( (abs(veckl(2)) .lt. eps) .and. .not. eff_equality(veckl(2),0.0_wp) ) then
      !   veckl(2) = 0.0_wp
      !   sq = sqrt( veckl(1)**2 + veckl(3)**2 )
      !   veckl(1) = veckl(1)/sq
      !   veckl(3) = veckl(3)/sq
      !else if ( (abs(veckl(3)) .lt. eps) .and. .not. eff_equality(veckl(3),0.0_wp) ) then
      !   veckl(3) = 0.0_wp
      !   sq = sqrt(veckl(1)**2 + veckl(2)**2)
      !   veckl(1) = veckl(1)/sq
      !   veckl(2) = veckl(2)/sq
      !endif

   end subroutine relvec

   logical pure function eff_equality(num1, num2)
      !> Numbers to compare
      real(wp), intent(in) :: num1, num2
      
      ! Logical deciding if numbers are (almost) equal or not
      eff_equality = (abs( num1 - num2 ) .le. 1.0e-12_wp)
   end function eff_equality

   pure subroutine harmtr(maxl,veckl,trafomat)
      !> Maximum angular momentum
      integer, intent(in)  :: maxl
      !> Normalized vector from atom k to atom l
      real(wp), intent(in) :: veckl(3)
      !> Transformation matrix
      real(wp), intent(out) :: trafomat(ndim(maxl+1),ndim(maxl+1))

      real(wp) :: cos2p, cos2t, cosp, cost, sin2p, sin2t, sinp, sint, sqrt3, len

      trafomat = 0.0_wp
      ! -----------------------------
      ! *** s functions (trafomat(1x1)) ***
      ! -----------------------------
      trafomat(1,1) = 1.0

      if ( maxl == 0 ) return
      
      len = sqrt(veckl(1)**2 + veckl(2)**2 + veckl(3)**2)
      ! Prepare spherical coordinats
      cost = veckl(3) / len
      if ( abs(cost) .eq. 1.0_wp ) then
         sint = 0.0_wp
         cosp = 1.0_wp
         sinp = 0.0_wp
      else if ( abs(cost) .eq. 0.0_wp ) then
         sint = 1.0_wp
         cosp = veckl(1)/len
         sinp = veckl(2)/len
      else
         !sint = SQRT(1.0_wp-COST**2)
         sint = SQRT(veckl(1)**2 + veckl(2)**2)
         cosp = veckl(1)/SINT
         sinp = veckl(2)/SINT
      endif

      ! -----------------------------
      ! *** p functions (trafomat(4x4)) ***
      ! -----------------------------

      ! tblite ordering with adapted column ordering
      ! 1st index:
      ! MSINDO defintion of p function ordering is converted to
      ! tblite definition of p function ordering. E.g. for first entry:
      ! trafomat(2,:)_MSINDO -> trafomat(px,:) -> trafomat(4:)_tblite
      ! 2nd index:
      ! Final ordering of p functions (see below) corresponds to the
      ! tblite ordering of p functions. For the second index, the ordering
      ! 3, 4, 2 holds, corresponding (in MSINDO convention)
      ! to the tblite ordering of p functions, i.e.
      ! y, z, x
      trafomat(4,3) = SINT*COSP
      trafomat(2,3) = SINT*SINP
      trafomat(3,3) = COST
      trafomat(4,4) = COST*COSP
      trafomat(2,4) = COST*SINP
      trafomat(3,4) = -SINT
      trafomat(4,2) = -SINP
      trafomat(2,2) = COSP
      trafomat(3,2) = 0.0_wp

      if ( maxl <= 1 ) return

      ! -----------------------------
      ! *** d functions (trafomat(9x9)) ***
      ! -----------------------------

      COS2T = COST**2 - SINT**2
      SIN2T = 2.0_wp * SINT*COST
      COS2P = COSP**2 - SINP**2
      SIN2P = 2.0_wp * SINP*COSP
      SQRT3 = SQRT(3.0_wp)

      ! Original MSINDO ordering
      ! The MSINDO d SAO ordering corresponds to the
      ! tblite ordering of d SAOs
      trafomat(5,5) = (3.0_wp * COST**2 - 1.0_wp) * 0.5_wp
      trafomat(6,5) = SQRT3*SIN2T*COSP*0.5_wp
      trafomat(7,5) = SQRT3*SIN2T*SINP*0.5_wp
      trafomat(8,5) = SQRT3*SINT**2*COS2P*0.5_wp
      trafomat(9,5) = SQRT3*SINT**2*SIN2P*0.5_wp
      trafomat(5,6) = -SQRT3*SIN2T*0.5_wp
      trafomat(6,6) = COS2T*COSP
      trafomat(7,6) = COS2T*SINP
      trafomat(8,6) = SIN2T*COS2P*0.5_wp
      trafomat(9,6) = SIN2T*SIN2P*0.5_wp
      trafomat(5,7) = 0.0_wp
      trafomat(6,7) = -COST*SINP
      trafomat(7,7) = COST*COSP
      trafomat(8,7) = -SINT*SIN2P
      trafomat(9,7) = SINT*COS2P
      trafomat(5,8) = SQRT3*SINT**2 * 0.5_wp
      trafomat(6,8) = -SIN2T*COSP*0.5_wp
      trafomat(7,8) = -SIN2T*SINP*0.5_wp
      trafomat(8,8) = (1.0_wp + COST**2) * COS2P * 0.5_wp
      trafomat(9,8) = (1.0_wp + COST**2) * SIN2P * 0.5_wp
      trafomat(5,9) = 0.0_wp
      trafomat(6,9) = SINT*SINP
      trafomat(7,9) = -SINT*COSP
      trafomat(8,9) = -COST*SIN2P
      trafomat(9,9) = COST*COS2P

      if ( maxl <= 2 ) return

      ! -----------------------------
      ! *** f functions (trafomat(16x16)) ***
      ! -----------------------------

      ! f-functions are not transformed to the diatomic frame!
      trafomat(10,10) = 1.0_wp
      trafomat(11,11) = 1.0_wp
      trafomat(12,12) = 1.0_wp
      trafomat(13,13) = 1.0_wp
      trafomat(14,14) = 1.0_wp
      trafomat(15,15) = 1.0_wp
      trafomat(16,16) = 1.0_wp

      if ( maxl <= 3 ) return

      ! -----------------------------
      ! *** g functions (trafomat(25x25)) ***
      ! -----------------------------

      ! g-functions are not transformed to the diatomic frame!
      trafomat(17,17) = 1.0_wp
      trafomat(18,18) = 1.0_wp
      trafomat(19,19) = 1.0_wp
      trafomat(20,20) = 1.0_wp
      trafomat(21,21) = 1.0_wp
      trafomat(22,22) = 1.0_wp
      trafomat(23,23) = 1.0_wp
      trafomat(24,24) = 1.0_wp
      trafomat(25,25) = 1.0_wp

      if ( maxl <= 4 ) return

   end subroutine harmtr

   !pure 
   subroutine d_harmtr(maxl,veckl,trafomat, dtrafomat)
      !> Maximum angular momentum
      integer, intent(in)  :: maxl
      !> Normalized vector from atom k to atom l
      real(wp), intent(in) :: veckl(3)
      !> Transformation matrix
      real(wp), intent(out) :: trafomat(ndim(maxl+1),ndim(maxl+1))
      !> Transformation matrix
      real(wp), intent(out) :: dtrafomat(3,ndim(maxl+1),ndim(maxl+1))
      
      !> Derivative of transformation matrix w.r.t. theta
      real(wp) :: trafomat_dt(ndim(maxl+1),ndim(maxl+1))
      !> Derivative of transformation matrix w.r.t. phi
      real(wp) :: trafomat_dp(ndim(maxl+1),ndim(maxl+1))
      !!> Derivative of transformation matrix w.r.t. phi
      !real(wp) :: trafomat_dpy(16,16)

      integer :: i,j

      real(wp) :: cos2p, cos2t, cosp, cost, sin2p, sin2t, sinp, sint, sqrt3, len
      real(wp) :: dcos2t, dsin2t, dcos2p, dsin2p, dpdx, dpdy, dpdz, dtdx, dtdy, dtdz

      trafomat = 0.0_wp
      trafomat_dt = 0.0_wp
      trafomat_dp = 0.0_wp
      !trafomat_dpy = 0.0_wp

      ! -----------------------------
      ! *** s functions (trafomat(1x1)) ***
      ! -----------------------------
      trafomat(1,1) = 1.0
      dtrafomat(1, 1, 1) = 0.0_wp

      if ( maxl == 0 ) return

      len = sqrt(veckl(1)**2 + veckl(2)**2 + veckl(3)**2)

      ! Prepare spherical coordinats
      cost = veckl(3)/len
      if ( abs(cost) .eq. 1.0_wp ) then
         sint = 0.0_wp
         cosp = 1.0_wp
         sinp = 0.0_wp
      else if ( abs(cost) .eq. 0.0_wp ) then
         sint = 1.0_wp
         cosp = veckl(1)
         sinp = veckl(2)
      else
         !sint = SQRT(1.0_wp-COST**2)
         sint = SQRT(veckl(1)**2 + veckl(2)**2)
         cosp = veckl(1)/SINT
         sinp = veckl(2)/SINT
      endif

      ! Prepare sperical coordinate derivative
      dpdx = -veckl(2) / (veckl(1)**2 + veckl(2)**2) ! SIGN????
      dpdy = veckl(1) / (veckl(1)**2 + veckl(2)**2)
      dpdz = 0.0_wp
      dtdx = veckl(1)*veckl(3) / (SQRT(veckl(1)**2 + veckl(2)**2)*len**2)
      dtdy = veckl(2)*veckl(3) / (SQRT(veckl(1)**2 + veckl(2)**2)*len**2)
      dtdz = SQRT(veckl(1)**2+veckl(2)**2) / len**2
      !dtdz = -1.0_wp / SQRT(1.0_wp - veckl(3)**2)
      
      ! -----------------------------
      ! *** p functions (trafomat(4x4)) ***
      ! -----------------------------

      ! tblite ordering with adapted column ordering
      ! 1st index:
      ! MSINDO defintion of p function ordering is converted to
      ! tblite definition of p function ordering. E.g. for first entry:
      ! trafomat(2,:)_MSINDO -> trafomat(px,:) -> trafomat(4:)_tblite
      ! 2nd index:
      ! Final ordering of p functions (see below) corresponds to the
      ! tblite ordering of p functions. For the second index, the ordering
      ! 3, 4, 2 holds, corresponding (in MSINDO convention)
      ! to the tblite ordering of p functions, i.e.
      ! y, z, x
      trafomat(4,3) = SINT*COSP
      trafomat(2,3) = SINT*SINP
      trafomat(3,3) = COST
      trafomat(4,4) = COST*COSP
      trafomat(2,4) = COST*SINP
      trafomat(3,4) = -SINT
      trafomat(4,2) = -SINP
      trafomat(2,2) = COSP
      trafomat(3,2) = 0.0_wp

      ! derivative w.r.t. theta
      trafomat_dt(4,3) = COST*COSP
      trafomat_dt(2,3) = COST*SINP
      trafomat_dt(3,3) = -SINT
      trafomat_dt(4,4) = -SINT*COSP
      trafomat_dt(2,4) = -SINT*SINP
      trafomat_dt(3,4) = -COST
      trafomat_dt(4,2) = 0.0_wp
      trafomat_dt(2,2) = 0.0_wp
      trafomat_dt(3,2) = 0.0_wp

      ! derivative w.r.t. phi
      trafomat_dp(4,3) = -SINT*SINP
      trafomat_dp(2,3) = SINT*COSP
      trafomat_dp(3,3) = 0.0_wp
      trafomat_dp(4,4) = -COST*SINP
      trafomat_dp(2,4) = COST*COSP
      trafomat_dp(3,4) = 0.0_wp
      trafomat_dp(4,2) = -COSP
      trafomat_dp(2,2) = -SINP
      trafomat_dp(3,2) = 0.0_wp
 
      if ( maxl <= 1 ) then 
         dtrafomat(1, 1:4, 1:4) = dpdx * trafomat_dp(1:4, 1:4) + dtdx * trafomat_dt(1:4, 1:4) 
         dtrafomat(2, 1:4, 1:4) = dpdy * trafomat_dp(1:4, 1:4) + dtdy * trafomat_dt(1:4, 1:4) 
         dtrafomat(3, 1:4, 1:4) = dpdz * trafomat_dp(1:4, 1:4) + dtdz * trafomat_dt(1:4, 1:4) 
         return
      end if

      ! -----------------------------
      ! *** d functions (trafomat(9x9)) ***
      ! -----------------------------

      COS2T = COST**2 - SINT**2
      SIN2T = 2.0_wp * SINT*COST
      COS2P = COSP**2 - SINP**2
      SIN2P = 2.0_wp * SINP*COSP
      SQRT3 = SQRT(3.0_wp)

      DCOS2T = -2.0_wp * SIN2T
      DSIN2T =  2.0_wp * COS2T
      DCOS2P = -2.0_wp * SIN2P
      DSIN2P =  2.0_wp * COS2P

      ! Original MSINDO ordering
      ! The MSINDO d SAO ordering corresponds to the
      ! tblite ordering of d SAOs
      trafomat(5,5) = (3.0_wp * COST**2 - 1.0_wp) * 0.5_wp
      trafomat(6,5) = SQRT3*SIN2T*COSP*0.5_wp
      trafomat(7,5) = SQRT3*SIN2T*SINP*0.5_wp
      trafomat(8,5) = SQRT3*SINT**2*COS2P*0.5_wp
      trafomat(9,5) = SQRT3*SINT**2*SIN2P*0.5_wp
      trafomat(5,6) = -SQRT3*SIN2T*0.5_wp
      trafomat(6,6) = COS2T*COSP
      trafomat(7,6) = COS2T*SINP
      trafomat(8,6) = SIN2T*COS2P*0.5_wp
      trafomat(9,6) = SIN2T*SIN2P*0.5_wp
      trafomat(5,7) = 0.0_wp
      trafomat(6,7) = -COST*SINP
      trafomat(7,7) = COST*COSP
      trafomat(8,7) = -SINT*SIN2P
      trafomat(9,7) = SINT*COS2P
      trafomat(5,8) = SQRT3*SINT**2 * 0.5_wp
      trafomat(6,8) = -SIN2T*COSP*0.5_wp
      trafomat(7,8) = -SIN2T*SINP*0.5_wp
      trafomat(8,8) = (1.0_wp + COST**2) * COS2P * 0.5_wp
      trafomat(9,8) = (1.0_wp + COST**2) * SIN2P * 0.5_wp
      trafomat(5,9) = 0.0_wp
      trafomat(6,9) = SINT*SINP
      trafomat(7,9) = -SINT*COSP
      trafomat(8,9) = -COST*SIN2P
      trafomat(9,9) = COST*COS2P
      
      ! derivative w.r.t. theta
      trafomat_dt(5,5) = -3.0_wp*SIN2T*0.5_wp
      trafomat_dt(6,5) = SQRT3*DSIN2T*COSP*0.5_wp
      trafomat_dt(7,5) = SQRT3*DSIN2T*SINP*0.5_wp
      trafomat_dt(8,5) = SQRT3*SIN2T*COS2P*0.5_wp
      trafomat_dt(9,5) = SQRT3*SIN2T*SIN2P*0.5_wp
      trafomat_dt(5,6) = -SQRT3*DSIN2T*0.5_wp
      trafomat_dt(6,6) = DCOS2T*COSP
      trafomat_dt(7,6) = DCOS2T*SINP
      trafomat_dt(8,6) = DSIN2T*COS2P*0.5_wp
      trafomat_dt(9,6) = DSIN2T*SIN2P*0.5_wp
      trafomat_dt(5,7) = 0.0_wp
      trafomat_dt(6,7) = SINT*SINP
      trafomat_dt(7,7) = -SINT*COSP
      trafomat_dt(8,7) = -COST*SIN2P
      trafomat_dt(9,7) = COST*COS2P
      trafomat_dt(5,8) = SQRT3*SIN2T*0.5_wp
      trafomat_dt(6,8) = -DSIN2T*COSP*0.5_wp
      trafomat_dt(7,8) = -DSIN2T*SINP*0.5_wp
      trafomat_dt(8,8) = -SIN2T*COS2P*0.5_wp
      trafomat_dt(9,8) = -SIN2T*SIN2P*0.5_wp
      trafomat_dt(5,9) = 0.0_wp
      trafomat_dt(6,9) = COST*SINP
      trafomat_dt(7,9) = -COST*COSP
      trafomat_dt(8,9) = SINT*SIN2P
      trafomat_dt(9,9) = -SINT*COS2P

      ! derivative w.r.t. phi
      trafomat_dp(5,5) = 0.0_wp
      trafomat_dp(6,5) = -SQRT3*SIN2T*SINP*0.5_wp
      trafomat_dp(7,5) = SQRT3*SIN2T*COSP*0.5_wp
      trafomat_dp(8,5) = SQRT3*SINT**2*DCOS2P*0.5_wp
      trafomat_dp(9,5) = SQRT3*SINT**2*DSIN2P*0.5_wp
      trafomat_dp(5,6) = 0.0_wp
      trafomat_dp(6,6) = -COS2T*SINP
      trafomat_dp(7,6) = COS2T*COSP
      trafomat_dp(8,6) = SIN2T*DCOS2P*0.5_wp
      trafomat_dp(9,6) = SIN2T*DSIN2P*0.5_wp
      trafomat_dp(5,7) = 0.0_wp
      trafomat_dp(6,7) = -COST*COSP
      trafomat_dp(7,7) = -COST*SINP
      trafomat_dp(8,7) = -SINT*DSIN2P
      trafomat_dp(9,7) = SINT*DCOS2P
      trafomat_dp(5,8) = 0.0_wp
      trafomat_dp(6,8) = SIN2T*SINP*0.5_wp
      trafomat_dp(7,8) = -SIN2T*COSP*0.5_wp
      trafomat_dp(8,8) = (1.0+COST**2)*DCOS2P*0.5_wp
      trafomat_dp(9,8) = (1.0+COST**2)*DSIN2P*0.5_wp
      trafomat_dp(5,9) = 0.0_wp
      trafomat_dp(6,9) = SINT*COSP
      trafomat_dp(7,9) = SINT*SINP
      trafomat_dp(8,9) = -COST*DSIN2P
      trafomat_dp(9,9) = COST*DCOS2P

      if ( maxl <= 2 ) then 
         ! Transform to cartesian coordinates
         dtrafomat(1, 1:9, 1:9) = dpdx * trafomat_dp(1:9, 1:9) + dtdx * trafomat_dt(1:9, 1:9)
         dtrafomat(2, 1:9, 1:9) = dpdy * trafomat_dp(1:9, 1:9) + dtdy * trafomat_dt(1:9, 1:9)
         dtrafomat(3, 1:9, 1:9) = dpdz * trafomat_dp(1:9, 1:9) + dtdz * trafomat_dt(1:9, 1:9)
         return
      end if


      ! -----------------------------
      ! *** f functions (trafomat(16x16)) ***
      ! -----------------------------

      ! f-functions are not transformed to the diatomic frame!
      ! Hence, the derivative remaints 0!
      trafomat(10,10) = 1.0_wp
      trafomat(11,11) = 1.0_wp
      trafomat(12,12) = 1.0_wp
      trafomat(13,13) = 1.0_wp
      trafomat(14,14) = 1.0_wp
      trafomat(15,15) = 1.0_wp
      trafomat(16,16) = 1.0_wp

      if ( maxl <= 3 ) then 
         ! Transform to cartesian coordinates
         dtrafomat(1, 1:16, 1:16) = dpdx * trafomat_dp(1:16, 1:16) + dtdx * trafomat_dt(1:16, 1:16)
         dtrafomat(2, 1:16, 1:16) = dpdy * trafomat_dp(1:16, 1:16) + dtdy * trafomat_dt(1:16, 1:16)
         dtrafomat(3, 1:16, 1:16) = dpdz * trafomat_dp(1:16, 1:16) + dtdz * trafomat_dt(1:16, 1:16)
         return
      end if


      ! -----------------------------
      ! *** g functions (trafomat(25x25)) ***
      ! -----------------------------

      ! g-functions are not transformed to the diatomic frame!
      trafomat(17,17) = 1.0_wp
      trafomat(18,18) = 1.0_wp
      trafomat(19,19) = 1.0_wp
      trafomat(20,20) = 1.0_wp
      trafomat(21,21) = 1.0_wp
      trafomat(22,22) = 1.0_wp
      trafomat(23,23) = 1.0_wp
      trafomat(24,24) = 1.0_wp
      trafomat(25,25) = 1.0_wp

      if ( maxl <= 4 ) then 
         ! Transform to cartesian coordinates
         dtrafomat(1, 1:25, 1:25) = dpdx * trafomat_dp(1:25, 1:25) + dtdx * trafomat_dt(1:25, 1:25)
         dtrafomat(2, 1:25, 1:25) = dpdy * trafomat_dp(1:25, 1:25) + dtdy * trafomat_dt(1:25, 1:25)
         dtrafomat(3, 1:25, 1:25) = dpdz * trafomat_dp(1:25, 1:25) + dtdz * trafomat_dt(1:25, 1:25)
         return
      end if

   end subroutine d_harmtr

   pure subroutine scale_diatomic_frame(diat_mat, ksig, kpi, kdel, maxl)
      !> Block matrix in the diatomic frame to be scaled
      real(wp),intent(inout)    :: diat_mat(:,:)
      !> Scaling parameters for different bonding contributions
      real(wp),intent(in)       :: ksig, kpi, kdel
      !> Highest angular momentum between the two shells
      integer,intent(in)        :: maxl

      diat_mat(1,1) = diat_mat(1,1)*ksig ! Sigma bond s   <-> s
      if (maxl > 0) then
         diat_mat(1,3) = diat_mat(1,3)*ksig ! Sigma bond s   <-> pz
         diat_mat(3,1) = diat_mat(3,1)*ksig ! Sigma bond pz  <-> s
         diat_mat(3,3) = diat_mat(3,3)*ksig ! Sigma bond pz  <-> pz
         diat_mat(4,4) = diat_mat(4,4)*kpi  ! Pi    bond px  <-> px
         diat_mat(2,2) = diat_mat(2,2)*kpi  ! Pi    bond py  <-> py
         if (maxl > 1) then
            diat_mat(5,1) = diat_mat(5,1)*ksig ! Sigma bond dz2 <-> s
            diat_mat(1,5) = diat_mat(1,5)*ksig ! Sigma bond s   <-> dz2
            diat_mat(3,5) = diat_mat(3,5)*ksig ! Sigma bond pz  <-> dz2
            diat_mat(5,3) = diat_mat(5,3)*ksig ! Sigma bond dz2 <-> pz
            diat_mat(5,5) = diat_mat(5,5)*ksig ! Sigma bond dz2 <-> dz2
            diat_mat(4,6) = diat_mat(4,6)*kpi  ! Pi    bond px  <-> dxz
            diat_mat(6,4) = diat_mat(6,4)*kpi  ! Pi    bond dxz <-> px
            diat_mat(6,6) = diat_mat(6,6)*kpi  ! Pi    bond dxz <-> dxz
            diat_mat(2,7) = diat_mat(2,7)*kpi  ! Pi    bond py  <-> dyz
            diat_mat(7,2) = diat_mat(7,2)*kpi  ! Pi    bond dyz <-> py
            diat_mat(7,7) = diat_mat(7,7)*kpi  ! Pi    bond dyz <-> dyz
            diat_mat(8,8) = diat_mat(8,8)*kdel ! Delta bond dx2-y2 <-> dx2-y2
            diat_mat(9,9) = diat_mat(9,9)*kdel ! Delta bond dxy <-> dxy
         endif
         ! f-functions remain unscaled
      endif

   end subroutine scale_diatomic_frame


















subroutine write_2d_matrix(matrix, name, unit, step)
   implicit none
   real(wp), intent(in) :: matrix(:, :)
   character(len=*), intent(in), optional :: name
   integer, intent(in), optional :: unit
   integer, intent(in), optional :: step
   integer :: d1, d2
   integer :: i, j, k, l, istep, iunit

   d1 = size(matrix, dim=1)
   d2 = size(matrix, dim=2)

   if (present(unit)) then
     iunit = unit
   else
     iunit = output_unit
   end if

   if (present(step)) then
     istep = step
   else
     istep = 6
   end if

   if (present(name)) write (iunit, '(/,"matrix printed:",1x,a)') name

   do i = 1, d2, istep
     l = min(i + istep - 1, d2)
     write (iunit, '(/,6x)', advance='no')
     do k = i, l
       write (iunit, '(6x,i7,1x)', advance='no') k
     end do
     write (iunit, '(a)')
     do j = 1, d1
       write (iunit, '(i6)', advance='no') j
       do k = i, l
         write (iunit, '(1x,f13.8)', advance='no') matrix(j, k)
       end do
       write (iunit, '(a)')
     end do
   end do

 end subroutine write_2d_matrix






end module tblite_integral_diat_trafo

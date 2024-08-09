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

!> @dir tblite/exchange/mulliken_kfock
!> Contains the implementation of the Mulliken approximated Fock exchange

!> @file tblite/exchange/matrix.f90
!> Provides an Mulliken approximated Fock exchange interaction

!> Mulliken approximated Fock exchange
module tblite_mulliken_kfock_matrix
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_blas, only: gemm
   use tblite_mulliken_kfock_type, only : mulliken_kfock_type
   use tblite_scf_potential, only : potential_type
   use tblite_container_cache, only : container_cache
   use tblite_exchange_cache, only : exchange_cache
   use tblite_basis_type, only : basis_type
   use tblite_coulomb_charge, only : average_interface, &
      & arithmetic_average, geometric_average, harmonic_average, general_average
   use tblite_wignerseitz, only : wignerseitz_cell
   implicit none
   private

   public :: new_mulliken_exchange_matrix

   type, public, extends(mulliken_kfock_type) :: mulliken_kfock_matrix
      !> Averaged Hubbard parameter for each shell and species
      real(wp), allocatable :: hubbard(:, :, :, :)
      !> Smooth interpolation of arithmetic (0), geometric (1), and harmonic (2) average
      real(wp) :: average
      !> Smoothening exponent (1 = Mataga-Nishimoto, 2 = Klopman-Ohno)
      real(wp) :: gexp
      !> Long-range cutoff
      real(wp) :: rcut
   contains
      !> Evaluate exchange gamma matrix
      procedure :: get_exchange_matrix
      !> Evaluate uncontracted derivatives of exchange gamma matrix
      procedure :: get_exchange_derivs
      !> Calculate exchange contribution to the Fock matrix
      procedure :: get_KFock
      !> Calculate exchange contribution to the gradient
      procedure :: get_KGrad
   end type mulliken_kfock_matrix

   character(len=*), parameter :: label = "Mulliken approximated exchange"

contains


!> Create a new Mulliken approximated exchange container
subroutine new_mulliken_exchange_matrix(self, mol, bas, hubbard, &
   & shell_resolved_hubbard, average, gexp, frscale, omega, lrscale)
   !> Instance of the Mulliken exchange container
   type(mulliken_kfock_matrix), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Description of the basis set
   type(basis_type), intent(in) :: bas
   !> Hubbard parameter for all shells and species
   real(wp), intent(in) :: hubbard(:, :)
   !> Whether the Hubbard parameters are shell-dependent
   logical :: shell_resolved_hubbard
   !> Averaging function for Hubbard parameter of a shell-pair
   real(wp), intent(in) :: average
   !> Smoothening exponent; 1 = Mataga-Nishimoto, 2 = Klopman-Ohno
   real(wp), intent(in) :: gexp
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Optional range separation parameter
   real(wp), intent(in), optional :: omega
   !> Optional long range scaling for range seperated exchange
   real(wp), intent(in), optional :: lrscale

   integer :: isp, jsp, ish, jsh, mshell
   procedure(average_interface), pointer :: averager

   self%label = label

   allocate(self%nsh_id(mol%nid), self%nao_sh(bas%nsh), self%ish_at(mol%nat), &
      & self%iao_sh(bas%nsh))
   self%nao = bas%nao
   self%nsh = bas%nsh
   self%nsh_id = bas%nsh_id
   self%nao_sh = bas%nao_sh
   self%ish_at = bas%ish_at
   self%iao_sh = bas%iao_sh

   self%frscale = frscale
   if (present(omega).and.present(lrscale)) then
      self%omega = omega
      self%lrscale = lrscale
   end if

   self%gexp = gexp
   self%rcut = 10.0_wp

   !> Select averaging scheme
   if(average == 0.0_wp) then
      averager => arithmetic_average
   else if(average == 1.0_wp) then
      averager => geometric_average
   else if(average == 2.0_wp) then
      averager => harmonic_average
   else
      averager => general_average
   end if

   self%shell_resolved_hubbard = shell_resolved_hubbard
   if (shell_resolved_hubbard) then
      mshell = maxval(bas%nsh_id)
      allocate(self%hubbard(mshell, mshell, mol%nid, mol%nid))
      do isp = 1, mol%nid
         do jsp = 1, mol%nid
            self%hubbard(:, :, jsp, isp) = 0.0_wp
            do ish = 1, bas%nsh_id(isp)
               do jsh = 1, bas%nsh_id(jsp)
                  self%hubbard(jsh, ish, jsp, isp) = &
                     & averager(hubbard(ish, isp), hubbard(jsh, jsp), average)
               end do
            end do
         end do
      end do
   else
      allocate(self%hubbard(1, 1, mol%nid, mol%nid))
      do isp = 1, mol%nid
         do jsp = 1, mol%nid
            self%hubbard(1, 1, jsp, isp) = &
               & averager(hubbard(1, isp), hubbard(1, jsp), average)
         end do
      end do
   end if

end subroutine new_mulliken_exchange_matrix


!> Calculate exchange contribution to the Fock matrix
subroutine get_KFock(self, mol, cache, density, overlap, fock)
   !> Instance of the exchange container
   class(mulliken_kfock_matrix), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container with exchange matrix and intermediates
   type(exchange_cache), intent(inout) :: cache
   !> Density matrix
   real(wp),intent(in) :: density(:, :, :)
   !> Overlap matrix
   real(wp),intent(in) :: overlap(:, :)
   !> Fock Matrix contribution
   real(wp),intent(inout) :: fock(:, :, :) 
   
   integer :: spin
   real(wp), allocatable :: tmpA(:,:), tmpB(:,:)

   allocate(tmpA(self%nao, self%nao), tmpB(self%nao, self%nao))

   ! Evaluate the Mulliken approximated exchange Fock matrix contribution 
   ! with a matrix-multiplication based algorithm for a symmetric density matrix
   ! Here we use prev_F as a temporary storage
   do spin = 1, size(density, 3)
      call gemm(amat=overlap, bmat=density(:, :, spin), &
         & cmat=tmpA)
      cache%SP(:, :, spin) = tmpA
      call gemm(amat=tmpA, bmat=overlap, alpha = 0.5_wp, &
         & cmat=cache%prev_F(:, :, spin))
      
      cache%prev_F(:, :, spin) = cache%gmat * cache%prev_F(:, :, spin)
      tmpA = cache%gmat * tmpA
      tmpB = cache%gmat * density(:, :, spin)

      call gemm(amat=overlap, bmat=tmpB, alpha = 0.5_wp, &
         & cmat=tmpA, beta = 1.0_wp)
      call gemm(amat=tmpA, bmat=overlap, alpha = 0.5_wp, &
         & cmat=cache%prev_F(:, :, spin), beta = 1.0_wp)

      tmpB = transpose(cache%prev_F(:, :, spin))
      cache%prev_F(:, :, spin) = -0.25_wp * (cache%prev_F(:, :, spin) + tmpB)

      fock(:, :, spin) = fock(:, :, spin) + 0.5_wp * cache%prev_F(:, :, spin)
   end do 

end subroutine get_KFock


!> Calculate exchange contribution to the gradient
subroutine get_KGrad(self, mol, cache, density, overlap, sigma, gradient)
   !> Instance of the exchange container
   class(mulliken_kfock_matrix), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container with exchange matrix and intermediates
   type(exchange_cache), intent(inout) :: cache
   !> Density matrix
   real(wp),intent(in) :: density(:, :, :)
   !> Overlap matrix
   real(wp),intent(in) :: overlap(:, :)
   !> Molecular gradient of the exchange energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   gradient = 0.0_wp
   sigma = 0.0_wp
end subroutine get_KGrad


!> Evaluate exchange gamma matrix
subroutine get_exchange_matrix(self, mol, cache, gmat)
   !> Instance of the exchange container
   class(mulliken_kfock_matrix), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(exchange_cache), intent(inout) :: cache
   !> Coulomb matrix
   real(wp), contiguous, intent(out) :: gmat(:, :)
   
   integer :: igexp, nao, nsh 
   
   gmat(:, :) = 0.0_wp

   if (any(mol%periodic)) then
      call get_gmat_3d(mol, self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, &
         & self%hubbard, self%gexp, self%frscale, self%omega, self%lrscale, &
         & self%rcut, cache%wsc, cache%alpha, gmat)
   else
      call get_gmat_0d(mol, self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, &
         & self%hubbard, self%gexp, self%frscale, self%omega, self%lrscale, gmat)
   end if

end subroutine get_exchange_matrix


!> Evaluate range separated exchange matrix for finite systems
subroutine get_gmat_0d(mol, nsh_id, nao_sh, ish_at, iao_sh, hubbard, gexp, &
   & frscale, omega, lrscale, gmat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, allocatable :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, allocatable :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, allocatable :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, allocatable :: iao_sh(:)
   !> Hubbard parameter parameter for each shell
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Exponent of exchange kernel
   real(wp), intent(in) :: gexp
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Range separation parameter
   real(wp), intent(in) :: omega
   !> Long-range scaling factor
   real(wp), intent(in) :: lrscale
   !> Exchange matrix
   real(wp), intent(inout) :: gmat(:, :)

   integer :: iat, jat, izp, jzp, is, js, ii, jj, ish, jsh, iao, jao
   real(wp) :: vec(3), r1, r1g, gam, tmp, rsh

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(gmat, mol, nsh_id, nao_sh, ish_at, iao_sh, hubbard) &
   !$omp shared(gexp, frscale, omega, lrscale) &
   !$omp private(iat, izp, is, ii, ish, iao, jat, jzp, js, jj, jsh, jao) &
   !$omp private(gam, vec, r1, r1g, tmp, rsh)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = ish_at(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         js = ish_at(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         r1g = r1**gexp
         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js+jsh)
               gam = hubbard(jsh, ish, jzp, izp) * frscale
               tmp = 1.0_wp/(r1g + gam**(-gexp))**(1.0_wp/gexp)
               rsh = (frscale+lrscale*erf(omega*r1)) * tmp 
               do iao = 1, nao_sh(is + ish)
                  do jao = 1, nao_sh(js + jsh)
                     !$omp atomic
                     gmat(jj+jao, ii+iao) = gmat(jj+jao, ii+iao) + rsh
                     !$omp atomic
                     gmat(ii+iao, jj+jao) = gmat(ii+iao, jj+jao) + rsh
                  end do
               end do
            end do
         end do
      end do
      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         do jsh = 1, ish-1
            jj = iao_sh(is+jsh)
            gam = hubbard(jsh, ish, izp, izp) * frscale
            do iao = 1, nao_sh(is + ish)
               do jao = 1, nao_sh(is + jsh)
                  !$omp atomic
                  gmat(jj+jao, ii+iao) = gmat(jj+jao, ii+iao) + gam
                  !$omp atomic
                  gmat(ii+iao, jj+jao) = gmat(ii+iao, jj+jao) + gam
               end do
            end do
         end do
         ! Diagonal elements
         do iao = 1, nao_sh(is + ish)
            !$omp atomic
            gmat(ii+iao, ii+iao) = gmat(ii+iao, ii+iao) + &
               & hubbard(ish, ish, izp, izp) * frscale
         end do
      end do
   end do 

end subroutine get_gmat_0d


!> Evaluate the coulomb matrix for 3D systems
subroutine get_gmat_3d(mol, nsh_id, nao_sh, ish_at, iao_sh, hubbard, gexp, &
   & frscale, omega, lrscale, rcut, wsc, alpha, gmat)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, allocatable :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, allocatable :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, allocatable :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, allocatable :: iao_sh(:)
   !> Hubbard parameter parameter for each shell
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Exponent of exchange kernel
   real(wp), intent(in) :: gexp
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Range separation parameter
   real(wp), intent(in) :: omega
   !> Long-range scaling factor
   real(wp), intent(in) :: lrscale
   !> Long-range cutoff
   real(wp), intent(in) :: rcut
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence factor
   real(wp), intent(in) :: alpha
   !> Exchange matrix
   real(wp), intent(inout) :: gmat(:, :)

   ! integer :: iat, jat, izp, jzp, img, ii, jj, ish, jsh
   ! real(wp) :: vec(3), gam, wsw, dtmp, rtmp, vol, aval
   ! real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   ! vol = abs(matdet_3x3(mol%lattice))
   ! call get_dir_trans(mol%lattice, alpha, conv, dtrans)
   ! call get_rec_trans(mol%lattice, alpha, vol, conv, rtrans)

   ! !$omp parallel do default(none) schedule(runtime) shared(amat) &
   ! !$omp shared(mol, nshell, offset, hubbard, gexp, wsc, dtrans, rtrans, alpha, vol, rcut) &
   ! !$omp private(iat, izp, jat, jzp, ii, jj, ish, jsh, gam, wsw, vec, dtmp, rtmp, aval)
   ! do iat = 1, mol%nat
   !    izp = mol%id(iat)
   !    ii = offset(iat)
   !    do jat = 1, iat-1
   !       jzp = mol%id(jat)
   !       jj = offset(jat)
   !       wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
   !       do img = 1, wsc%nimg(jat, iat)
   !          vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
   !          call get_amat_rec_3d(vec, vol, alpha, 0.0_wp, rtrans, rtmp)
   !          do ish = 1, nshell(iat)
   !             do jsh = 1, nshell(jat)
   !                gam = hubbard(jsh, ish, jzp, izp)
   !                call get_amat_dir_3d(vec, gam, gexp, rcut, alpha, dtrans, dtmp)
   !                aval = (dtmp + rtmp) * wsw
   !                !$omp atomic
   !                amat(jj+jsh, ii+ish) = amat(jj+jsh, ii+ish) + aval
   !                !$omp atomic
   !                amat(ii+ish, jj+jsh) = amat(ii+ish, jj+jsh) + aval
   !             end do
   !          end do
   !       end do
   !    end do

   !    wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
   !    do img = 1, wsc%nimg(iat, iat)
   !       vec = wsc%trans(:, wsc%tridx(img, iat, iat))
   !       call get_amat_rec_3d(vec, vol, alpha, 0.0_wp, rtrans, rtmp)
   !       rtmp = rtmp - 2 * alpha / sqrtpi
   !       do ish = 1, nshell(iat)
   !          do jsh = 1, ish-1
   !             gam = hubbard(jsh, ish, izp, izp)
   !             call get_amat_dir_3d(vec, gam, gexp, rcut, alpha, dtrans, dtmp)
   !             aval = (dtmp + rtmp + gam) * wsw
   !             !$omp atomic
   !             amat(ii+jsh, ii+ish) = amat(ii+jsh, ii+ish) + aval
   !             !$omp atomic
   !             amat(ii+ish, ii+jsh) = amat(ii+ish, ii+jsh) + aval
   !          end do
   !          gam = hubbard(ish, ish, izp, izp)
   !          call get_amat_dir_3d(vec, gam, gexp, rcut, alpha, dtrans, dtmp)
   !          aval = (dtmp + rtmp + gam) * wsw
   !          !$omp atomic
   !          amat(ii+ish, ii+ish) = amat(ii+ish, ii+ish) + aval
   !       end do
   !    end do

   ! end do

   gmat(:, :) = 0.0_wp

end subroutine get_gmat_3d


!> Evaluate uncontracted derivatives of exchange gamma matrix
subroutine get_exchange_derivs(self, mol, cache, dgdr, dgdL)
   !> Instance of the exchange container
   class(mulliken_kfock_matrix), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(exchange_cache), intent(inout) :: cache
   !> Derivative of interactions with respect to cartesian displacements
   real(wp), contiguous, intent(out) :: dgdr(:, :, :)
   !> Derivative of interactions with respect to strain deformations
   real(wp), contiguous, intent(out) :: dgdL(:, :, :)
end subroutine get_exchange_derivs

end module tblite_mulliken_kfock_matrix

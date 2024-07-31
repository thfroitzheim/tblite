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

!> @file tblite/wavefunction/mulliken.f90
!> Provides Mulliken population analysis

!> Wavefunction analysis via Mulliken populations
module tblite_wavefunction_mulliken
   use iso_fortran_env, only: output_unit
   
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_blas, only: gemm
   use tblite_wavefunction_spin, only : updown_to_magnet
   implicit none
   private

   public :: get_mulliken_shell_charges, get_mulliken_atomic_multipoles
   public :: get_molecular_dipole_moment, get_molecular_quadrupole_moment
   public :: get_mayer_bond_orders, get_mayer_bond_orders_uhf
   public :: get_mulliken_atomic_charges_gradient

contains


subroutine get_mulliken_shell_charges(bas, smat, pmat, n0sh, qsh)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: smat(:, :)
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(in) :: n0sh(:)
   real(wp), intent(out) :: qsh(:, :)

   integer :: iao, jao, spin
   real(wp) :: pao

   qsh(:, :) = 0.0_wp
   !$omp parallel do default(none) collapse(2) schedule(runtime) reduction(+:qsh) &
   !$omp shared(bas, pmat, smat) private(spin, iao, jao, pao)
   do spin = 1, size(pmat, 3)
      do iao = 1, bas%nao
         pao = 0.0_wp
         do jao = 1, bas%nao
            pao = pao + pmat(jao, iao, spin) * smat(jao, iao)
         end do
         qsh(bas%ao2sh(iao), spin) = qsh(bas%ao2sh(iao), spin) - pao
      end do
   end do

   call updown_to_magnet(qsh)
   qsh(:, 1) = qsh(:, 1) + n0sh

end subroutine get_mulliken_shell_charges


!> Evaluate the derivative of the CEH Mulliken charges values with respect to the nuclear coordinates
subroutine get_mulliken_atomic_charges_gradient(bas, mol, smat, pmat, dsmat, dpmatdr, dpmatdL, dqatdr, dqatdL)
   type(basis_type), intent(in) :: bas
   !> Molecular structure data
   type(structure_type), intent(in)  :: mol
   real(wp), intent(in) :: smat(:, :)
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(in) :: dsmat(:, :, :)
   real(wp), intent(in) :: dpmatdr(:, :, :, :)
   real(wp), intent(in) :: dpmatdL(:, :, :, :, :)
   real(wp), intent(out) :: dqatdr(:, :, :, :)
   real(wp), intent(out) :: dqatdL(:, :, :, :)

   integer :: ic, iao, jao, spin, kao, jat, kat, iat
   real(wp) :: dpao, vec(3), r2, r
   real(wp), allocatable :: tmp(:,:), dpsmat(:,:), pdsmat(:,:), tmp1(:,:)

   allocate(tmp(size(dqatdr,2), size(dqatdr,3)), tmp1(bas%nao, bas%nao), dpsmat(bas%nao,bas%nao), pdsmat(bas%nao,bas%nao))
   write(*,*) "in get_mulliken_atomic_charges_gradient"
   dqatdr = 0.0_wp

   ! + ! + ! $omp parallel do default(none) collapse(2) schedule(runtime) reduction(+:qsh) &
   ! + ! + ! $omp shared(bas, pmat, smat) private(spin, iao, jao, pao)
   do spin = 1, size(pmat, 3)
      do ic = 1, 3
         ! write(*,*) "doing dimension", ic
         ! call write_2d_matrix(smat, "smat")
         ! call write_2d_matrix(pmat(:, :, 1), "pmat")      
         ! call write_2d_matrix(dsmat(ic,:,:), "dsmat")
         ! call write_2d_matrix(dpmatdr(ic,:,:,1), "dpmat")
         ! tmp = 0.0_wp
         ! tmp1 = 0.0_wp

         ! tmp1 = smat * dpmatdr(ic,:,:,1)
         ! call write_2d_matrix(tmp1, "s*dp")
         ! tmp1 = 0.0_wp

         ! tmp1 = dsmat(ic,:,:) * pmat(:,:,1)
         ! call write_2d_matrix(tmp1, "ds*p")

         tmp1 = 0.0_wp

         tmp1 = smat * dpmatdr(ic,:,:,1) + dsmat(ic,:,:) * pmat(:,:,1) 
         call write_2d_matrix(tmp1, "s*dp + ds*p")

         do iao = 1, bas%nao
            iat = bas%ao2at(iao)
            dpao = 0.0_wp
            do jao = 1, bas%nao
               jat = bas%ao2at(jao)
               dqatdr(ic,iat,jat,spin) = dqatdr(ic,iat,jat,spin) - tmp1(iao, jao)
               !write(*,*) iat, jat, -tmp1(iao, jao)

            end do 
            do jat = 1, mol%nat
               do kao = 1, bas%nao
                  kat = bas%ao2at(kao)

                  if (jat /= kat .and. kat /= iat ) then
                     !write(*,*) iat, jat, kat, tmp1(iao, kao)
                     dqatdr(ic,iat,jat,spin) = dqatdr(ic,iat,jat,spin) + tmp1(iao, kao)
                  end if
               end do 
            end do 
         end do

      end do
   end do

end subroutine get_mulliken_atomic_charges_gradient



! !> Evaluate the derivative of the CEH Mulliken charges values with respect to the nuclear coordinates
! subroutine get_mulliken_atomic_charges_gradient(bas, mol, smat, pmat, dsmat, dpmatdr, dpmatdL, dqatdr, dqatdL)
!    type(basis_type), intent(in) :: bas
!    !> Molecular structure data
!    type(structure_type), intent(in)  :: mol
!    real(wp), intent(in) :: smat(:, :)
!    real(wp), intent(in) :: pmat(:, :, :)
!    real(wp), intent(in) :: dsmat(:, :, :)
!    real(wp), intent(in) :: dpmatdr(:, :, :, :)
!    real(wp), intent(in) :: dpmatdL(:, :, :, :)
!    real(wp), intent(out) :: dqatdr(:, :, :, :)
!    real(wp), intent(out) :: dqatdL(:, :, :, :)

!    integer :: ic, iao, jao, spin, kao, jat, kat, iat
!    real(wp) :: dpao, vec(3), r2, r
!    real(wp), allocatable :: tmp(:,:), dpsmat(:,:), pdsmat(:,:), tmp1(:,:)

!    allocate(tmp(size(dqatdr,2), size(dqatdr,3)), tmp1(bas%nao, bas%nao), dpsmat(bas%nao,bas%nao), pdsmat(bas%nao,bas%nao))
!    write(*,*) "in get_mulliken_atomic_charges_gradient"
!    dqatdr = 0.0_wp

!    ! + ! + ! $omp parallel do default(none) collapse(2) schedule(runtime) reduction(+:qsh) &
!    ! + ! + ! $omp shared(bas, pmat, smat) private(spin, iao, jao, pao)
!    do spin = 1, size(pmat, 3)
!       do ic = 1, 3
!          write(*,*) "doing dimension", ic
!          call write_2d_matrix(smat, "smat")
!          call write_2d_matrix(pmat(:, :, 1), "pmat")      
!          call write_2d_matrix(dsmat(ic,:,:), "dsmat")
!          call write_2d_matrix(dpmatdr(ic,:,:,1), "dpmat")
!          tmp = 0.0_wp
!          tmp1 = 0.0_wp

!          tmp1 = smat * dpmatdr(ic,:,:,1)
!          call write_2d_matrix(tmp1, "s*dp")
!          tmp1 = 0.0_wp

!          tmp1 = dsmat(ic,:,:) * pmat(:,:,1)
!          call write_2d_matrix(tmp1, "ds*p")

!          tmp1 = 0.0_wp

!          tmp1 = smat * dpmatdr(ic,:,:,1) + dsmat(ic,:,:) * pmat(:,:,1) 
!          call write_2d_matrix(tmp1, "s*dp + ds*p")


!          do iao = 1, bas%nao
!             iat = bas%ao2at(iao)
!             dpao = 0.0_wp
!             do jao = 1, bas%nao ! iat + 1, bas%nao
!                jat = bas%ao2at(jao)
!                dpao = dpao + tmp1(jao, iao)
!                write(*,*) jao, iao, "tmp1", tmp1(jao, iao)

!                !if(iat == jat) then
!                !   dqatdr(ic,jat,iat,spin) = dqatdr(ic,jat,iat,spin) + tmp1(iao, jao)
!                !   dqatdr(ic,jat,iat:,spin) = dqatdr(ic,jat,iat:,spin) - tmp1(iao, jao)
!                !   !& + dsmat(ic, jao, iao) * pmat(iao, jao, spin) + smat(iao, jao) * dpmatdr(ic, jao, iao, spin)
!                !else
!                !   !dqatdr(ic,jat,:,spin) = dqatdr(ic,jat,:,spin) - tmp1(iao, jao) 
!                !   !& - dsmat(ic, jao, iao) * pmat(iao, jao, spin) - smat(iao, jao) * dpmatdr(ic, jao, iao, spin)   
!                !end if 

!                ! dqatdr(ic,iat,jat,spin) = dqatdr(ic,iat,jat,spin) &
!                ! & - dsmat(ic, jao, iao) * pmat(iao, jao, spin) - smat(iao, jao) * dpmatdr(ic, jao, iao, spin)

!                !write(*,*) iat, jat, + dsmat(ic, jao, iao) * pmat(iao, jao, spin), smat(iao, jao) * dpmatdr(ic, jao, iao, spin)
!                if(iat /= jat) then
!                   dqatdr(ic,jat,iat,spin) = dqatdr(ic,jat,iat,spin) + tmp1(jao, iao)
!                end if 
!             end do 
!             dqatdr(ic,iat,iat,spin) = dqatdr(ic,iat,iat,spin) - dpao
!             ! dqatdr(ic,iat,iat,spin) = dqatdr(ic,iat,iat,spin) &
!             ! & - dsmat(ic, iao, iao) * pmat(iao, iao, spin) - smat(iao, iao) * dpmatdr(ic, iao, iao, spin)
!             ! write(*,*) iat, -dsmat(ic, iao, iao) * pmat(iao, iao, spin),  -smat(iao, iao) * dpmatdr(ic, iao, iao, spin)
!             ! do jao = 1, iao-1! iat + 1, bas%nao
!             !    jat = bas%ao2at(jao)
!             !    dqatdr(ic,iat,iat,spin) = dqatdr(ic,iat,iat,spin) &
!             !    & + dsmat(ic, jao, iao) * pmat(iao, jao, spin) + smat(iao, jao) * dpmatdr(ic, jao, iao, spin)
!             !    write(*,*) iat, jat, + dsmat(ic, jao, iao) * pmat(iao, jao, spin), smat(iao, jao) * dpmatdr(ic, jao, iao, spin)
!             ! end do 
!          end do
!          ! write(*,*) "done diagonal"
!          ! do iao = 1, bas%nao
!          !    iat = bas%ao2at(iao)
!          !    do jao = iao + 1, bas%nao
!          !       jat = bas%ao2at(jao)
!          !       dqatdr(ic,iat,jat,spin) = dqatdr(ic,iat,jat,spin) &
!          !       & + dsmat(ic, iao, jao) * pmat(iao, jao, spin) &
!          !       & + smat(iao, jao) * dpmatdr(ic, iao, jao, spin)
!          !       write(*,*) iat, jat, (dsmat(ic, iao, jao)) * pmat(iao, jao, spin) , &
!          !       & smat(iao, jao) * (dpmatdr(ic, iao, jao, spin))

!          !    end do
!          ! end do


!          ! do iao = 1, bas%nao
!          !    dpao = 0.0_wp
!          !    do jao = 1, bas%nao
!          !       !tmp(bas%ao2at(iao), bas%ao2at(jao)) = tmp(bas%ao2at(iao), bas%ao2at(jao)) & 
!          !       tmp1(iao, jao) = &
!          !       & + (dsmat(ic, iao, jao) - dsmat(ic, jao, iao)) * pmat(iao, jao, spin) &
!          !       & + smat(iao, jao) * ((dpmatdr(ic, iao, jao, spin)-dpmatdr(ic, jao, iao, spin)))
!          !    end do
!          ! end do

!          ! do iat = 1, size(dqatdr,2)
!          !    do jat = iat, size(dqatdr,2)
!          !       vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
!          !       r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
!          !       r = sqrt(r2) 
!          !       !do kat = 1,size(dqatdr,3)
!          !          if(iat == jat) then
!          !             !if(kat == jat) then
!          !             dqatdr(ic, iat, jat, 1) = dqatdr(ic, iat, jat, 1) + tmp(iat, jat)
!          !             !else
!          !             !   dqatdr(ic, iat, jat) = dqatdr(ic, iat, jat) - tmp(kat, jat) !* 
!          !             !end if
!          !          else
!          !             !if(kat == jat) then
!          !             dqatdr(ic, iat, jat, 1) = dqatdr(ic, iat, jat, 1) + tmp(iat, jat) !* vec(ic)/r
!          !             dqatdr(ic, jat, iat, 1) = dqatdr(ic, jat, iat, 1) + tmp(jat, iat) !* vec(ic)/r
!          !             !else
!          !             !   dqatdr(ic, iat, jat) = dqatdr(ic, iat, jat) + tmp(kat, jat) !* 
!          !             !end if
!          !          end if 
!          !       !end do
!          !    end do
!          ! end do 
         
!       end do
!    end do

! end subroutine get_mulliken_atomic_charges_gradient


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

subroutine get_mulliken_atomic_multipoles(bas, mpmat, pmat, mpat)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: mpmat(:, :, :)
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(out) :: mpat(:, :, :)

   integer :: iao, jao, spin
   real(wp) :: pao(size(mpmat, 1))

   mpat(:, :, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) reduction(+:mpat) &
   !$omp shared(bas, pmat, mpmat) private(spin, iao, jao, pao)
   do spin = 1, size(pmat, 3)
      do iao = 1, bas%nao
         pao(:) = 0.0_wp
         do jao = 1, bas%nao
            pao(:) = pao + pmat(jao, iao, spin) * mpmat(:, jao, iao)
         end do
         mpat(:, bas%ao2at(iao), spin) = mpat(:, bas%ao2at(iao), spin) - pao
      end do
   end do

   call updown_to_magnet(mpat)

end subroutine get_mulliken_atomic_multipoles


subroutine get_molecular_dipole_moment(mol, qat, dpat, dpmom)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: dpat(:, :)
   real(wp), intent(out) :: dpmom(:)

   integer :: iat

   dpmom(:) = 0.0_wp
   do iat = 1, mol%nat
      dpmom(:) = dpmom + mol%xyz(:, iat) * qat(iat) + dpat(:, iat)
   end do
end subroutine get_molecular_dipole_moment

subroutine get_molecular_quadrupole_moment(mol, qat, dpat, qpat, qpmom)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: dpat(:, :)
   real(wp), intent(in) :: qpat(:, :)
   real(wp), intent(out) :: qpmom(:)

   integer :: iat
   real(wp) :: vec(3), cart(6), tr

   qpmom(:) = 0.0_wp
   do iat = 1, mol%nat
      vec(:) = mol%xyz(:, iat)*qat(iat)
      cart([1, 3, 6]) = mol%xyz(:, iat) * (vec + 2*dpat(:, iat))
      cart(2) = mol%xyz(1, iat) * (vec(2) + dpat(2, iat)) + dpat(1, iat)*mol%xyz(2, iat)
      cart(4) = mol%xyz(1, iat) * (vec(3) + dpat(3, iat)) + dpat(1, iat)*mol%xyz(3, iat)
      cart(5) = mol%xyz(2, iat) * (vec(3) + dpat(3, iat)) + dpat(2, iat)*mol%xyz(3, iat)
      tr = 0.5_wp * (cart(1) + cart(3) + cart(6))
      cart(1) = 1.5_wp * cart(1) - tr
      cart(2) = 3.0_wp * cart(2)
      cart(3) = 1.5_wp * cart(3) - tr
      cart(4) = 3.0_wp * cart(4)
      cart(5) = 3.0_wp * cart(5)
      cart(6) = 1.5_wp * cart(6) - tr
      qpmom(:) = qpmom(:) + qpat(:, iat) + cart
   end do
end subroutine get_molecular_quadrupole_moment


!> Evaluate Wiberg/Mayer bond orders
subroutine get_mayer_bond_orders(bas, smat, pmat, mbo)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(in) :: smat(:, :)
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :, :)
   !> Wiberg/Mayer bond orders
   real(wp), intent(out) :: mbo(:, :, :)

   integer :: iao, jao, iat, jat, spin
   real(wp) :: pao
   real(wp), allocatable :: psmat(:, :)

   allocate(psmat(bas%nao, bas%nao))

   mbo(:, :, :) = 0.0_wp
   do spin = 1, size(pmat, 3)
      call gemm(pmat(:, :, spin), smat, psmat)
      !$omp parallel do default(none) collapse(2) &
      !$omp shared(bas, psmat, mbo, spin) private(iao, jao, iat, jat, pao)
      do iao = 1, bas%nao
         do jao = 1, bas%nao
            iat = bas%ao2at(iao)
            jat = bas%ao2at(jao)
            pao = merge(psmat(iao, jao) * psmat(jao, iao), 0.0_wp, iat /= jat)
            !$omp atomic
            mbo(jat, iat, spin) = mbo(jat, iat, spin) + pao
         end do
      end do
   end do

   call updown_to_magnet(mbo)
end subroutine get_mayer_bond_orders

!> Evaluate Wiberg/Mayer bond orders
subroutine get_mayer_bond_orders_uhf(bas, smat, pmat, mbo)
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(in) :: smat(:, :)
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :, :)
   !> Wiberg/Mayer bond orders
   real(wp), intent(out) :: mbo(:, :, :)

   integer :: iao, jao, iat, jat, spin
   real(wp) :: pao
   real(wp), allocatable :: psmat(:, :)

   allocate(psmat(bas%nao, bas%nao))

   mbo(:, :, :) = 0.0_wp
   do spin = 1, size(pmat, 3)
      call gemm(pmat(:, :, spin), smat, psmat)
      !$omp parallel do default(none) collapse(2) &
      !$omp shared(bas, psmat, mbo, spin) private(iao, jao, iat, jat, pao)
      do iao = 1, bas%nao
         do jao = 1, bas%nao
            iat = bas%ao2at(iao)
            jat = bas%ao2at(jao)
            pao = merge(psmat(iao, jao) * psmat(jao, iao), 0.0_wp, iat /= jat)
            !$omp atomic
            mbo(jat, iat, 1) = mbo(jat, iat, 1) + pao
         end do
      end do
   end do

   call updown_to_magnet(mbo)
end subroutine get_mayer_bond_orders_uhf


end module tblite_wavefunction_mulliken

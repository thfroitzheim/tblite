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

module test_qvszp
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_basis_qvszp, only : qvszp_basis_type, qvszp_cgto_type, new_basis, add_qvszp_cgtos
   use tblite_integral_overlap, only : overlap_cgto

   implicit none
   private

   public :: collect_qvszp

   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains

   !> Collect all exported unit tests
   subroutine collect_qvszp(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("norm-qvszp-h2", test_norm_qvszp_cgto_h2), &
         new_unittest("norm-qvszp-lih", test_norm_qvszp_cgto_lih), &
         new_unittest("norm-qvszp-s2", test_norm_qvszp_cgto_s2), &
         new_unittest("norm-qvszp-sih4", test_norm_qvszp_cgto_sih4) &
         & ]

   end subroutine collect_qvszp


   subroutine test_norm_qvszp_cgto(error, mol)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(structure_type), intent(in) :: mol

      type(qvszp_cgto_type), allocatable :: cgto(:, :)
      integer, allocatable :: nsh_id(:)

      real(wp), parameter :: vec(3) = 0.0_wp, r2 = 0.0_wp
      real(wp) :: overlap(9, 9)

      integer :: iat, isp, ish, iao, l

      ! Setup qvszp basis
      call add_qvszp_cgtos(cgto, mol, nsh_id)

      do iat = 1, mol%nat
         isp = mol%id(iat)
         do ish = 1, nsh_id(isp)
            overlap = 0.0_wp
            call overlap_cgto(cgto(ish,iat), cgto(ish,iat), r2, vec, 100.0_wp, overlap)
            l = cgto(ish,iat)%ang
            do iao = 1, l*(l+1)
               call check(error, overlap(iao, iao), 1.0_wp, thr=thr)
               if (allocated(error)) then
                  print*, (overlap(iao, iao) - 1.0_wp)
               end if
            end do
         end do
      end do

   end subroutine test_norm_qvszp_cgto

   subroutine test_norm_qvszp_cgto_h2(error)
   
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol
   
      call get_structure(mol, "MB16-43", "H2")
      call test_norm_qvszp_cgto(error, mol)
   
   end subroutine test_norm_qvszp_cgto_h2
   
   subroutine test_norm_qvszp_cgto_lih(error)
   
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol
   
      call get_structure(mol, "MB16-43", "LiH")
      call test_norm_qvszp_cgto(error, mol)
   
   end subroutine test_norm_qvszp_cgto_lih
   
   subroutine test_norm_qvszp_cgto_s2(error)
   
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol
   
      call get_structure(mol, "MB16-43", "S2")
      call test_norm_qvszp_cgto(error, mol)
   
   end subroutine test_norm_qvszp_cgto_s2

   subroutine test_norm_qvszp_cgto_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_norm_qvszp_cgto(error, mol)

   end subroutine test_norm_qvszp_cgto_sih4


end module test_qvszp

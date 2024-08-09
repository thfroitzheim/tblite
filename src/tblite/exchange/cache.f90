
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

!> @file tblite/exchange/cache.f90
!> Provides a cache specific for all exchange interactions

!> Data container for mutable data in exchange calculations
module tblite_exchange_cache
   use mctc_env, only : wp, dp
   use mctc_io, only : structure_type
   use tblite_coulomb_ewald, only : get_alpha
   use tblite_wignerseitz, only : wignerseitz_cell, new_wignerseitz_cell
   implicit none
   private

   public :: exchange_cache


   type :: exchange_cache
      real(wp) :: alpha
      type(wignerseitz_cell) :: wsc
      !> Exchange matrix
      real(wp), allocatable :: gmat(:, :)
      !> Overlap times density matrix intermediate
      real(wp), allocatable :: SP(:, :, :)

      !> Note whether previous Fock matrix is available 
      logical :: has_prev_F = .false. 
      !> Reference density matrix for incremental build
      real(dp), allocatable :: ref_P(:, :, :)
      !> Reference Fock matrix for incremental build
      real(dp), allocatable :: ref_F(:, :, :)
      !> Previously calculated Fock matrix contribution 
      real(dp), allocatable :: prev_F(:, :, :)
   contains
      procedure :: update
   end type exchange_cache


contains


subroutine update(self, mol)
   !> Instance of the electrostatic container
   class(exchange_cache), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   if (any(mol%periodic)) then
      call new_wignerseitz_cell(self%wsc, mol)
      call get_alpha(mol%lattice, self%alpha, .false.)
   end if

end subroutine update

end module tblite_exchange_cache

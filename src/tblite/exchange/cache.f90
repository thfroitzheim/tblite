
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
      real(dp), allocatable :: curr_D(:, :)
      real(dp), allocatable :: ref_D(:, :)
      real(dp), allocatable :: prev_F(:, :)
      real(dp), allocatable :: gamma_(:,:)
   contains
      procedure :: update
   end type exchange_cache


contains


subroutine update(self, mol)
   !> Instance of the electrostatic container
   class(exchange_cache), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol


end subroutine update

end module tblite_exchange_cache

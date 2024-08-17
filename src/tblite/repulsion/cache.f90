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

!> @file tblite/reulsion/cache.f90
!> Provides a cache for repulsion interactions

!> Reusable data container for repulsion related calculations
module tblite_repulsion_cache
   use mctc_env, only : wp
   implicit none
   private
   public :: repulsion_cache

   type :: repulsion_cache
      real(wp), allocatable :: scaled_zeff(:, :)
      real(wp), allocatable :: scaled_dzeffdr(:, :, :, :)
      real(wp), allocatable :: scaled_dzeffdL(:, :, :, :)
   end type repulsion_cache

end module tblite_repulsion_cache

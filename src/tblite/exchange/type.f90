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

!> @file tblite/exchange/type.f90
!> Provides a base class for all exchange interactions
module tblite_exchange_type
   use tblite_container_type, only : container_type
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_container_cache, only : container_cache
   use mctc_env, only : wp

   implicit none
   private

   !> General base class for Coulombic interactions
   type, public, extends(container_type), abstract :: exchange_type
   contains
      !> Evaluate charge dependent potential shift from the interaction
      procedure :: get_potential_w_overlap
      !> Evaluate charge dependent potential shift from the interaction
      procedure :: get_gradient_w_overlap
   end type exchange_type
contains

   !> Evaluate charge dependent potential shift from the interaction
subroutine get_potential_w_overlap(self, mol, cache, wfn, pot, overlap)
   !> Instance of the interaction container
   class(exchange_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot
   !> Overlap integral matrix
   real(wp), intent(in) :: overlap(:,:)
end subroutine get_potential_w_overlap

subroutine get_gradient_w_overlap(self, mol, cache, wfn, gradient, sigma, overlap)
   !> Instance of the interaction container
   class(exchange_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
    !> Molecular gradient of the exchange energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)
   !> Overlap integral matrix
   real(wp), intent(in) :: overlap(:,:)
end subroutine get_gradient_w_overlap

end module tblite_exchange_type

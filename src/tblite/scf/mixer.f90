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

!> @dir tblite/scf/mixer
!> Routines for implementing electronic mixing

!> @file tblite/scf/mixer.f90
!> Proxy module for electronic mixing routines

!> Provides an electronic mixer implementation
module tblite_scf_mixer
   use mctc_env, only : wp
   use tblite_scf_mixer_broyden, only : broyden_mixer, broyden_input, new_broyden
   use tblite_scf_mixer_type, only : mixer_type
   implicit none
   private

   public :: mixer_type, new_mixer


   !> Input for selecting electronic mixer
   type, public :: mixer_input
      !> Input for Broyden mixer
      type(broyden_input), allocatable :: broyden
   end type mixer_input

contains

!> Create a new instance of the mixer
subroutine new_mixer(self, memory, ndim, damp, start)
   !> Instance of the mixer on exit
   class(mixer_type), allocatable, intent(out) :: self
   !> Dimension of the subspace used for extrapolation
   integer, intent(in) :: memory
   !> Size of the extrapolated quantity
   integer, intent(in) :: ndim
   !> Damping parameter
   real(wp), intent(in) :: damp
   !> Iteration to start mixing
   integer, intent(in) :: start

   block
      type(broyden_mixer), allocatable :: mixer
      allocate(mixer)
      call new_broyden(mixer, ndim, broyden_input(memory, damp, start))
      call move_alloc(mixer, self)
   end block
end subroutine new_mixer

end module tblite_scf_mixer

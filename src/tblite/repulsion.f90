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

!> @dir tblite/repulsion
!> Contains the repulsive interaction implementations.

!> @file tblite/repulsion.f90
!> Provides a repulsion container base class

!> Proxy module for repulsion interactions.
module tblite_repulsion
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_repulsion_gfn, only : gfn_repulsion, new_repulsion_gfn
   use tblite_repulsion_gp3, only : gp3_repulsion, new_repulsion_gp3
   use tblite_repulsion_type, only : repulsion_type
   implicit none
   private

   public :: repulsion_type, new_repulsion

contains

!> Create a new repulsion container
subroutine new_repulsion(self, mol, rep_type, zeff, alpha, rexp, &
   & kexp, kexp_light, kcn, kq, rcov_rep, rcov_cn, exp_cn, cutoff) 
   !> Instance of the repulsion container
   class(repulsion_type), allocatable, intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Repulsion type
   character(len=*), intent(in) :: rep_type
   !> Exponent of the repulsion polynomial
   real(wp), intent(in) :: rexp
   !> Effective nuclear charge
   real(wp), intent(in) :: zeff(:)
   !> Exponent for the damping of the repulsion
   real(wp), intent(in) :: alpha(:)
   !> Distance exponent for repulsion damping
   real(wp), intent(in), optional :: kexp
   !> Distance exponent for repulsion damping of light atom pairs
   real(wp), intent(in), optional :: kexp_light
   !> Coefficient for CN dependence
   real(wp), intent(in), optional :: kcn(:)
   !> Coefficient for atomic charge dependence
   real(wp), intent(in), optional :: kq(:)
   !> Fitted covalent atomic radii for repulsion
   real(wp), intent(in), optional :: rcov_rep(:)
   !> Fitted covalent atomic radii for repulsion CN
   real(wp), intent(in), optional :: rcov_cn(:)
   !> Exponent of the repulsion CN
   real(wp), intent(in), optional :: exp_cn
   !> Real-space cutoff
   real(wp), intent(in), optional :: cutoff

   select case(rep_type)
   case("gfn")
      block
         type(gfn_repulsion), allocatable :: tmp
         allocate(tmp)
         call new_repulsion_gfn(tmp, mol, alpha, zeff, kexp, kexp_light, rexp, cutoff)
         call move_alloc(tmp, self)
      end block
   case("gp3")
      block
         type(gp3_repulsion), allocatable :: tmp
         allocate(tmp)
         call new_repulsion_gp3(tmp, mol, alpha, zeff, kcn, kq, rexp,&
            & rcov_rep, rcov_cn, exp_cn, cutoff)
         call move_alloc(tmp, self)
      end block
   end select

end subroutine new_repulsion

end module tblite_repulsion

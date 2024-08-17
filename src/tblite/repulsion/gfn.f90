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

!> @file tblite/repulsion/gfn.f90
!> Provides a screened Coulomb repulsion interaction

!> Classical repulsion interaction as used with the GFNn-xTB Hamiltonian
module tblite_repulsion_gfn
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_container, only : container_cache
   use tblite_repulsion_cache, only : repulsion_cache
   use tblite_repulsion_type, only : repulsion_type, view, taint
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_repulsion_gfn


   !> Container to evaluate classical repulsion interactions for the GFN Hamiltonian
   type, public, extends(repulsion_type) :: gfn_repulsion
      !> Exponent for the damping of the repulsion
      real(wp), allocatable :: alpha(:, :)
      !> Effective nuclear charge
      real(wp), allocatable :: zeff(:, :)
      !> Scaling of the repulsion exponents
      real(wp), allocatable :: kexp(:, :)
   contains
      !> Update repulsion cache
      procedure :: update
      !> Evaluates the damping for classical repulsion
      procedure :: get_damping
      !> Evaluates the derivative of the damping for classical repulsion
      procedure :: get_damping_derivs
   end type gfn_repulsion

   character(len=*), parameter :: label = "screened Coulomb repulsion"

contains


subroutine new_repulsion_gfn(self, mol, alpha, zeff, kexp, kexp_light, rexp, cutoff)
   !> Instance of the repulsion container
   type(gfn_repulsion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Exponent for the damping of the repulsion
   real(wp), intent(in) :: alpha(:)
   !> Effective nuclear charge
   real(wp), intent(in) :: zeff(:)
   !> Distance exponent for repulsion damping
   real(wp), intent(in) :: kexp
   !> Distance exponent for repulsion damping of light atom pairs
   real(wp), intent(in) :: kexp_light
   !> Exponent of the repulsion polynomial
   real(wp), intent(in) :: rexp
   !> Real-space cutoff
   real(wp), intent(in), optional :: cutoff

   integer :: isp, izp, jsp, jzp

   self%label = label
   allocate(self%alpha(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%alpha(jsp, isp) = sqrt(alpha(isp)*alpha(jsp))
      end do
   end do

   allocate(self%zeff(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%zeff(jsp, isp) = zeff(isp)*zeff(jsp)
      end do
   end do

   allocate(self%kexp(mol%nid, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, mol%nid
         jzp = mol%num(jsp)
         self%kexp(jsp, isp) = merge(kexp, kexp_light, izp > 2 .or. jzp > 2)
      end do
   end do

   allocate(self%rexp(mol%nid, mol%nid))
   self%rexp(:, :) = rexp
   if (present(cutoff)) self%cutoff = cutoff

end subroutine new_repulsion_gfn


!> Update repulsion cache
subroutine update(self, mol, cache, wfn)
   !> Instance of the classical repulsion
   class(gfn_repulsion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Optional auxiliary wavefunction data
   type(wavefunction_type), intent(in), optional :: wfn

   real(wp), allocatable :: lattr(:, :)
   type(repulsion_cache), pointer :: ptr
   integer :: iat, izp, jat, jzp

   call taint(cache, ptr)

   if (.not.allocated(ptr%scaled_zeff)) allocate(ptr%scaled_zeff(mol%nat, mol%nat))
   if (.not.allocated(ptr%scaled_dzeffdr)) then 
      allocate(ptr%scaled_dzeffdr(3, mol%nat, mol%nat, mol%nat), source=0.0_wp)
   end if
   if (.not.allocated(ptr%scaled_dzeffdL)) then 
      allocate(ptr%scaled_dzeffdL(3, 3, mol%nat, mol%nat), source=0.0_wp)
   end if 

   ! Distribute the effective nuclear charges to be atom resolved
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         ptr%scaled_zeff(jat, iat) = self%zeff(jzp, izp)
      end do
   end do 

end subroutine update


!> Evaluates the damping for GFN repulsion
elemental function get_damping(self, izp, jzp, r) result(damp)
   !> Instance of the repulsion container
   class(gfn_repulsion), intent(in) :: self
   !> Atom i index
   integer, intent(in)  :: izp
   !> Atom j index
   integer, intent(in)  :: jzp
   !> Current distance
   real(wp), intent(in) :: r

   real(wp) :: damp, r1k

   r1k = r**self%kexp(jzp, izp)
   damp = exp(-self%alpha(jzp, izp)*r1k)

end function get_damping

!> Evaluates the derivative of the damping for GFN repulsion
elemental function get_damping_derivs(self, izp, jzp, r) result(ddamp)
   !> Instance of the repulsion container
   class(gfn_repulsion), intent(in) :: self
   !> Atom i index
   integer, intent(in)  :: izp
   !> Atom j index
   integer, intent(in)  :: jzp
   !> Current distance
   real(wp), intent(in) :: r

   real(wp) :: ddamp, prefactor, r1kmin

   r1kmin = r**(self%kexp(jzp, izp) - 1)
   prefactor = - self%alpha(jzp, izp) * r1kmin * self%kexp(jzp, izp)
   ddamp = prefactor * exp(-self%alpha(jzp, izp) * r1kmin * r)

end function get_damping_derivs

end module tblite_repulsion_gfn

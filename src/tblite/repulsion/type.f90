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

!> @file tblite/repulsion/type.f90
!> Provides a general implementation of the screened Coulomb repulsion interaction

!> Classical repulsion interaction as used with the xTB Hamiltonian
module tblite_repulsion_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_classical_type, only : classical_type
   use tblite_container, only : container_cache
   use tblite_repulsion_cache, only : repulsion_cache
   use tblite_cutoff, only : get_lattice_points
   implicit none
   private

   public :: view, taint

   !> Container to evaluate classical repulsion interactions for the xTB Hamiltonian
   type, public, extends(classical_type), abstract :: repulsion_type
      !> Exponent of the repulsion polynomial
      real(wp), allocatable :: rexp(:, :)
      !> Real-space cutoff
      real(wp) :: cutoff = 25.0_wp
   contains
      !> Evaluate non-selfconsistent repulsion energy and gradient
      procedure :: get_engrad
      !> Evaluates the damping for classical repulsion
      procedure(get_damping), deferred :: get_damping
      !> Evaluates the derivative of the damping for classical repulsion
      procedure(get_damping_derivs), deferred :: get_damping_derivs
   end type repulsion_type

   abstract interface

      !> Abstract damping function
      elemental function get_damping(self, izp, jzp, r) result(damp)
         import :: repulsion_type, wp
         !> Instance of the repulsion container
         class(repulsion_type), intent(in) :: self
         !> Atom i index
         integer, intent(in)  :: izp
         !> Atom j index
         integer, intent(in)  :: jzp
         !> Current distance
         real(wp), intent(in) :: r

         real(wp) :: damp
      end function get_damping

      !> Abstract derivative of damping function
      elemental function get_damping_derivs(self, izp, jzp, r) result(ddamp)
         import :: repulsion_type, wp
         !> Instance of the repulsion container
         class(repulsion_type), intent(in) :: self
         !> Atom i index
         integer, intent(in)  :: izp
         !> Atom j index
         integer, intent(in)  :: jzp
         !> Current distance
         real(wp), intent(in) :: r

         real(wp) :: ddamp
      end function get_damping_derivs

   end interface

contains


!> Evaluate classical repulsion interaction for energy and derivatives
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the repulsion container
   class(repulsion_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Repulsion energy
   real(wp), intent(inout) :: energies(:)
   !> Molecular gradient of the repulsion energy
   real(wp), contiguous, intent(inout), optional :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), contiguous, intent(inout), optional :: sigma(:, :)

   real(wp), allocatable :: trans(:, :)
   type(repulsion_cache), pointer :: ptr
   
   call view(cache, ptr)

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff, trans)

   if (present(gradient) .and. present(sigma)) then
      call get_repulsion_derivs(self, mol, trans, ptr%scaled_zeff, &
         & ptr%scaled_dzeffdr, ptr%scaled_dzeffdL, energies, gradient, sigma)
   else
      call get_repulsion_energy(self, mol, trans, ptr%scaled_zeff, &
         & energies)
   end if

end subroutine get_engrad


!> Evaluates the repulsion energy
subroutine get_repulsion_energy(self, mol, trans, zeff, energies)
   !> Instance of the repulsion container
   class(repulsion_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Effective nuclear charges for all element pairs
   real(wp), intent(in) :: zeff(:, :)
   !> Repulsion energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r1, r2, rij(3), r1r, damp, cutoff2, dE

   cutoff2 = self%cutoff**2

   !$omp parallel do default(none) schedule(runtime) reduction(+:energies) &
   !$omp shared(self, mol, trans, cutoff2, zeff) &
   !$omp private(iat, jat, izp, jzp, itr, r1, r2, rij, r1r, damp, dE)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)
            
            ! Damping function
            damp = self%get_damping(izp, jzp, r1)
            ! Coulomb type (1/r^n) interaction 
            r1r = r1**self%rexp(jzp, izp)
            dE = zeff(jat, iat) * damp/r1r
            energies(iat) = energies(iat) + 0.5_wp * dE
            if (iat /= jat) then
               energies(jat) = energies(jat) + 0.5_wp * dE
            end if
         end do
      end do
   end do

end subroutine get_repulsion_energy

!> Evaluates the repulsion energy and gradient
subroutine get_repulsion_derivs(self, mol, trans, zeff, dzeffdr, dzeffdL, &
   & energies, gradient, sigma)
   !> Instance of the repulsion container
   class(repulsion_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Effective nuclear charges for all element pairs
   real(wp), intent(in) :: zeff(:, :)
   !> Derivative of the effective nuclear charge w.r.t. coordinates
   real(wp), intent(in) :: dzeffdr(:, :, :, :)
   !> Derivative of the effective nuclear charge w.r.t. lattice vectors
   real(wp), intent(in) :: dzeffdL(:, :, :, :)
   !> Repulsion energy
   real(wp), intent(inout) :: energies(:)
   !> Molecular gradient of the repulsion energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r1, r2, rij(3), r1k, r1r, damp, ddamp(3), cutoff2
   real(wp) :: dE, dG(3), dS(3, 3)

   cutoff2 = self%cutoff**2

   !$omp parallel do default(none) schedule(runtime) reduction(+:energies, gradient, sigma) &
   !$omp shared(self, mol, trans, cutoff2, zeff, dzeffdr, dzeffdL) &
   !$omp private(iat, jat, izp, jzp, itr, r1, r2, rij, r1r, damp, ddamp, dE, dG, dS)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            ! Damping function
            damp = self%get_damping(izp, jzp, r1)
            ddamp = self%get_damping_derivs(izp, jzp, r1)
            ! Coulomb type (1/r^n) interaction 
            r1r = r1**self%rexp(jzp, izp)
            dE = zeff(jat, iat) * damp/r1r
            ! Explicit derivative of coulomb type (1/r^n) interaction and damping
            dG = (zeff(jat, iat) * ddamp / r1r) * rij/r1
            dG = dG - (self%rexp(jzp, izp) * dE) * rij/r2
            dS = spread(dG, 1, 3) * spread(rij, 2, 3)
            energies(iat) = energies(iat) + 0.5_wp * dE
            if (iat /= jat) then
               energies(jat) = energies(jat) + 0.5_wp * dE
               gradient(:, iat) = gradient(:, iat) + dG
               gradient(:, jat) = gradient(:, jat) - dG
               sigma(:, :) = sigma + dS
            else
               sigma(:, :) = sigma + 0.5_wp * dS
            end if
            ! Implicit derivative of the effective nuclear charge dependence
            gradient(:, :) = gradient(:, :) + dzeffdr(:, :, jat, iat) * damp/r1r
            sigma(:, :) = sigma(:, :) + dzeffdL(:, :, jat, iat) * damp/r1r
         end do
      end do
   end do

end subroutine get_repulsion_derivs


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(repulsion_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(repulsion_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to container cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(repulsion_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(repulsion_cache)
      ptr => target
   end select
end subroutine view

end module tblite_repulsion_type

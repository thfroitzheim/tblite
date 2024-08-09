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

!> @file tblite/exchange/mulliken_kfock.f90
!> Provides a base class for calculating Mulliken approximated Fock exchange. 
!> Child classes contain different algorithms. 

!> General Mulliken approximated Fock exchange
module tblite_mulliken_kfock_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_blas, only : gemm
   use tblite_exchange_type, only : exchange_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_scf_potential, only : potential_type
   use tblite_container_cache, only : container_cache
   use tblite_exchange_cache, only : exchange_cache
   use tblite_basis_type, only : basis_type
   implicit none
   private


   !> General Mulliken approximated Fock exchange
   type, public, extends(exchange_type), abstract :: mulliken_kfock_type
      !> Whether the Hubbard parameters are shell-dependent
      logical :: shell_resolved_hubbard
      !> Full-range admixture of exchange
      real(wp) :: frscale
      !> Range separation parameter 
      real(wp) :: omega
      !> Scaling of the range separated exchange
      real(wp) :: lrscale
      !> Number of AOS
      integer :: nao
      !> Number of shells
      integer :: nsh
      !> Number of shells for each species
      integer, allocatable :: nsh_id(:)
      !> Number of spherical atomic orbitals for each shell
      integer, allocatable :: nao_sh(:)
      !> Index offset for each atom in the shell space
      integer, allocatable :: ish_at(:)
      !> Index offset for each shell in the atomic orbital space
      integer, allocatable :: iao_sh(:)
      !> Mapping from shells to the respective atom
      integer, allocatable :: sh2at(:)
   contains
      !> Update container cache
      procedure :: update
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Evaluate selfconsistent energy of the interaction
      procedure :: get_energy
      ! !> Evaluate density dependent potential
      ! procedure :: get_potential_w_overlap
      !> Evaluate gradient contributions from the selfconsistent interaction
      procedure :: get_gradient_w_overlap
      !> Information on container
      procedure :: info
      !> Evaluate exchange gamma matrix
      procedure(get_exchange_matrix), deferred :: get_exchange_matrix
      !> Evaluate uncontracted derivatives of exchange gamma matrix
      procedure(get_exchange_derivs), deferred :: get_exchange_derivs
      !> Calculate exchange contribution to the Fock matrix
      procedure(get_KFock), deferred :: get_KFock
      !> Calculate exchange contribution to the gradient
      procedure(get_KGrad), deferred :: get_KGrad
   end type mulliken_kfock_type

   abstract interface
      !> Evaluate exchange gamma matrix
      subroutine get_exchange_matrix(self, mol, cache, gmat)
         import :: mulliken_kfock_type, structure_type, exchange_cache, wp
         !> Instance of the exchange container
         class(mulliken_kfock_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container
         type(exchange_cache), intent(inout) :: cache
         !> Exchange matrix
         real(wp), contiguous, intent(out) :: gmat(:, :)
      end subroutine get_exchange_matrix

      !> Evaluate uncontracted derivatives of exchange gamma matrix
      subroutine get_exchange_derivs(self, mol, cache, dgdr, dgdL)
         import :: mulliken_kfock_type, structure_type, exchange_cache, wp
         !> Instance of the exchange container
         class(mulliken_kfock_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container
         type(exchange_cache), intent(inout) :: cache
         !> Derivative of interactions with respect to cartesian displacements
         real(wp), contiguous, intent(out) :: dgdr(:, :, :)
         !> Derivative of interactions with respect to strain deformations
         real(wp), contiguous, intent(out) :: dgdL(:, :, :)
      end subroutine get_exchange_derivs

      !> Calculate exchange contribution to the Fock matrix
      subroutine get_KFock(self, mol, cache, density, overlap, fock)
         import :: mulliken_kfock_type, structure_type, exchange_cache, wp
         !> Instance of the exchange container
         class(mulliken_kfock_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container with exchange matrix and intermediates
         type(exchange_cache), intent(inout) :: cache
         !> Density matrix
         real(wp),intent(in) :: density(:, :, :)
         !> Overlap matrix
         real(wp),intent(in) :: overlap(:, :)
         !> Fock Matrix contribution
         real(wp),intent(inout) :: fock(:, :, :)      
      end subroutine get_KFock

      !> Calculate exchange contribution to the gradient
      subroutine get_KGrad(self, mol, cache, density, overlap, sigma, gradient)
         import :: mulliken_kfock_type, structure_type, exchange_cache, wp
         !> Instance of the exchange container
         class(mulliken_kfock_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container with exchange matrix and intermediates
         type(exchange_cache), intent(inout) :: cache
         !> Density matrix
         real(wp),intent(in) :: density(:, :, :)
         !> Overlap matrix
         real(wp),intent(in) :: overlap(:, :)
         !> Molecular gradient of the exchange energy
         real(wp), contiguous, intent(inout) :: gradient(:, :)
         !> Strain derivatives of the exchange energy
         real(wp), contiguous, intent(inout) :: sigma(:, :)
      end subroutine get_KGrad

   end interface

contains


!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the multipole container
   class(mulliken_kfock_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(exchange_cache), pointer :: ptr

   call taint(cache, ptr)
   call ptr%update(mol)

   if (.not.allocated(ptr%gmat)) then
      allocate(ptr%gmat(self%nao, self%nao))
   end if
   call self%get_exchange_matrix(mol, ptr, ptr%gmat)

end subroutine update


!> Evaluate selfconsistent energy of the interaction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the exchange container
   class(mulliken_kfock_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(exchange_cache), pointer :: ptr
   integer :: iao, jao, spin

   call view(cache, ptr)

   !$omp parallel do collapse(3) schedule(runtime) default(none) &
   !$omp reduction(+:energies) shared(ptr, wfn) private(spin, iao, jao)
   do spin = 1, size(wfn%density, 3)
      do iao = 1, size(wfn%density, 2)
         do jao = 1, size(wfn%density, 1)
            energies(iao) = energies(iao) &
               & + 0.5_wp * ptr%prev_F(jao, iao, spin) * wfn%density(jao, iao, spin)
         end do
      end do
   end do
   
end subroutine get_energy


!> Evaluate density dependent potential 
subroutine get_potential_w_overlap(self, mol, cache, wfn, pot, overlap)
   !> Instance of the exchange container
   class(mulliken_kfock_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential shifts
   type(potential_type), intent(inout) :: pot
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Overlap matrix
   real(wp), intent(in) :: overlap(:,:)

   type(exchange_cache), pointer :: ptr

   call view(cache, ptr)

   if (.not.allocated(ptr%prev_F)) then 
      allocate(ptr%prev_F(self%nao, self%nao, wfn%nspin), source = 0.0_wp)
   end if
   if (.not.allocated(pot%kao)) then 
      allocate(pot%kao(self%nao, self%nao, wfn%nspin), source = 0.0_wp)
   end if

   call self%get_KFock(mol, ptr, wfn%density, overlap, pot%kao)

end subroutine get_potential_w_overlap


!> Evaluate gradient contributions from the selfconsistent interaction
subroutine get_gradient_w_overlap(self, mol, cache, wfn, gradient, sigma, overlap)
   !> Instance of the exchange container
   class(mulliken_kfock_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Molecular gradient of the exchange energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)
   !> Overlap
   real(wp), intent(in) :: overlap(:,:)

   type(exchange_cache), pointer :: ptr
   integer :: i

   call view(cache, ptr)

   !call KGradSymSQM(0, self%nao, mol%nat, self%nsh, self%aonum, self%sh2at, ptr%gamma_, &
   !   & wfn%density(:, :, 1), overlap, sigma, gradient)

end subroutine get_gradient_w_overlap


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, orbital_resolved, not_used
   !> Instance of the electrostatic container
   class(mulliken_kfock_type), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=not_used ,dipole=not_used, quadrupole=not_used, density=orbital_resolved)
end function variable_info


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(exchange_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(exchange_cache), allocatable :: tmp
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
   type(exchange_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(exchange_cache)
      ptr => target
   end select
end subroutine view


!> Information on container
pure function info(self, verbosity, indent) result(str)
   use tblite_output_format, only : format_string
   !> Instance of the interaction container
   class(mulliken_kfock_type), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   character(len=*), parameter :: nl = new_line('a'), marker = " * "

   if (allocated(self%label)) then
      str = self%label
   else
      str = "Mulliken approximated Fock exchange"
   end if

   if (verbosity > 1) then
      str = str // nl // indent // "Using the shell resolved chemical hardness"
      if (self%omega.ne.0.0_wp .and. self%lrscale.ne.0.0_wp) then
         str= str // nl // indent // "Range separted exchange is used:" // &
            & nl // indent // marker // "Full-range scale: " // format_string(self%frscale,'(f5.2)') // &
            & nl // indent // marker // "Omega: " // format_string(self%omega,'(f5.2)') // &
            & nl // indent // marker // "Long-range scale: " // format_string(self%lrscale,'(f5.2)')
      else
         str = str // nl // indent // "Full range exchange is used" // &
            & nl // indent // marker // "Full-range scale: " // format_string(self%frscale,'(f5.2)')
      end if
   end if 
end function info

end module tblite_mulliken_kfock_type

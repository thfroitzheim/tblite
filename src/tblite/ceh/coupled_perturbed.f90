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

!> @file tblite/ceh/singlepoint.f90
!> Provides main entry point for performing a CEH calculation
!> follows in close analogy the xtb singlepoint

!> Implementation of the coupled-perturbed equations for the
!> density matrix and Mulliken charge derivatives
module tblite_ceh_coupled_perturbed
   use mctc_env, only : error_type, wp
   use mctc_io, only: structure_type
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_cutoff, only : get_lattice_points
   use tblite_basis_type, only : basis_type, get_cutoff
   use tblite_context, only : context_type
   use tblite_output_format, only: format_string
   use tblite_wavefunction, only : wavefunction_type, &
   & get_alpha_beta_occupation
   use tblite_wavefunction_mulliken, only: get_mulliken_shell_charges, &
   & get_mulliken_atomic_multipoles
   use tblite_blas, only : gemv
   implicit none
   private

   public :: c

   real(wp), parameter :: cn_cutoff = 25.0_wp


   character(len=*), parameter :: real_format = "(es20.13)"
   character(len=25), parameter :: label_cutoff = "integral cutoff"
contains


   !> Determine the density matrix gradient based on the coupled-perturbed formalism
   subroutine get_density_matrix_gradient(mol,bas,wfn,list,dh0dr,dh0dL,doverlap,ddensity)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol      
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Wavefunction data
      type(wavefunction_type), intent(in) :: wfn
      !> Neighbour list
      type(adjacency_list), intent(in) :: list
      !> Derivative of the electronic energy w.r.t. coordinate displacements
      real(wp), intent(in) :: dh0dr(:, :, :)
      !> Derivative of the electronic energy w.r.t. the lattice vector
      real(wp), intent(in) :: dh0dL(:, :, :)
      !> Derivative of the electronic energy w.r.t. coordinate displacements
      real(wp), intent(in) :: doverlap(:, :, :)
      !> Derivative of the electronic energy w.r.t. coordinate displacements
      real(wp), intent(in) :: ddensity(:, :, :)
      

      ! Transform the derivatives of Fock and overlap matrix into the MO basis

      ! Determine the unitary transformation matrix from the unperturbed 
      ! to the perturbed basis coefficient matrix

      ! Determine the occupation number derivative based on the Fermi-Smearing

      ! Calculate the first derivative of the density matrix in the MO basis

      ! Transform the first derivative of the density matrix back to the AO basis

    end subroutine ceh_guess

end module tblite_ceh_coupled_perturbed

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

!> @dir tblite/exchange
!> Contains the exchange related interactions.

!> @file tblite/exchange.f90
!> Reexports of the exchange interaction containers.

!> Proxy module for handling exchange-type interactions
module tblite_exchange
   use mctc_env, only : wp
   use tblite_mulliken_kfock, only : mulliken_kfock_type
   use tblite_mulliken_kfock_matrix, only : mulliken_kfock_matrix, new_mulliken_exchange_matrix
   use tblite_mulliken_kfock_sqmbox, only : mulliken_kfock_sqmbox, new_mulliken_exchange_sqmbox
   use mctc_io, only : structure_type
   use tblite_exchange_type, only : exchange_type
   use tblite_basis_type, only : basis_type
   implicit none
   private

   public :: new_mulliken_exchange, exchange_type
   logical, parameter :: sqmbox = .false.

contains

!> Create a new Mulliken approximated exchange container
!> Screen the choice of the exchange implementation
subroutine new_mulliken_exchange(self, mol, bas, hubbard, shell_resolved_hubbard, &
   & average, gexp, frscale, omega, lrscale, incremental, allowsingle)
   !> Instance of the exchange container
   class(exchange_type), allocatable, intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Description of the basis set
   type(basis_type), intent(in) :: bas
   !> Hubbard parameter for all shells and species
   real(wp), intent(in) :: hubbard(:, :)
   !> Whether the Hubbard parameters are shell-dependent
   logical :: shell_resolved_hubbard
   !> Averaging scheme for Hubbard: 0 = arith., 1 = geom., 2 = arith.
   real(wp), intent(in) :: average
   !> Smoothening exponent; 1 = Mataga-Nishimoto, 2 = Klopman-Ohno
   real(wp), intent(in) :: gexp
   !> fullrange scale for K
   real(wp), intent(in) :: frscale
   !> omega if range seperated exchange is used
   real(wp), intent(in), optional :: omega
   !> long range sclae for range seperated exchange treatment
   real(wp), intent(in), optional :: lrscale
   !> allow single precision , incremental Fock Build
   logical, intent(in), optional :: allowsingle, incremental

   if (sqmbox) then
      block
         type(mulliken_kfock_sqmbox), allocatable :: tmp
         allocate(tmp)
         call new_mulliken_exchange_sqmbox(tmp, mol, bas, hubbard, &
            & allowsingle, incremental, average, gexp, frscale, omega, lrscale)
            call move_alloc(tmp, self)
      end block
   else
      block 
         type(mulliken_kfock_matrix), allocatable :: tmp
         allocate(tmp)
         call new_mulliken_exchange_matrix(tmp, mol, bas, hubbard, &
            & shell_resolved_hubbard, average, gexp, frscale, omega, lrscale)
         call move_alloc(tmp, self)
      end block
   end if

end subroutine new_mulliken_exchange

end module tblite_exchange

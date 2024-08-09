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

!> @dir tblite/exchange/mulliken_kfock
!> Contains the implementation of the Mulliken approximated Fock exchange

!> @file tblite/exchange/mulliken_kfock.f90
!> Provides an Mulliken approximated Fock exchange interaction

!> Proxy module for defining Mulliken approximated Fock exchange
module tblite_mulliken_kfock
   use tblite_mulliken_kfock_matrix, only : new_mulliken_exchange_matrix
   use tblite_mulliken_kfock_sqmbox, only : new_mulliken_exchange_sqmbox
   use tblite_mulliken_kfock_type, only : mulliken_kfock_type
   implicit none
   private

   public :: new_mulliken_exchange_matrix, new_mulliken_exchange_sqmbox
   public :: mulliken_kfock_type

end module tblite_mulliken_kfock

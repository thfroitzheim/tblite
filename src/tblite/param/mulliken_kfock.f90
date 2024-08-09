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

!> @file tblite/param/exchange.f90
!> Provides parameter record for the exchange

!> Defines model for the exchange parameters
module tblite_param_mulliken_kfock
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_toml, only : toml_table, get_value, set_value, add_table
   use tblite_param_serde, only : serde_record
   implicit none
   private

   public :: count

   character(len=*), parameter :: exchange_label = "mulliken_kfock", &
   & average_label = "average", exp_label = "gexp" , frscale_label = "frscale", &
   & omega_label = "omega", lrscale_label = "lrscale"

   !> Parametrization record for the exchange interaction electrostatics
   type, public, extends(serde_record) :: mulliken_kfock_record
      !> Full-range scaling for the K exchange
      real(wp) :: frscale 
      !> Range separation parameter
      real(wp) :: omega
      !> Scaling of the long range exchange
      real(wp) :: lrscale
      !> Averaging scheme for hardness with smooth interpolation
      !> (0) arithmetic, (1) geometric, (2) harmonic
      real(wp) :: average
      !> Smoothening exponent 
      !> (1) Mataga-Nishimoto, (2) Klopman-Ohno
      real(wp) :: gexp
   contains
      generic :: load => load_from_array
      generic :: dump => dump_to_array
      !> Read parametrization data from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization data to TOML data structure
      procedure :: dump_to_toml
      !> Read parametrization data from parameter array
      procedure, private :: load_from_array
      !> Write parametrization data to parameter array
      procedure, private :: dump_to_array
   end type

   !> Masking for the Mulliken exchange model
   type, public :: mulliken_kfock_mask
   end type mulliken_kfock_mask

   interface count
      module procedure :: count_mask
   end interface count

contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(mulliken_kfock_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   real(wp), allocatable :: last
   integer :: l, stat

   call get_value(table, exchange_label, child, requested=.false.)
   if (.not.associated(child)) then
      call fatal_error(error, "No entry for Mulliken exchange found")
      return
   end if

   call get_value(child, average_label, self%average, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for hardness averaging")
      return
   end if

   call get_value(child, exp_label, self%gexp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for effective interaction exponent")
      return
   end if

   call get_value(child, frscale_label, self%frscale, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid entry for fullrange scale of exchange")
      return
   end if

   if (child%has_key(omega_label)) then
      call get_value(child, omega_label, self%omega, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for range separation parameter")
         return
      end if

      call get_value(child, lrscale_label, self%lrscale, stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Invalid entry for long range scaling")
         return
      end if
   else
      self%omega = 0.0_wp
      self%lrscale = 0.0_wp
   endif

end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(mulliken_kfock_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   integer :: l

   call add_table(table, exchange_label, child)
   call set_value(child, average_label, self%average)
   call set_value(child, exp_label, self%gexp)
   call set_value(child, frscale_label, self%frscale)

   if (self%omega .ne. 0.0_wp) then
      call set_value(child, omega_label, self%omega)
      call set_value(child, lrscale_label, self%lrscale)
   end if

end subroutine dump_to_toml

 
!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(mulliken_kfock_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(mulliken_kfock_record), intent(in) :: base
   type(mulliken_kfock_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (mulliken_kfock_record)
      self = base
   end select

end subroutine load_from_array
 
!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(mulliken_kfock_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(mulliken_kfock_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

end subroutine dump_to_array
 
elemental function count_mask(mask) result(ncount)
   type(mulliken_kfock_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
end function count_mask

end module tblite_param_mulliken_kfock
 
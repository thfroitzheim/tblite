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
module tblite_param_exchange
    use mctc_env, only : wp, error_type, fatal_error
    use tblite_toml, only : toml_table, get_value, set_value, add_table
    use tblite_param_serde, only : serde_record
    implicit none
    private

   public :: count

    character(len=*), parameter :: exchange_label = "mulliken", average_label = "average_scheme" , &
    & exp_label = "expsmooth" , frscale_label = "frscale", omega_label = "omega", lrscale_label = "lrscale"

    type, public, extends(serde_record) :: exchange_record
        !> fullrange scale for the K scale
        real(wp) :: frscale 
        real(wp), allocatable :: omega, lrscale
        !> Averaging scheme for hardness:
        !> 0 = arith.
        !> 1 = geom.
        !> 2 = arith.
        integer :: average
        !> Smoothening exponent 
        !> 1 =  Mataga
        !> 2 = Klopman-type
        integer :: expsmooth
        logical :: mulliken

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

   type, public :: exchange_mask
   end type

   interface count
      module procedure :: count_mask
   end interface count

contains

   subroutine load_from_toml(self, table, error)
    !> Instance of the parametrization data
    class(exchange_record), intent(inout) :: self
    !> Data structure
    type(toml_table), intent(inout) :: table
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
 
    type(toml_table), pointer :: child
    real(wp), allocatable :: last
    integer :: l, stat
    real(wp) :: omega, lrscale
    call get_value(table, exchange_label, child, requested=.false.)

    if (.not.associated(child)) then
        call fatal_error(error, "No entry for mulliken exchange found, there is no other form of exchange integrated")
        return
    else
      self%mulliken = .true.
    end if
    
    call get_value(child, average_label, self%average, stat=stat)
    if (stat /= 0) then
        call fatal_error(error, "Invalid entry for averaging scheme of chemical hardness")
        return
    end if

    call get_value(child, exp_label, self%expsmooth, stat=stat)
    if (stat /= 0) then
        call fatal_error(error, "Invalid entry for smoothing scheme of chemical hardness")
        return
    end if

    call get_value(child, frscale_label, self%frscale, stat=stat)
    if (stat /= 0) then
        call fatal_error(error, "Invalid entry for fullrange scale of exchange")
        return
    end if
    
    if (child%has_key(omega_label)) then
      call get_value(child, omega_label, omega, stat=stat)
      allocate(self%omega, self%lrscale)
      self%omega = omega 
      call get_value(child, lrscale_label, lrscale, stat=stat)
      self%lrscale = lrscale
    endif


 end subroutine load_from_toml
 
 
 !> Write parametrization data to TOML datastructure
 subroutine dump_to_toml(self, table, error)
    !> Instance of the parametrization data
    class(exchange_record), intent(in) :: self
    !> Data structure
    type(toml_table), intent(inout) :: table
    !> Error handling
    type(error_type), allocatable, intent(out) :: error
 
    type(toml_table), pointer :: child
    integer :: l
 
    call add_table(table, exchange_label, child)
    call set_value(child, average_label, self%average)
    call set_value(child, exp_label, self%expsmooth)
    call set_value(child, frscale_label, self%frscale)

    if (allocated(self%omega)) then
        call set_value(child, omega_label, self%omega)
        call set_value(child, lrscale_label, self%lrscale)
     end if



 end subroutine dump_to_toml
 
 
 !> Read parametrization data from parameter array
 subroutine load_from_array(self, array, offset, base, mask, error)
    class(exchange_record), intent(inout) :: self
    real(wp), intent(in) :: array(:)
    integer, intent(inout) :: offset
    type(exchange_record), intent(in) :: base
    type(exchange_mask), intent(in) :: mask
    type(error_type), allocatable, intent(out) :: error
 
    select type(self)
    type is (exchange_record)
       self = base
    end select
 
 end subroutine load_from_array
 
 !> Write parametrization data to parameter array
 subroutine dump_to_array(self, array, offset, mask, error)
    class(exchange_record), intent(in) :: self
    real(wp), intent(inout) :: array(:)
    integer, intent(inout) :: offset
    type(exchange_mask), intent(in) :: mask
    type(error_type), allocatable, intent(out) :: error
 
 end subroutine dump_to_array
 
 elemental function count_mask(mask) result(ncount)
    type(exchange_mask), intent(in) :: mask
    integer :: ncount
    ncount = 0
 end function count_mask

end module
 
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

!> @file tblite/exchange/sqmbox.f90
!> Interface to the SQMBox for calculating the exchange matrix with GPU support

!> Mulliken approximated Fock exchange
module tblite_mulliken_kfock_sqmbox
   use mctc_env, only : wp, dp
   use mctc_io, only : structure_type
   use tblite_mulliken_kfock_type, only : mulliken_kfock_type
   use tblite_scf_potential, only : potential_type
   use tblite_container_cache, only : container_cache
   use tblite_exchange_cache, only : exchange_cache
   use tblite_basis_type, only : basis_type
   implicit none
   private

   public :: new_mulliken_exchange_sqmbox

   type, public, extends(mulliken_kfock_type) :: mulliken_kfock_sqmbox
      !> Hubbard parameter per shell with duplicates for reappearing atom types
      real(wp), allocatable :: hubbard(:)
      !> Arithmetic (0), geometric (1), and harmonic (2) average
      integer :: average
      !> Smoothening exponent (1 = Mataga-Nishimoto, 2 = Klopman-Ohno)
      integer :: gexp
      !> Allow single precision, 0 false, 1 true
      integer :: allowsingle
      !> Whether the fock matrix is build in incremental fashion, 0 false, 1 true
      integer :: incremental
   contains
      !> Evaluate exchange gamma matrix
      procedure :: get_exchange_matrix
      !> Evaluate uncontracted derivatives of exchange gamma matrix
      procedure :: get_exchange_derivs
      !> Calculate exchange contribution to the Fock matrix
      procedure :: get_KFock
      !> Calculate exchange contribution to the gradient
      procedure :: get_KGrad
   end type mulliken_kfock_sqmbox

   character(len=*), parameter :: label = "Mulliken approximated exchange - SQMBox"

   !> Interface to the SQMBox for calculating the exchange matrix with GPU support
   !> Should be handled as a subproject with an import
   interface  KFockSymSQM
      subroutine KFockSymSQM(allowsingle,incremental,nao,nat,nsh,aonum,sh2at,&
         temp,density,overlap,fock)
      import :: dp
      !> do single precision using integers, incremental Fock Build
      integer,intent(in) :: allowsingle, incremental
      !> Number of AOS
      integer,intent(in) :: nao
      !> Number of atoms
      integer,intent(in) :: nat
      !> Number of shells
      integer,intent(in) :: nsh
      !> Number of AOs in shell
      integer,intent(in) :: aonum(nsh)
      !> Convert from sh to at
      integer,intent(in) :: sh2at(nsh)
      !> density matrix
      real(kind=dp),intent(in) :: density(nao,nao)
      !> Overlap
      real(kind=dp),intent(in) :: overlap(nao,nao)
      !> Fock Matrix
      real(kind=dp),intent(inout) :: fock(nao,nao)
      real(kind=dp), intent(inout) :: temp(nao,nao)
      end subroutine
   end interface KFockSymSQM

   interface compute_gamma_rs
      subroutine compute_gamma_rs(nao,nat,nsh,aonum,sh2at,average,expsmooth,xyz,hubbard,&
         & frscale, omega, lrscale, temp)
      import :: dp
      !> Number of AOS
      integer,intent(in) :: nao
      !> Number of atoms
      integer,intent(in) :: nat
      !> Number of shells
      integer,intent(in) :: nsh
      !> Number of AOs in shell
      integer,intent(in) :: aonum(nsh)
      !> Convert from sh to at
      integer,intent(in) :: sh2at(nsh)
      !> Averaging scheme for Hubbard:
      !> 0 = arith.
      !> 1 = geom.
      !> 2 = arith.
      integer,intent(in) :: average
      !> Smoothening exponent
      !> 1 = Mataga
      !> 2 = Klopman-type
      integer,intent(in) :: expsmooth
      !> Cartesian coordinates
      real(kind=dp),intent(in) :: xyz(3,nat)
      !> Hubbard parameters 
      real(kind=dp),intent(in) :: hubbard(nsh)
      !>full-range K scale
      real(kind=dp),intent(in) :: frscale
      !> omega value for range seperated
      real(kind=dp):: omega
      !> longrange scale for range seperated functionals
      real(kind=dp) :: lrscale
      real(kind=dp), intent(out) :: temp(nao,nao)
      end subroutine

   end interface compute_gamma_rs

   interface compute_gamma_fr
      subroutine compute_gamma_fr(nao,nat,nsh,aonum,sh2at,average,expsmooth,xyz,hubbard,&
         & frscale, temp)
      import :: dp
      !> Number of AOS
      integer,intent(in) :: nao
      !> Number of atoms
      integer,intent(in) :: nat
      !> Number of shells
      integer,intent(in) :: nsh
      !> Number of AOs in shell
      integer,intent(in) :: aonum(nsh)
      !> Convert from sh to at
      integer,intent(in) :: sh2at(nsh)
      !> Averaging scheme for Hubbard:
      !> 0 = arith.
      !> 1 = geom.
      !> 2 = arith.
      integer,intent(in) :: average
      !> Smoothening exponent
      !> 1 =  Mataga
      !> 2 = Klopman-type
      integer,intent(in) :: expsmooth
      !> Cartesian coordinates
      real(kind=dp),intent(in) :: xyz(3,nat)
      !> Hubbard parameters
      real(kind=dp),intent(in) :: hubbard(nsh)
      !>full-range K scale
      real(kind=dp),intent(in) :: frscale

      real(kind=dp), intent(out) :: temp(nao,nao)
      end subroutine
   end interface compute_gamma_fr

contains


!> Create a new Mulliken approximated exchange container
subroutine new_mulliken_exchange_sqmbox(self, mol, bas, hubbard, &
   & allowsingle, incremental, average, gexp, frscale, omega, lrscale)
   !> Instance of the Mulliken exchange container
   type(mulliken_kfock_sqmbox), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Description of the basis set
   type(basis_type), intent(in) :: bas
   !> Hubbard parameter for all shells and species
   real(wp), intent(in) :: hubbard(:, :)
   !> allow single precision , incremental Fock Build
   logical, intent(in) :: allowsingle, incremental
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

   integer :: iish, ish, isp

   self%label = label

   allocate(self%nao_sh(bas%nsh), self%sh2at(bas%nsh))
   self%nao = bas%nao
   self%nsh = bas%nsh
   self%nao_sh = bas%nao_sh
   self%sh2at = bas%sh2at

   if (allowsingle .eqv. .true.) then
      self%allowsingle = 1
   else
      self%allowsingle= 0
   end if
   if (incremental .eqv. .true.) then
      self%incremental = 1
   else
      self%incremental= 0
   end if

   self%frscale = frscale
   if (present(omega).and.present(lrscale)) then
      self%omega = omega
      self%lrscale = lrscale
   else
      self%omega = 0.0_wp
      self%lrscale = 0.0_wp
   end if

   ! Convert the exponent to integer
   if(self%gexp == 1.0_wp .or. self%gexp == 2.0_wp) then
      self%gexp = int(gexp)
   else
      write(*,*) 'Error: Invalid exponent for smoothening'
   end if

   ! Convert the averaging scheme indicator to integer
   if(self%average == 1.0_wp .or. self%average == 2.0_wp) then
      self%average = int(average)
   else 
      write(*,*) 'Error: Invalid averaging scheme'
   end if

   allocate(self%hubbard(bas%nsh))
   iish = 0
   do isp = 1, mol%nid
      do ish = 1, bas%nsh_id(isp)
         iish = iish + 1
         self%hubbard(iish) = hubbard(ish, isp)
      end do
   end do 

end subroutine new_mulliken_exchange_sqmbox


!> Calculate exchange contribution to the Fock matrix
subroutine get_KFock(self, mol, cache, density, overlap, fock)
   !> Instance of the exchange container
   class(mulliken_kfock_sqmbox), intent(in) :: self
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

   integer :: spin
   real(wp), allocatable :: diff_P(:, :)

   if (self%incremental == 1) then
      if (.not.allocated(cache%ref_P).and..not.allocated(cache%ref_F)) then
         do spin = 1, size(density, 3)
            ! Do full Fock build in double precision
            call KFockSymSQM(0, 0, self%nao, mol%nat, self%nsh, self%nao_sh, self%sh2at, &
               & cache%gmat, density(:, :, spin), overlap, cache%prev_F(:, :, spin))
            fock(:, :, spin) = fock(:, :, spin) + 0.5_wp * cache%prev_F(:, :, spin)
         end do 

         ! Save the double precision reference density and Fock matrix
         ! Use the first variational iteration after the guess
         if (cache%has_prev_F) then 
            allocate(cache%ref_P(self%nao, self%nao, size(density, 3)), source=density)
            allocate(cache%ref_F(self%nao, self%nao, size(density, 3)), source=cache%prev_F)
         end if 
      else
         allocate(diff_P(self%nao, self%nao))
         do spin = 1, size(density, 3)
            ! Form the difference density w.r.t. the double precision reference
            diff_P(:, :)  = density(:, :, spin) - cache%ref_P(:, :, spin)
            ! Calculate the difference Fock matrix in single precision
            call KFockSymSQM(1, 0, self%nao, mol%nat, self%nsh, self%nao_sh, self%sh2at, &
               & cache%gmat, diff_P(:, :) , overlap, cache%prev_F(:, :, spin))
            ! Reconstruct the full Fock matrix
            cache%prev_F(:, :, spin) = cache%prev_F(:, :, spin) + cache%ref_F(:, :, spin) 
            fock(:, :, spin) = fock(:, :, spin) + 0.5_wp * cache%prev_F(:, :, spin)
         end do 
      end if
   else
      ! Do full Fock build in chosen precision
      do spin = 1, size(density, 3)
         call KFockSymSQM(0, 0, self%nao, mol%nat, self%nsh, self%nao_sh, self%sh2at, &
            & cache%gmat, density(:, :, spin), overlap, cache%prev_F(:, :, spin))
         fock(:, :, spin) = fock(:, :, spin) + 0.5_wp * cache%prev_F(:, :, spin)
      end do 
   endif

   !> Indicate that a previous Fock matrix is available
   cache%has_prev_F = .true.

end subroutine get_KFock


!> Calculate exchange contribution to the gradient
subroutine get_KGrad(self, mol, cache, density, overlap, sigma, gradient)
   !> Instance of the exchange container
   class(mulliken_kfock_sqmbox), intent(in) :: self
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

   gradient = 0.0_wp
   sigma = 0.0_wp
end subroutine get_KGrad


!> Evaluate exchange gamma matrix
subroutine get_exchange_matrix(self, mol, cache, gmat)
   !> Instance of the exchange container
   class(mulliken_kfock_sqmbox), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(exchange_cache), intent(inout) :: cache
   !> Coulomb matrix
   real(wp), contiguous, intent(out) :: gmat(:, :)
   
   integer :: igexp, nao, nsh 
   
   gmat(:, :) = 0.0_wp

   if (self%omega.ne.0.0_wp.and. self%lrscale.ne.0.0_wp) then
      call compute_gamma_rs(self%nao, mol%nat, self%nsh, self%nao_sh, self%sh2at,&
         & self%average, igexp, mol%xyz, self%hubbard, self%frscale, &
         & self%omega, self%frscale, gmat)
   else
      call compute_gamma_fr(self%nao, mol%nat, self%nsh, self%nao_sh, self%sh2at,&
         & self%average, igexp, mol%xyz, self%hubbard, self%frscale, gmat)
   end if

end subroutine get_exchange_matrix


!> Evaluate uncontracted derivatives of exchange gamma matrix
subroutine get_exchange_derivs(self, mol, cache, dgdr, dgdL)
   !> Instance of the exchange container
   class(mulliken_kfock_sqmbox), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(exchange_cache), intent(inout) :: cache
   !> Derivative of interactions with respect to cartesian displacements
   real(wp), contiguous, intent(out) :: dgdr(:, :, :)
   !> Derivative of interactions with respect to strain deformations
   real(wp), contiguous, intent(out) :: dgdL(:, :, :)
end subroutine get_exchange_derivs

end module tblite_mulliken_kfock_sqmbox

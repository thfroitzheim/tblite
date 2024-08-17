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

module test_gp3_xtb
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeq_guess
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gp3, only : new_gp3_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_ceh_singlepoint, only : ceh_singlepoint
   use tblite_ceh_ceh, only : new_ceh_calculator
   implicit none
   private

   public :: collect_gp3_xtb

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 0.0_wp * 3.166808578545117e-06_wp
   real(wp), parameter :: kt_ceh = 4000.0_wp * 3.166808578545117e-06_wp
   logical, parameter :: derive_basis = .true.

contains


!> Collect all exported unit tests
subroutine collect_gp3_xtb(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-lih", test_e_lih), &
      new_unittest("energy-mb01", test_e_mb01), &
      new_unittest("gradient-lih", test_g_lih), &
      new_unittest("gradient-mb02", test_g_mb02), &
      new_unittest("sigma-lih", test_s_lih), &
      new_unittest("sigma-mb03", test_s_mb03) &
      ]

end subroutine collect_gp3_xtb


subroutine test_e_gen(error, mol, ref)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reference GP3 energy
   real(wp), intent(in) :: ref

   type(context_type) :: ctx
   type(xtb_calculator) :: calc, calc_ceh
   type(wavefunction_type) :: wfn, wfn_ceh
   real(wp) :: energy
   real(wp), allocatable :: cn(:)
   real(wp), parameter :: accuracy = 1e-8_wp

   energy = 0.0_wp

   ! Setup GP3 calculator and wavefunction
   call new_gp3_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .false.)

   ! Setup CEH calculator and wavefunction
   call new_ceh_calculator(calc_ceh, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn_ceh, mol%nat, calc_ceh%bas%nsh, calc_ceh%bas%nao, 1,&
      & kt_ceh, .false.)
   ! Perform required CEH calculation
   call ceh_singlepoint(ctx, calc_ceh, mol, wfn_ceh, accuracy, .false., 0)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if
   ! Calculate the coordination number required for basis set scaling
   if (allocated(calc%bas%ncoord)) then
      allocate(cn(mol%nat))
      call calc%bas%ncoord%get_cn(mol, cn) 
   end if
   ! Scale the basis set with its charge and CN dependence
   call calc%bas%scale_basis(mol, wfn_ceh%qat(:, 1), cn, .true.)

   ! Use CEH as a guess
   wfn%qat(:, 1) = wfn_ceh%qat(:, 1)

   ! Perform GP3 calculation
   !call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, wfn_ceh=wfn_ceh)

   ! Check result
   call check(error, energy, ref, thr=1e-7_wp)

end subroutine test_e_gen

subroutine test_e_lih(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol
   ! Reference energy calculated with GP3 development code
   real(wp), parameter :: ref = -8.07788354_wp 

   call get_structure(mol, "MB16-43", "LiH")
   call test_e_gen(error, mol, ref)

end subroutine test_e_lih

subroutine test_e_mb01(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol
   ! Reference energy calculated with GP3 development code
   real(wp), parameter :: ref = -1278.37023512_wp 

   call get_structure(mol, "MB16-43", "01")
   call test_e_gen(error, mol, ref)

end subroutine test_e_mb01


subroutine test_g_gen(error, mol)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   type(context_type) :: ctx
   type(xtb_calculator) :: calc, calc_ceh
   type(wavefunction_type) :: wfn, wfn_ceh
   real(wp) :: energy, er, el 
   real(wp), allocatable :: gradient(:, :), numgrad(:, :), sigma(:, :)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: dqatnumdr(:, :, :)
   real(wp), allocatable :: ql(:), qr(:)
   real(wp), parameter :: accuracy = 1e-8_wp
   real(wp), parameter :: step = 1.0e-6_wp
   integer :: iat, ic

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat), sigma(3, 3))
   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   ! Setup GP3 calculator and wavefunction
   call new_gp3_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .false.)

   ! Setup CEH calculator and wavefunction
   allocate(ql(mol%nat), qr(mol%nat), dqatnumdr(3, mol%nat, mol%nat))
   call new_ceh_calculator(calc_ceh, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn_ceh, mol%nat, calc_ceh%bas%nsh, calc_ceh%bas%nao, 1,&
      & kt_ceh, .false.)

   if(.not.derive_basis) then
      ! Perform required CEH calculation 
      call ceh_singlepoint(ctx, calc_ceh, mol, wfn_ceh, accuracy, .false., 0)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
      ! Calculate the coordination number required for basis set scaling
      if (allocated(calc%bas%ncoord)) then
         allocate(cn(mol%nat))
         call calc%bas%ncoord%get_cn(mol, cn) 
      end if
      ! Scale the basis set only once with its charge and CN dependence
      call calc%bas%scale_basis(mol, wfn_ceh%qat(:, 1), cn, .true.)   
   end if 

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         
         ! Right hand
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         ! Obtain CEH charges 
         call ceh_singlepoint(ctx, calc_ceh, mol, wfn_ceh, accuracy, .false.)
         if (ctx%failed()) then
            call ctx%get_error(error)
            return
         end if
         qr = wfn_ceh%qat(:,1)
         if(derive_basis) then
            ! Calculate the coordination number required for basis set scaling
            call calc%bas%ncoord%get_cn(mol, cn) 
            ! Scale the basis set with its charge and CN dependence
            call calc%bas%scale_basis(mol, wfn_ceh%qat(:, 1), cn, .true.)
         end if
         ! Use CEH as a guess
         wfn%qat(:, 1) = wfn_ceh%qat(:, 1)
         ! Perform GP3 calculation
         !call xtb_singlepoint(ctx, mol, calc, wfn, acc, er, verbosity=0, wfn_ceh=wfn_ceh)
         
         ! Left hand
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         ! Obtain CEH charges 
         call ceh_singlepoint(ctx, calc_ceh, mol, wfn_ceh, accuracy, .false.)
         if (ctx%failed()) then
            call ctx%get_error(error)
            return
         end if
         ql = wfn_ceh%qat(:,1)
         if(derive_basis) then
            ! Calculate the coordination number required for basis set scaling
            call calc%bas%ncoord%get_cn(mol, cn) 
            ! Scale the basis set with its charge and CN dependence
            call calc%bas%scale_basis(mol, wfn_ceh%qat(:, 1), cn, .true.)
         end if 
         ! Use CEH as a guess
         wfn%qat(:, 1) = wfn_ceh%qat(:, 1)
         ! Perform GP3 calculation
         !call xtb_singlepoint(ctx, mol, calc, wfn, acc, el, verbosity=0, wfn_ceh=wfn_ceh)
         
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
         ! Do numerical CEH gradient calculation (intermediate solution)
         dqatnumdr(ic, iat, :) = 0.5_wp*(qr - ql)/step
      end do
   end do

   ! Setup CEH calculator and wavefunction with gradient
   call new_wavefunction(wfn_ceh, mol%nat, calc_ceh%bas%nsh, calc_ceh%bas%nao, 1,&
   & kt_ceh, .true.)
   call ceh_singlepoint(ctx, calc_ceh, mol, wfn_ceh, accuracy, .true., 0)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if
   wfn_ceh%dqatdr(:, :, :, 1) = dqatnumdr
   wfn_ceh%dqatdL(:, :, :, 1) = 0.0_wp 
 
   if(derive_basis) then
      ! Calculate the coordination number with gradient required for basis set scaling
      call calc%bas%ncoord%get_cn(mol, cn, dcndr, dcndL)
      ! Scale the basis set with its charge and CN dependence
      call calc%bas%scale_basis(mol, wfn_ceh%qat(:, 1), cn, .true., &
      & dqatdr=wfn_ceh%dqatdr(:, :, :, 1), dqatdL=wfn_ceh%dqatdL(:, :, :, 1), &
      & dcndr=dcndr, dcndL=dcndL)
   else
      ! Calculate the coordination number required for basis set scaling
      call calc%bas%ncoord%get_cn(mol, cn)
      ! Scale the basis set with its charge and CN dependence
      call calc%bas%scale_basis(mol, wfn_ceh%qat(:, 1), cn, .true.)
   end if

   ! GP3 calculation
   !call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, 0, wfn_ceh=wfn_ceh)

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, "Gradient of energy does not match")
      print'(3es21.14)', gradient
      print'("---")'
      print'(3es21.14)', numgrad
      print'("---")'
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_g_gen

subroutine test_g_lih(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_g_gen(error, mol)

end subroutine test_g_lih


subroutine test_g_mb02(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_g_gen(error, mol)

end subroutine test_g_mb02


subroutine test_s_gen(error, mol)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   type(context_type) :: ctx
   type(xtb_calculator) :: calc, calc_ceh
   type(wavefunction_type) :: wfn, wfn_ceh
   real(wp) :: energy, er, el
   real(wp), allocatable :: gradient(:, :), numsigma(:, :), sigma(:, :)
   real(wp), allocatable :: eps(:, :), xyz(:, :), lat(:, :)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: dqatnumdL(:, :, :)
   real(wp), allocatable :: ql(:), qr(:)
   real(wp), parameter :: unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: accuracy = 1e-8_wp
   real(wp), parameter :: step = 1.0e-6_wp
   integer :: ic, jc

   allocate(gradient(3, mol%nat), numsigma(3, 3), sigma(3, 3), &
      & eps(3, 3), xyz(3, mol%nat), lat(3, 3))
   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   ! Setup GP3 calculator and wavefunction
   call new_gp3_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .false.)

   ! Setup CEH calculator and wavefunction
   allocate(ql(mol%nat), qr(mol%nat), dqatnumdL(3, 3, mol%nat))
   call new_ceh_calculator(calc_ceh, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn_ceh, mol%nat, calc_ceh%bas%nsh, calc_ceh%bas%nao, 1,&
      & kt_ceh, .false.)

   if(.not.derive_basis) then
      ! Perform required CEH calculation 
      call ceh_singlepoint(ctx, calc_ceh, mol, wfn_ceh, accuracy, .false., 0)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
      ! Calculate the coordination number required for basis set scaling
      if (allocated(calc%bas%ncoord)) then
         allocate(cn(mol%nat))
         call calc%bas%ncoord%get_cn(mol, cn) 
      end if
      ! Scale the basis set only once with its charge and CN dependence
      call calc%bas%scale_basis(mol, wfn_ceh%qat(:, 1), cn, .true.)   
   end if 

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   if (any(mol%periodic)) &
   lat(:, :) = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         er = 0.0_wp
         el = 0.0_wp

         ! Right hand
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (any(mol%periodic)) &
         mol%lattice(:, :) = matmul(eps, lat)
         ! Obtain CEH charges 
         call ceh_singlepoint(ctx, calc_ceh, mol, wfn_ceh, accuracy, .false.)
         if (ctx%failed()) then
            call ctx%get_error(error)
            return
         end if
         qr = wfn_ceh%qat(:,1)
         if(derive_basis) then
            ! Calculate the coordination number required for basis set scaling
            call calc%bas%ncoord%get_cn(mol, cn) 
            ! Scale the basis set with its charge and CN dependence
            call calc%bas%scale_basis(mol, wfn_ceh%qat(:, 1), cn, .true.)
         end if
         ! Use CEH as a guess
         wfn%qat(:, 1) = wfn_ceh%qat(:, 1)
         ! Perform GP3 calculation
         !call xtb_singlepoint(ctx, mol, calc, wfn, acc, er, verbosity=0, wfn_ceh=wfn_ceh)
         
         ! Left hand
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (any(mol%periodic)) &
         mol%lattice(:, :) = matmul(eps, lat)
         ! Obtain CEH charges 
         call ceh_singlepoint(ctx, calc_ceh, mol, wfn_ceh, accuracy, .false.)
         if (ctx%failed()) then
            call ctx%get_error(error)
            return
         end if
         ql = wfn_ceh%qat(:,1)
         if(derive_basis) then
            ! Calculate the coordination number required for basis set scaling
            call calc%bas%ncoord%get_cn(mol, cn) 
            ! Scale the basis set with its charge and CN dependence
            call calc%bas%scale_basis(mol, wfn_ceh%qat(:, 1), cn, .true.)
         end if 
         ! Use CEH as a guess
         wfn%qat(:, 1) = wfn_ceh%qat(:, 1)
         ! Perform GP3 calculation
         !call xtb_singlepoint(ctx, mol, calc, wfn, acc, el, verbosity=0, wfn_ceh=wfn_ceh)
         
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (any(mol%periodic)) &
         mol%lattice(:, :) = lat
         numsigma(jc, ic)  = 0.5_wp*(er - el)/step
         ! Do numerical CEH gradient calculation (intermediate solution)
         dqatnumdL(ic, jc, :) = 0.5_wp*(qr - ql)/step
      end do
   end do

   ! Setup CEH calculator and wavefunction with gradient
   call new_wavefunction(wfn_ceh, mol%nat, calc_ceh%bas%nsh, calc_ceh%bas%nao, 1,&
   & kt_ceh, .true.)
   call ceh_singlepoint(ctx, calc_ceh, mol, wfn_ceh, accuracy, .true., 0)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if
   wfn_ceh%dqatdr(:, :, :, 1) = 0.0_wp
   wfn_ceh%dqatdL(:, :, :, 1) = dqatnumdL
 
   if(derive_basis) then
      ! Calculate the coordination number with gradient required for basis set scaling
      call calc%bas%ncoord%get_cn(mol, cn, dcndr, dcndL)
      ! Scale the basis set with its charge and CN dependence
      call calc%bas%scale_basis(mol, wfn_ceh%qat(:, 1), cn, .true., &
      & dqatdr=wfn_ceh%dqatdr(:, :, :, 1), dqatdL=wfn_ceh%dqatdL(:, :, :, 1), &
      & dcndr=dcndr, dcndL=dcndL)
   else
      ! Calculate the coordination number required for basis set scaling
      call calc%bas%ncoord%get_cn(mol, cn)
      ! Scale the basis set with its charge and CN dependence
      call calc%bas%scale_basis(mol, wfn_ceh%qat(:, 1), cn, .true.)
   end if

   ! GP3 calculation
   !call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, 0, wfn_ceh=wfn_ceh)

   if (any(abs(sigma - numsigma) > thr2)) then
      call test_failed(error, "Sigma of energy does not match")
      print'(3es21.14)', sigma
      print'("---")'
      print'(3es21.14)', numsigma
      print'("---")'
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_s_gen

subroutine test_s_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_s_gen(error, mol)

end subroutine test_s_lih

subroutine test_s_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_s_gen(error, mol)

end subroutine test_s_mb03

end module test_gp3_xtb
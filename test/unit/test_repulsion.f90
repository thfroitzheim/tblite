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

module test_repulsion
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_container_cache, only : container_cache
   use tblite_repulsion, only : repulsion_type, new_repulsion
   use tblite_context_type, only : context_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_ceh_singlepoint, only : ceh_singlepoint
   use tblite_ceh_ceh, only : new_ceh_calculator
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   implicit none
   private

   public :: collect_repulsion

   ! CEH temperature
   real(wp), parameter :: kt = 4000.0_wp * 3.166808578545117e-06_wp
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine repulsion_maker(rep, mol)
         import :: repulsion_type, structure_type
         class(repulsion_type), allocatable, intent(out) :: rep
         type(structure_type), intent(in) :: mol
      end subroutine repulsion_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_repulsion(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-gfn1", test_e_effective_m01), &
      new_unittest("energy-gfn2", test_e_effective_m02), &
      new_unittest("energy-gp3", test_e_effective_m03), &
      new_unittest("energy-pbc-gfn2", test_e_effective_uracil), &
      new_unittest("gradient-gfn1", test_g_effective_m03), &
      new_unittest("gradient-gfn2", test_g_effective_m04), &
      new_unittest("gradient-gp3", test_g_effective_m05), &
      new_unittest("gradient-pbc-gfn2", test_g_effective_urea), &
      new_unittest("sigma-gfn1", test_s_effective_m05), &
      new_unittest("sigma-gfn2", test_s_effective_m06), &
      new_unittest("sigma-gp3", test_s_effective_m07), &
      new_unittest("sigma-pbc-gfn2", test_s_effective_succinic) &
      ]

end subroutine collect_repulsion


!> Factory to create repulsion objects based on GFN1-xTB
subroutine make_repulsion_gfn1(rep, mol)
   class(repulsion_type), allocatable, intent(out) :: rep
   type(structure_type), intent(in) :: mol

   real(wp), parameter :: alpha_gfn1(1:20) = [&
      & 2.209700_wp, 1.382907_wp, 0.671797_wp, 0.865377_wp, 1.093544_wp, &
      & 1.281954_wp, 1.727773_wp, 2.004253_wp, 2.507078_wp, 3.038727_wp, &
      & 0.704472_wp, 0.862629_wp, 0.929219_wp, 0.948165_wp, 1.067197_wp, &
      & 1.200803_wp, 1.404155_wp, 1.323756_wp, 0.581529_wp, 0.665588_wp]
   real(wp), parameter :: zeff_gfn1(1:20) = [&
      &  1.116244_wp,  0.440231_wp,  2.747587_wp,  4.076830_wp,  4.458376_wp, &
      &  4.428763_wp,  5.498808_wp,  5.171786_wp,  6.931741_wp,  9.102523_wp, &
      & 10.591259_wp, 15.238107_wp, 16.283595_wp, 16.898359_wp, 15.249559_wp, &
      & 15.100323_wp, 17.000000_wp, 17.153132_wp, 20.831436_wp, 19.840212_wp]
   real(wp), allocatable :: alpha(:), zeff(:), kexp

   alpha = alpha_gfn1(mol%num)
   zeff = zeff_gfn1(mol%num)
   kexp = 1.5_wp
   call new_repulsion(rep, mol, rep_type="gfn", zeff=zeff, alpha=alpha, rexp=1.0_wp, &
      & kexp=1.5_wp, kexp_light=kexp)
end subroutine make_repulsion_gfn1

!> Factory to create repulsion objects based on GFN2-xTB
subroutine make_repulsion_gfn2(rep, mol)
   class(repulsion_type), allocatable, intent(out) :: rep
   type(structure_type), intent(in) :: mol

   real(wp), parameter :: alpha_gfn2(1:20) = [&
      & 2.213717_wp, 3.604670_wp, 0.475307_wp, 0.939696_wp, 1.373856_wp, &
      & 1.247655_wp, 1.682689_wp, 2.165712_wp, 2.421394_wp, 3.318479_wp, &
      & 0.572728_wp, 0.917975_wp, 0.876623_wp, 1.187323_wp, 1.143343_wp, &
      & 1.214553_wp, 1.577144_wp, 0.896198_wp, 0.482206_wp, 0.683051_wp]
   real(wp), parameter :: zeff_gfn2(1:20) = [&
      &  1.105388_wp,  1.094283_wp,  1.289367_wp,  4.221216_wp,  7.192431_wp, &
      &  4.231078_wp,  5.242592_wp,  5.784415_wp,  7.021486_wp, 11.041068_wp, &
      &  5.244917_wp, 18.083164_wp, 17.867328_wp, 40.001111_wp, 19.683502_wp, &
      & 14.995090_wp, 17.353134_wp,  7.266606_wp, 10.439482_wp, 14.786701_wp]
   real(wp), allocatable :: alpha(:), zeff(:), kexp

   alpha = alpha_gfn2(mol%num)
   zeff = zeff_gfn2(mol%num)
   kexp = 1.0_wp
   call new_repulsion(rep, mol, rep_type="gfn", zeff=zeff, alpha=alpha, rexp=1.0_wp, &
      & kexp=1.5_wp, kexp_light=kexp)
end subroutine make_repulsion_gfn2

!> Factory to create repulsion objects based on GP3-xTB
subroutine make_repulsion_gp3(rep, mol)
   class(repulsion_type), allocatable, intent(out) :: rep
   type(structure_type), intent(in) :: mol

   real(wp), parameter :: alpha_gp3(1:20) = [&
      &  0.2342286016_wp,  0.8803317221_wp,  0.5146531158_wp,  0.4280997510_wp, &
      &  0.5058996423_wp,  0.4868704255_wp,  0.4608909036_wp,  0.3755559418_wp, &
      &  0.4975431361_wp,  0.5733862867_wp,  0.4461313165_wp,  0.4822649381_wp, &
      &  0.4096852133_wp,  0.4439802705_wp,  0.4668800640_wp,  0.4902942052_wp, &
      &  0.4513441213_wp,  0.5404520531_wp,  0.5123109723_wp,  0.4594484167_wp]
   real(wp), parameter :: zeff_gp3(1:20) = [&
      &  0.8349280246_wp,  0.5164194548_wp,  3.3186588270_wp,  2.7483927859_wp, &
      &  3.1514033253_wp,  3.2715032000_wp,  3.0825487458_wp,  2.9582589564_wp, &
      &  2.5414641410_wp,  2.1050638756_wp,  6.2584193721_wp,  5.6489204656_wp, &
      &  4.5532494572_wp,  5.5733558960_wp,  6.3218383567_wp,  5.5237011332_wp, &
      &  4.8839050775_wp,  5.8546075966_wp,  8.9371689515_wp,  7.0245743182_wp]
   real(wp), parameter :: kcn_gp3(1:20) = [&
      &  0.1888951308_wp, -0.0150613309_wp, -0.1711013242_wp, -0.0049403131_wp, &
      &  0.0360342527_wp,  0.0233077978_wp,  0.0218716480_wp,  0.0859251997_wp, &
      & -0.0164621829_wp,  0.0108864040_wp, -0.0356858269_wp, -0.1139920806_wp, &
      & -0.0159103005_wp,  0.0336457887_wp,  0.0227112823_wp,  0.0673624372_wp, &
      &  0.0412699810_wp,  0.0130095413_wp, -0.1619818551_wp, -0.2030414817_wp]
   real(wp), parameter :: kq_gp3(1:20) = [&
      &  0.0008271593_wp, -0.3522609502_wp,  0.1353380841_wp,  0.0450706434_wp, &
      &  0.0525063232_wp,  0.0007776371_wp,  0.0005586691_wp,  0.0007880504_wp, &
      &  0.0101403049_wp, -0.0852263855_wp,  0.3882216335_wp,  0.2065255534_wp, &
      &  0.1501898531_wp,  0.0116197630_wp,  0.0249841827_wp,  0.0076204783_wp, &
      &  0.0075788259_wp,  0.1054809329_wp,  0.4208321354_wp,  0.1771453922_wp]
   real(wp), parameter :: rcov_gp3(1:20) = [&
      &  0.6719079170_wp,  1.7019143884_wp,  1.2509057985_wp,  1.1064591188_wp, &
      &  1.0462436357_wp,  0.9852758482_wp,  0.8489038165_wp,  0.6513920820_wp, &
      &  0.7811154636_wp,  1.1152976802_wp,  1.1127693304_wp,  1.3352730420_wp, &
      &  1.3915878165_wp,  1.2195716461_wp,  1.1057052621_wp,  1.0988287529_wp, &
      &  0.9886349126_wp,  1.1352768544_wp,  1.6450241053_wp,  1.4869197542_wp]
   real(wp), parameter :: rcov_cn_gp3(1:20) = [&
      &  0.2806743503_wp, -0.2367094066_wp,  1.3199466509_wp,  1.0928398320_wp, &
      &  1.2163297359_wp,  1.2892723926_wp,  1.1093215479_wp,  0.9522993236_wp, &
      &  0.7601560354_wp,  0.9091963275_wp,  1.7896424369_wp,  1.5939753643_wp, &
      &  1.1783917691_wp,  1.1091054929_wp,  1.6808311996_wp,  1.6602686396_wp, &
      &  1.5106000154_wp,  1.3188860170_wp,  2.5074795871_wp,  1.7241829531_wp]
   
   real(wp), allocatable :: alpha(:), zeff(:), rcov(:), kcn(:), kq(:), rcov_cn(:)
   real(wp) :: exp_cn 

   alpha = alpha_gp3(mol%num)
   zeff = zeff_gp3(mol%num)
   kcn = kcn_gp3(mol%num)
   kq = kq_gp3(mol%num)
   rcov = rcov_gp3(mol%num)
   rcov_cn = rcov_cn_gp3(mol%num)
   exp_cn = 1.812_wp
   call new_repulsion(rep, mol, rep_type="gp3", zeff=zeff, alpha=alpha, rexp=2.0_wp, &
      & kcn=kcn, kq=kq, rcov_rep=rcov, rcov_cn=rcov_cn, exp_cn=exp_cn)
end subroutine make_repulsion_gp3

subroutine test_generic(error, mol, make_repulsion, ref, chargedep)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Factory to produce repulsion objects
   procedure(repulsion_maker) :: make_repulsion

   !> Reference value to compare against
   real(wp), intent(in) :: ref

   !> Flag to signal charge dependence of the repulsion
   logical, intent(in) :: chargedep

   class(repulsion_type), allocatable :: rep
   type(container_cache) :: cache
   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type), allocatable :: wfn
   real(wp) :: energy(mol%nat)
   real(wp), parameter :: accuracy = 1e-8_wp
   
   if(chargedep) then
      ! Obtain CEH charges for charge dependence
      call new_ceh_calculator(calc, mol, error)
      if (allocated(error)) return
      allocate(wfn)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .false.)
      ctx%verbosity = 0
      call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .false.)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
   end if

   energy = 0.0_wp
   call make_repulsion(rep, mol)
   call rep%update(mol, cache, wfn)
   call rep%get_engrad(mol, cache, energy)

   call check(error, sum(energy), ref, thr=thr2)
   if (allocated(error)) then
      print*,sum(energy), ref, sum(energy) - ref, thr, thr2
   end if

end subroutine test_generic


subroutine test_numgrad(error, mol, make_repulsion, chargedep)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Factory to create a repulsion object
   procedure(repulsion_maker) :: make_repulsion

   !> Flag to signal charge dependence of the repulsion
   logical, intent(in) :: chargedep

   integer :: iat, ic
   class(repulsion_type), allocatable :: rep
   type(container_cache) :: cache
   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type), allocatable :: wfn
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :), dqatnumdr(:, :, :)
   real(wp), allocatable :: ql(:), qr(:)
   real(wp), parameter :: accuracy = 1e-8_wp
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call make_repulsion(rep, mol)

   if(chargedep) then
      ! Setup CEH for charge dependence
      allocate(ql(mol%nat), qr(mol%nat), dqatnumdr(3, mol%nat, mol%nat))
      call new_ceh_calculator(calc, mol, error)
      if (allocated(error)) return
      allocate(wfn)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .false.)
      ctx%verbosity = 0
   end if

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp

         ! Right hand
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         if(chargedep) then
            ! Obtain CEH charges for charge dependence
            call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .false.)
            if (ctx%failed()) then
               call ctx%get_error(error)
               return
            end if
            qr = wfn%qat(:,1)
         end if
         ! Update and calculate repulsion
         call rep%update(mol, cache, wfn)
         call rep%get_engrad(mol, cache, er)
         
         ! Left hand
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         if(chargedep) then
            ! Obtain CEH charges for charge dependence
            call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .false.)
            if (ctx%failed()) then
               call ctx%get_error(error)
               return
            end if
            ql = wfn%qat(:,1)
         end if
         ! Update and calculate repulsion
         call rep%update(mol, cache, wfn)
         call rep%get_engrad(mol, cache, el)
      
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
         if(chargedep) then
            ! CEH gradient (intermediate solution)
            dqatnumdr(ic, iat, :) = 0.5_wp*(qr - ql)/step
         end if
      end do
   end do

   if(chargedep) then
      ! Obtain CEH charges and derivatives for charge adaptation
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .true.)
      call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .true.)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
      wfn%dqatdr(:, :, :, 1) = dqatnumdr
      wfn%dqatdL(:, :, :, 1) = 0.0_wp 
   end if

   call rep%update(mol, cache, wfn)
   call rep%get_engrad(mol, cache, energy, gradient, sigma)

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, "Gradient of repulsion energy does not match")
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol, make_repulsion, chargedep)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Factory to create a repulsion object
   procedure(repulsion_maker) :: make_repulsion

   !> Flag to signal charge dependence of the repulsion
   logical, intent(in) :: chargedep

   integer :: ic, jc
   class(repulsion_type), allocatable :: rep
   type(container_cache) :: cache
   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type), allocatable :: wfn
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3), eps(3, 3), numsigma(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :), lat(:, :), dqatnumdL(:, :, :)
   real(wp), allocatable :: ql(:), qr(:)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: accuracy = 1e-8_wp
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), xyz(3, mol%nat), lat(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call make_repulsion(rep, mol)

   
   if(chargedep) then
      ! Setup CEH for charge dependence
      allocate(ql(mol%nat), qr(mol%nat), dqatnumdL(3, 3, mol%nat))
      call new_ceh_calculator(calc, mol, error)
      if (allocated(error)) return
      allocate(wfn)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .false.)
      ctx%verbosity = 0
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
         if(chargedep) then
            ! Obtain CEH charges for charge dependence
            call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .false.)
            if (ctx%failed()) then
               call ctx%get_error(error)
               return
            end if
            qr = wfn%qat(:,1)
         end if 
         ! Update and calculate repulsion
         call rep%update(mol, cache, wfn)
         call rep%get_engrad(mol, cache, er)

         ! Left hand
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (any(mol%periodic)) &
         mol%lattice(:, :) = matmul(eps, lat)
         if(chargedep) then
            ! Obtain CEH charges for charge dependence
            call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .false.)
            if (ctx%failed()) then
               call ctx%get_error(error)
               return
            end if
            ql = wfn%qat(:,1)
         end if
         ! Update and calculate repulsion
         call rep%update(mol, cache, wfn)
         call rep%get_engrad(mol, cache, el)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (any(mol%periodic)) &
         mol%lattice(:, :) = lat
         numsigma(jc, ic) = 0.5_wp*(sum(er) - sum(el))/step
         if(chargedep) then
            ! CEH gradient (intermediate solution)
            dqatnumdL(ic, jc, :) = 0.5_wp*(qr - ql)/step
         end if
      end do
   end do

   if(chargedep) then
      ! Obtain CEH charges and derivatives for charge adaptation
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .true.)
      call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .true.)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
      wfn%dqatdr(:, :, :, 1) = 0.0_wp
      wfn%dqatdL(:, :, :, 1) = dqatnumdL 
   end if

   call rep%update(mol, cache, wfn)
   call rep%get_engrad(mol, cache, energy, gradient, sigma)

   if (any(abs(sigma - numsigma) > thr2)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma


subroutine test_e_effective_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, make_repulsion_gfn1, 0.16777923624986593_wp, .false.)

end subroutine test_e_effective_m01


subroutine test_e_effective_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, make_repulsion_gfn2, 0.10745931926703985_wp, .false.)

end subroutine test_e_effective_m02

subroutine test_e_effective_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_generic(error, mol, make_repulsion_gp3, 1.69091292_wp, .true.)

end subroutine test_e_effective_m03

subroutine test_e_effective_uracil(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "uracil")
   call test_generic(error, mol, make_repulsion_gfn2, 1.0401472262740301_wp, .false.)

end subroutine test_e_effective_uracil


subroutine test_g_effective_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, make_repulsion_gfn1, .false.)

end subroutine test_g_effective_m03


subroutine test_g_effective_m04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol, make_repulsion_gfn2, .false.)

end subroutine test_g_effective_m04


subroutine test_g_effective_m05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numgrad(error, mol, make_repulsion_gp3, .true.)

end subroutine test_g_effective_m05


subroutine test_g_effective_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "urea")
   call test_numgrad(error, mol, make_repulsion_gfn2, .false.)

end subroutine test_g_effective_urea


subroutine test_s_effective_m05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numsigma(error, mol, make_repulsion_gfn1, .false.)

end subroutine test_s_effective_m05


subroutine test_s_effective_m06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol, make_repulsion_gfn2, .false.)

end subroutine test_s_effective_m06


subroutine test_s_effective_m07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "07")
   call test_numsigma(error, mol, make_repulsion_gp3, .true.)

end subroutine test_s_effective_m07


subroutine test_s_effective_succinic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "succinic")
   call test_numsigma(error, mol, make_repulsion_gfn2, .false.)

end subroutine test_s_effective_succinic


end module test_repulsion

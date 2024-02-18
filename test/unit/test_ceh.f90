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

module test_ceh
   use iso_fortran_env, only: output_unit
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
   & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_lapack_solver, only : lapack_solver, lapack_algorithm
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_cutoff, only : get_lattice_points
   use tblite_basis_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_xtb_h0, only : tb_hamiltonian, new_hamiltonian
   use tblite_data_covrad, only : get_covalent_rad
   use tblite_data_paulingen, only : get_pauling_en
   use tblite_ncoord_erf
   use tblite_ncoord_erf_en
   use tblite_ncoord_type, only : get_coordination_number
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_ceh_singlepoint, only : ceh_guess
   use tblite_ceh_ceh, only : ceh_h0spec, new_ceh_calculator
   use tblite_ceh_h0, only : get_scaled_selfenergy, get_hamiltonian, get_hamiltonian_gradient
   use tblite_scf, only: new_potential, potential_type
   use tblite_blas, only: gemv
   use tblite_container, only : container_type, container_cache
   use tblite_external_field, only : electric_field
   implicit none
   private

   public :: collect_ceh

   real(wp), parameter :: kt = 5000.0_wp * 3.166808578545117e-06_wp
   real(wp), parameter :: thr2 = 1.0e2_wp*sqrt(epsilon(1.0_wp))

contains

   !> Collect all exported unit tests
   subroutine collect_ceh(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("scaled-selfenergy-H2", test_scaled_selfenergy_h2), &
         new_unittest("scaled-selfenergy-LiH", test_scaled_selfenergy_lih), &
         new_unittest("scaled-selfenergy-S2", test_scaled_selfenergy_s2), &
         new_unittest("scaled-selfenergy-SiH4", test_scaled_selfenergy_sih4), &
         new_unittest("scaled-selfenergy_grad-H2", test_scaled_selfenergy_numgrad_h2), &
         new_unittest("scaled-selfenergy_grad-LiH", test_scaled_selfenergy_numgrad_lih), &
         new_unittest("scaled-selfenergy_grad-S2", test_scaled_selfenergy_numgrad_s2), &
         new_unittest("scaled-selfenergy_grad-SiH4", test_scaled_selfenergy_numgrad_sih4), &
         new_unittest("hamiltonian-H2", test_hamiltonian_h2), &
         new_unittest("hamiltonian-LiH", test_hamiltonian_lih), &
         new_unittest("hamiltonian-S2", test_hamiltonian_s2), &
         new_unittest("hamiltonian-SiH4", test_hamiltonian_sih4), &
         new_unittest("hamiltonian_grad-H2", test_hamiltonian_numgrad_h2), &
         new_unittest("hamiltonian_grad-LiH", test_hamiltonian_numgrad_lih), &
         new_unittest("hamiltonian_grad-S2", test_hamiltonian_numgrad_s2), &
         new_unittest("hamiltonian_grad-PCl", test_hamiltonian_numgrad_pcl), &
         new_unittest("hamiltonian_grad-SiH4", test_hamiltonian_numgrad_sih4), &
         new_unittest("hamiltonian_grad-CeCl3", test_hamiltonian_numgrad_cecl3), &
         new_unittest("overlap_diat-H2", test_overlap_diat_h2), &
         new_unittest("overlap_diat-LiH", test_overlap_diat_lih), &
         new_unittest("overlap_diat-S2", test_overlap_diat_s2), &
         new_unittest("overlap_diat-SiH4", test_overlap_diat_sih4), &
         new_unittest("q-mol-h2", test_q_h2), &
         new_unittest("q-mol-lih", test_q_lih), &
         new_unittest("q-mol-1", test_q_mb01), &
         new_unittest("q-mol-2", test_q_mb02), &
         new_unittest("q-mol-3", test_q_mb03), &
         new_unittest("q-mol-4", test_q_mb04), &
         new_unittest("q-chrgd-efield-mol", test_q_ef_chrg_mb01), &
         new_unittest("d-mol", test_d_mb01), &
         new_unittest("d-field-mol", test_d_field_mb04) &
         ]

   end subroutine collect_ceh

   !> Testing on the CEH basis
   subroutine make_basis(bas, mol, ng)
      type(basis_type), intent(out) :: bas
      type(structure_type), intent(in) :: mol
      integer, intent(in) :: ng

      integer, parameter :: nsh(86) = [&
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, & ! 1-20
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & ! 21-40
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, & ! 41-60
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, & ! 61-80
      & 3, 3, 3, 3, 3, 3]
      integer, parameter :: lsh(3, 86) = reshape([&
      & 0, 0, 0,  0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, & ! 1-7
      & 0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 8-14
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2, & ! 15-21
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 22-28
      & 0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 29-35
      & 0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 36-42
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2, & ! 43-49
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2, & ! 50-56
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 57-63
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 64-70
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 71-77
      & 0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, & ! 78-84
      & 0, 1, 2,  0, 1, 2], &
      & shape(lsh))
      integer, parameter :: pqn(3, 86) = reshape([&
      & 1, 0, 0,  1, 0, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0, & ! 1-7
      & 2, 2, 0,  2, 2, 0,  2, 2, 0,  3, 3, 0,  3, 3, 3,  3, 3, 3,  3, 3, 3, & ! 8-14
      & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  4, 4, 0,  4, 4, 3,  4, 4, 3, & ! 15-21
      & 4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3,  4, 4, 3, & ! 22-28
      & 4, 4, 3,  4, 4, 0,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4,  4, 4, 4, & ! 29-35
      & 4, 4, 4,  5, 5, 0,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4, & ! 36-42
      & 5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 4,  5, 5, 0,  5, 5, 5, & ! 43-49
      & 5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  5, 5, 5,  6, 6, 0,  6, 6, 5, & ! 50-56
      & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, & ! 57-63
      & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, & ! 64-70
      & 6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, & ! 71-77
      & 6, 6, 5,  6, 6, 5,  6, 6, 0,  6, 6, 5,  6, 6, 5,  6, 6, 5,  6, 6, 5, & ! 78-84
      & 6, 6, 5,  6, 6, 5], &
      & shape(pqn))
      real(wp), parameter :: zeta(3, 86) = reshape([&
      & 1.23363166_wp, 0.00000000_wp, 0.00000000_wp, 2.27004605_wp, 0.00000000_wp, 0.00000000_wp, &
      & 0.86185456_wp, 1.42017184_wp, 0.00000000_wp, 1.76817995_wp, 1.44095844_wp, 0.00000000_wp, &
      & 2.06339837_wp, 1.52051807_wp, 0.00000000_wp, 2.56058582_wp, 1.86484737_wp, 0.00000000_wp, &
      & 2.71233631_wp, 2.19848968_wp, 0.00000000_wp, 3.21585650_wp, 2.41309737_wp, 0.00000000_wp, &
      & 3.82146807_wp, 2.63063636_wp, 0.00000000_wp, 4.62721228_wp, 2.53599954_wp, 0.00000000_wp, &
      & 0.93221172_wp, 1.55333839_wp, 0.00000000_wp, 1.77220557_wp, 1.59942632_wp, 2.98596647_wp, &
      & 2.26040231_wp, 1.78718151_wp, 2.00990188_wp, 1.85259089_wp, 1.81733349_wp, 1.65269988_wp, &
      & 2.65701241_wp, 2.03189759_wp, 2.03883661_wp, 2.60609998_wp, 2.16530440_wp, 2.41888232_wp, &
      & 2.78818934_wp, 2.24732894_wp, 1.99081182_wp, 2.55424399_wp, 2.20946190_wp, 1.93619550_wp, &
      & 1.73713827_wp, 1.33788617_wp, 0.00000000_wp, 2.47982574_wp, 1.07250770_wp, 2.11920764_wp, &
      & 2.22449249_wp, 1.55418319_wp, 2.00953578_wp, 2.58879616_wp, 0.99441077_wp, 1.88561781_wp, &
      & 3.04370654_wp, 4.03007600_wp, 1.66329169_wp, 2.25012727_wp, 2.70681556_wp, 1.67501904_wp, &
      & 2.20605319_wp, 2.82019792_wp, 1.86102254_wp, 1.57297015_wp, 1.98621494_wp, 2.83790684_wp, &
      & 1.80826602_wp, 1.73675835_wp, 2.79767448_wp, 2.00758945_wp, 2.25075692_wp, 2.98291663_wp, &
      & 2.18159986_wp, 2.38459096_wp, 3.09502522_wp, 2.26376756_wp, 2.20362977_wp, 0.00000000_wp, &
      & 2.63822153_wp, 2.06752328_wp, 2.11361643_wp, 2.52891955_wp, 2.19441794_wp, 1.77661998_wp, &
      & 3.55667605_wp, 2.42075463_wp, 1.46579772_wp, 2.89652631_wp, 2.45421858_wp, 2.27883625_wp, &
      & 3.28921099_wp, 2.56526915_wp, 1.64501640_wp, 5.20988189_wp, 2.84336725_wp, 2.75838814_wp, &
      & 1.26972917_wp, 1.88730596_wp, 0.00000000_wp, 1.86880714_wp, 1.78546342_wp, 2.16012236_wp, &
      & 0.92001877_wp, 1.45732462_wp, 2.22901395_wp, 6.50647305_wp, 1.43202338_wp, 2.11971490_wp, &
      & 2.10973371_wp, 2.79944781_wp, 2.01897369_wp, 2.58413333_wp, 3.02795359_wp, 2.08733665_wp, &
      & 2.62141555_wp, 3.13487625_wp, 2.13259872_wp, 2.73984499_wp, 2.18167834_wp, 2.54609647_wp, &
      & 1.84057176_wp, 2.97482636_wp, 3.10693700_wp, 1.75622839_wp, 3.39424756_wp, 3.20265306_wp, &
      & 3.05018811_wp, 2.34951987_wp, 3.35332952_wp, 2.41999128_wp, 2.28892954_wp, 0.00000000_wp, &
      & 2.87813961_wp, 2.44659724_wp, 2.75773502_wp, 3.03823214_wp, 2.32082155_wp, 1.77513328_wp, &
      & 2.68750711_wp, 2.38565373_wp, 2.12596190_wp, 2.81071790_wp, 2.45274786_wp, 2.01871821_wp, &
      & 2.90686956_wp, 2.49377102_wp, 1.90073732_wp, 4.17531340_wp, 2.86937955_wp, 2.96894812_wp, &
      & 1.24299361_wp, 1.99142040_wp, 0.00000000_wp, 1.31400366_wp, 1.16438481_wp, 2.12759606_wp, &
      & 2.81737350_wp, 1.69863323_wp, 2.27369715_wp, 2.84503901_wp, 1.46018192_wp, 2.53498936_wp, &
      & 2.81697107_wp, 1.47545307_wp, 2.54350275_wp, 2.78890313_wp, 1.49072422_wp, 2.55201615_wp, &
      & 2.76083520_wp, 1.50599537_wp, 2.56052955_wp, 2.73276726_wp, 1.52126652_wp, 2.56904294_wp, &
      & 2.70469932_wp, 1.53653767_wp, 2.57755634_wp, 2.67663138_wp, 1.55180881_wp, 2.58606974_wp, &
      & 2.64856345_wp, 1.56707996_wp, 2.59458313_wp, 2.62049551_wp, 1.58235111_wp, 2.60309653_wp, &
      & 2.59242757_wp, 1.59762226_wp, 2.61160992_wp, 2.56435964_wp, 1.61289341_wp, 2.62012332_wp, &
      & 2.53629170_wp, 1.62816456_wp, 2.62863672_wp, 2.50822376_wp, 1.64343571_wp, 2.63715011_wp, &
      & 2.48015583_wp, 1.65870685_wp, 2.64566351_wp, 3.19537752_wp, 2.24853837_wp, 2.41492177_wp, &
      & 3.14122020_wp, 2.48723489_wp, 2.21933576_wp, 3.17661283_wp, 3.39538568_wp, 2.37502789_wp, &
      & 3.14538352_wp, 2.58361113_wp, 2.47139347_wp, 1.81565647_wp, 2.48106221_wp, 3.18585355_wp, &
      & 2.11798490_wp, 2.85857032_wp, 3.47048400_wp, 2.71241232_wp, 3.37886078_wp, 3.64124964_wp, &
      & 2.80572458_wp, 2.82570220_wp, 3.72064445_wp, 2.61951362_wp, 2.69607886_wp, 0.00000000_wp, &
      & 3.05383193_wp, 2.61683803_wp, 3.32179612_wp, 3.02135073_wp, 2.59250246_wp, 4.24674489_wp, &
      & 3.16405210_wp, 2.63238785_wp, 3.04625573_wp, 2.96133467_wp, 2.71388453_wp, 2.31022562_wp, &
      & 2.98240599_wp, 2.95960758_wp, 2.43778345_wp, 3.07936232_wp, 2.68589775_wp, 2.10311395_wp],&
      & shape(zeta))

      integer :: isp, izp, ish, stat
      integer, allocatable :: nshell(:)
      type(cgto_type), allocatable :: cgto(:, :)

      nshell = nsh(mol%num)
      allocate(cgto(maxval(nshell), mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, nshell(isp)
            call slater_to_gauss(ng, pqn(ish, izp), lsh(ish, izp), zeta(ish, izp), &
            & cgto(ish, isp), .true., stat)
         end do
      end do

      call new_basis(bas, mol, nshell, cgto, 1.0_wp)

   end subroutine make_basis

   subroutine test_scaled_selfenergy_mol(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: ref(:)

      type(basis_type) :: bas
      type(tb_hamiltonian) :: h0
      type(erf_ncoord_type) :: ncoord
      type(erf_en_ncoord_type) :: ncoord_en
      type(adjacency_list) :: list
      real(wp), parameter :: cn_cutoff = 30.0_wp
      real(wp), allocatable :: lattr(:, :), cn(:), cn_en(:), rcov(:), en(:)
      real(wp), allocatable :: selfenergy(:)
      real(wp) :: cutoff
      integer :: ii, jj

      call make_basis(bas, mol, 6)

      call check(error, bas%nsh, size(ref, 1))
      if (allocated(error)) return

      call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))

      allocate(cn(mol%nat), cn_en(mol%nat), rcov(mol%nid), en(mol%nid))
      ! test with the standard Pyykko radii and Pauling EN (not as in CEH parametrization)
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)
      call new_erf_ncoord(ncoord, mol, cn_cutoff, rcov)
      call new_erf_en_ncoord(ncoord_en, mol, cn_cutoff, rcov)
      call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
      call get_coordination_number(ncoord , mol, lattr, cn_cutoff, cn)
      call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)

      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      allocate(selfenergy(bas%nsh))
      call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy)

      do ii = 1, size(selfenergy, 1)
         call check(error, selfenergy(ii), ref(ii), thr=thr2)
         if (allocated(error)) then
            print '(2es20.13)', selfenergy(ii), ref(ii)
            return
         end if
      end do

   end subroutine test_scaled_selfenergy_mol


   subroutine test_scaled_selfenergy_numgrad_mol(error, mol)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(inout) :: mol


      type(basis_type) :: bas
      type(tb_hamiltonian) :: h0
      type(erf_ncoord_type) :: ncoord
      type(erf_en_ncoord_type) :: ncoord_en
      type(adjacency_list) :: list      
      integer :: iat, ic
      real(wp), allocatable :: cn(:), cn_en(:), rcov(:), en(:)
      real(wp), allocatable :: selfenergy(:), selfenergyr(:), selfenergyl(:)
      real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :), dcn_endr(:, :, :), dcn_endL(:, :, :)
      real(wp), allocatable :: dsedr(:, :, :), dsedL(:, :, :)
      real(wp), allocatable :: numdr(:, :, :)
      real(wp), allocatable :: lattr(:, :)
      real(wp), parameter :: step = 1.0e-6_wp
      real(wp), parameter :: cn_cutoff = 30.0_wp
      real(wp) :: cutoff

      call make_basis(bas, mol, 6)
      
      allocate(cn(mol%nat), cn_en(mol%nat), rcov(mol%nid), en(mol%nid), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & dcn_endr(3, mol%nat, mol%nat), dcn_endL(3, 3, mol%nat), &
      & numdr(3, mol%nat, bas%nsh))

      ! test with the standard Pyykko radii and Pauling EN (not as in CEH parametrization)
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)

      call new_erf_ncoord(ncoord, mol, cn_cutoff, rcov)
      call new_erf_en_ncoord(ncoord_en, mol, cn_cutoff, rcov)
      cutoff = get_cutoff(bas)

      allocate(selfenergy(bas%nsh), selfenergyr(bas%nsh), selfenergyl(bas%nsh))
      do iat = 1, mol%nat
         do ic = 1, 3
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))
      
            call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
            call get_coordination_number(ncoord , mol, lattr, cn_cutoff, cn)
            call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)

            call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
            & selfenergy=selfenergyr)
            
            mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
            call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))
      
            call get_coordination_number(ncoord , mol, lattr, cn_cutoff, cn)
            call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)
            
            call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
            & selfenergy=selfenergyl)

            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            numdr(ic, iat, :) = 0.5_wp*(selfenergyr - selfenergyl)/step
         end do
      end do

      call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))
      
      call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
      call get_coordination_number(ncoord, mol, lattr, cutoff, cn, dcndr, dcndL)
      call get_coordination_number(ncoord_en, mol, lattr, cutoff, cn_en, dcn_endr, dcn_endL)

      allocate(dsedr(3, mol%nat,bas%nsh), dsedL(3, 3, bas%nsh))
      call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, &
      & cn=cn, cn_en=cn_en, dcndr=dcndr, dcndL=dcndL, dcn_endr=dcn_endr, dcn_endL=dcn_endL, &
      & selfenergy=selfenergy, dsedr=dsedr, dsedL=dsedL)

      if (any(abs(dsedr - numdr) > thr2)) then
         call test_failed(error, "Derivative of selfenergies does not match")
      end if

   end subroutine test_scaled_selfenergy_numgrad_mol



   subroutine test_hamiltonian_mol(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: ref(:, :)

      type(basis_type) :: bas
      type(tb_hamiltonian) :: h0
      type(erf_ncoord_type) :: ncoord
      type(erf_en_ncoord_type) :: ncoord_en
      type(adjacency_list) :: list
      real(wp), parameter :: cn_cutoff = 30.0_wp
      real(wp), allocatable :: lattr(:, :), cn(:), cn_en(:), rcov(:), en(:)
      real(wp), allocatable :: overlap(:, :), overlap_diat(:, :), dpint(:, :, :)
      real(wp), allocatable :: hamiltonian(:, :), selfenergy(:)
      real(wp) :: cutoff
      integer :: ii, jj

      call make_basis(bas, mol, 6)

      call check(error, bas%nao, size(ref, 1))
      if (allocated(error)) return

      call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))

      allocate(cn(mol%nat), cn_en(mol%nat), rcov(mol%nid), en(mol%nid))
      ! test with the standard Pyykko radii and Pauling EN (not as in CEH parametrization)
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)
      call new_erf_ncoord(ncoord, mol, cn_cutoff, rcov)
      call new_erf_en_ncoord(ncoord_en, mol, cn_cutoff, rcov)
      call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
      call get_coordination_number(ncoord , mol, lattr, cn_cutoff, cn)
      call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)

      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      allocate(selfenergy(bas%nsh))
      call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy)

      allocate(overlap(bas%nao, bas%nao), overlap_diat(bas%nao, bas%nao), &
      & dpint(3, bas%nao, bas%nao), hamiltonian(bas%nao, bas%nao))

      call get_hamiltonian(mol, lattr, list, bas, h0, selfenergy, &
      & overlap, overlap_diat, dpint, hamiltonian)

      do ii = 1, size(hamiltonian, 2)
         do jj = 1, size(hamiltonian, 1)
            call check(error, hamiltonian(jj, ii), ref(jj, ii), thr=thr2)
            if (allocated(error)) then
               print '(3es21.13)', hamiltonian(jj, ii), ref(jj, ii), &
               & hamiltonian(jj, ii) - ref(jj, ii)
               return
            end if
         end do
      end do

   end subroutine test_hamiltonian_mol


   subroutine test_hamiltonian_numgrad(error, mol)
   
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      !> Molecular structure data
      type(structure_type), intent(inout) :: mol
   
      type(basis_type) :: bas
      type(tb_hamiltonian) :: h0
      type(erf_ncoord_type) :: ncoord
      type(erf_en_ncoord_type) :: ncoord_en
      type(adjacency_list) :: list
      type(potential_type) :: pot

      real(wp), parameter :: cn_cutoff = 30.0_wp
      real(wp), parameter :: step = 1.0e-6_wp
      real(wp), allocatable :: lattr(:, :), cn_lattr(:, :), cn(:), cn_en(:), rcov(:), en(:)
      real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :), dcn_endr(:, :, :), dcn_endL(:, :, :)
      real(wp), allocatable :: overlapl(:, :), overlapr(:, :), overlap_diatl(:, :)
      real(wp), allocatable :: overlap_diatr(:, :), overlap_diat(:, :), dpint(:, :, :) 
      real(wp), allocatable :: hamiltonian(:,:), hamiltonianl(:, :), hamiltonianr(:, :)
      real(wp), allocatable :: selfenergy(:), dsedr(:,:,:), dsedL(:,:,:)
      real(wp), allocatable :: numdr(:, :, :), dh0dr(:, :, :), dh0dL(:, :, :), doverlap(:, :, :), doverlap_diat(:, :, :)  
      real(wp), allocatable :: dummy_pmat(:, :, :)
      real(wp) :: cutoff
      integer :: iat, ic, ii, jj, is, ish, izp, iao, jat, js, jsh, jzp, jao, kat
   
      call make_basis(bas, mol, 6)
   
      call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))
      
      !> Get initial potential
      call new_potential(pot, mol, bas, 1)
      !> Set potential to zero
      call pot%reset

      allocate(cn(mol%nat), cn_en(mol%nat), rcov(mol%nid), en(mol%nid), source=0.0_wp)
      allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), source=0.0_wp)
      allocate(dcn_endr(3, mol%nat, mol%nat), dcn_endL(3, 3, mol%nat), source=0.0_wp)
      allocate(selfenergy(bas%nsh), dsedr(3, mol%nat,bas%nsh), dsedL(3, 3, bas%nsh), source=0.0_wp)
      
      allocate(overlapr(bas%nao, bas%nao),overlapl(bas%nao, bas%nao), &
        & overlap_diatr(bas%nao, bas%nao), overlap_diatl(bas%nao, bas%nao), dpint(3, bas%nao, bas%nao), &
        & hamiltonianr(bas%nao, bas%nao), hamiltonianl(bas%nao, bas%nao), source=0.0_wp)
      
      allocate(numdr(3, bas%nao, bas%nao), dummy_pmat(bas%nao,bas%nao,1), &
      & doverlap(3,bas%nao,bas%nao), doverlap_diat(3,bas%nao,bas%nao), dh0dr(3, bas%nao, bas%nao), &
      & dh0dL(3, bas%nao, bas%nao), source=0.0_wp)

      ! test with the standard Pyykko radii and Pauling EN (not as in CEH parametrization)
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)

      call new_erf_ncoord(ncoord, mol, cn_cutoff, rcov)
      call new_erf_en_ncoord(ncoord_en, mol, cn_cutoff, rcov)


      do ic = 1, 3
         do iat = 1, mol%nat
            izp = mol%id(iat)
            is = bas%ish_at(iat)
            ! right hand
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            
            ! CN
            call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
            call get_coordination_number(ncoord, mol, lattr, cn_cutoff, cn)
            call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)

            ! Adjacency list
            cutoff = get_cutoff(bas)
            call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
            call new_adjacency_list(list, mol, lattr, cutoff)
            
            ! Self energy
            call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, & 
               &cn=cn, cn_en=cn_en, selfenergy=selfenergy)

            ! Hamiltonian
            call get_hamiltonian(mol, lattr, list, bas, h0, selfenergy, &
               & overlapr, overlap_diatr, dpint, hamiltonianr)
   
            ! left hand
            mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step

           ! CN
            call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
            call get_coordination_number(ncoord, mol, lattr, cn_cutoff, cn)
            call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)

            ! Adjacency list
            cutoff = get_cutoff(bas)
            call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
            call new_adjacency_list(list, mol, lattr, cutoff)
            
            ! Self energy
            call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, & 
               &cn=cn, cn_en=cn_en, selfenergy=selfenergy)

            ! Hamiltonian
            call get_hamiltonian(mol, lattr, list, bas, h0, selfenergy, &
               & overlapl, overlap_diatl, dpint, hamiltonianl)

            ! Geometry reset 
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

            ! Numerical gradient of the hamiltonian matrix
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do iao = 1, bas%nao_sh(is+ish)
                  ! Use only the upper triangular matrix to not sum different elements
                  do jat = 1, iat
                     jzp = mol%id(jat)
                     js = bas%ish_at(jat)
                     do jsh = 1, bas%nsh_id(jzp) 
                        jj = bas%iao_sh(js+jsh)
                        do jao = 1, bas%nao_sh(js+jsh)
                           ! Upper triangular matrix and diagonal
                           numdr(ic, jj+jao, ii+iao) = & 
                              & + 0.5_wp*(hamiltonianr(jj+jao, ii+iao) - hamiltonianl(jj+jao, ii+iao))/step

                           ! Lower triangular matrix
                           if(jat /= iat) then
                              numdr(ic, ii+iao, jj+jao) = &
                                 & - 0.5_wp*(hamiltonianr(ii+iao, jj+jao) - hamiltonianl(ii+iao, jj+jao))/step
                           end if
                        end do
                     end do 
                  end do
               end do
            end do
         end do
      end do

      ! CN
      call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
      call get_coordination_number(ncoord , mol, lattr, cn_cutoff, cn, dcndr, dcndL)
      call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en, dcn_endr, dcn_endL)

      ! Adjacency list
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, &
      & cn=cn, cn_en=cn_en, dcndr=dcndr, dcndL=dcndL, dcn_endr=dcn_endr, dcn_endL=dcn_endL, &
      & selfenergy=selfenergy, dsedr=dsedr, dsedL=dsedL)    
      
      call get_hamiltonian_gradient(mol, lattr, list, bas, h0, selfenergy, &
         & dsedr, dsedL, pot, dummy_pmat, dh0dr, dh0dL, doverlap, doverlap_diat)

      num: do ic = 1, 3
         do ii = 1, size(numdr,2)
            do jj = 1, size(numdr,3)
               call check(error, numdr(ic, ii, jj), dh0dr(ic, ii, jj), thr=thr2)
               if (allocated(error)) then 
                  call test_failed(error, "Hamiltonian derivative of does not match")
                  exit num
               end if
            end do           
         end do
      end do num

      if (any(abs(numdr - dh0dr) > thr2)) then
         call test_failed(error, "Derivative of does not match")
      end if
   
   end subroutine test_hamiltonian_numgrad

  subroutine write_2d_matrix(matrix, name, unit, step)
    implicit none
    real(wp), intent(in) :: matrix(:, :)
    character(len=*), intent(in), optional :: name
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: step
    integer :: d1, d2
    integer :: i, j, k, l, istep, iunit

    d1 = size(matrix, dim=1)
    d2 = size(matrix, dim=2)

    if (present(unit)) then
      iunit = unit
    else
      iunit = output_unit
    end if

    if (present(step)) then
      istep = step
    else
      istep = 6
    end if

    if (present(name)) write (iunit, '(/,"matrix printed:",1x,a)') name

    do i = 1, d2, istep
      l = min(i + istep - 1, d2)
      write (iunit, '(/,6x)', advance='no')
      do k = i, l
        write (iunit, '(6x,i7,3x)', advance='no') k
      end do
      write (iunit, '(a)')
      do j = 1, d1
        write (iunit, '(i6)', advance='no') j
        do k = i, l
          write (iunit, '(1x,f15.8)', advance='no') matrix(j, k)
        end do
        write (iunit, '(a)')
      end do
    end do

  end subroutine write_2d_matrix
! subroutine test_numsigma(error, mol, ncoord, cutoff)

!    !> Error handling
!    type(error_type), allocatable, intent(out) :: error

!    !> Molecular structure data
!    type(structure_type), intent(inout) :: mol

!    !> Coordination number type
!    class(ncoord_type)   :: ncoord

!    real(wp), intent(in) :: cutoff

!    integer :: ic, jc
!    real(wp) :: eps(3, 3)
!    real(wp), allocatable :: cn(:), cnr(:), cnl(:), xyz(:, :)
!    real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :)
!    real(wp), allocatable :: numdL(:, :, :)
!    real(wp), allocatable :: lattr(:, :), trans(:, :)
!    real(wp), parameter :: unity(3, 3) = reshape(&
!       & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
!    real(wp), parameter :: step = 1.0e-6_wp

!    allocate(cn(mol%nat), cnr(mol%nat), cnl(mol%nat), &
!       & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), xyz(3, mol%nat), &
!       & numdL(3, 3, mol%nat))
!    call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

!    eps(:, :) = unity
!    xyz(:, :) = mol%xyz
!    trans = lattr
!    do ic = 1, 3
!       do jc = 1, 3
!          eps(jc, ic) = eps(jc, ic) + step
!          mol%xyz(:, :) = matmul(eps, xyz)
!          lattr(:, :) = matmul(eps, trans)
!          call get_coordination_number(ncoord, mol, lattr, cutoff, cnr)
!          eps(jc, ic) = eps(jc, ic) - 2*step
!          mol%xyz(:, :) = matmul(eps, xyz)
!          lattr(:, :) = matmul(eps, trans)
!          call get_coordination_number(ncoord, mol, lattr, cutoff, cnl)
!          eps(jc, ic) = eps(jc, ic) + step
!          mol%xyz(:, :) = xyz
!          lattr(:, :) = trans
!          numdL(jc, ic, :) = 0.5_wp*(cnr - cnl)/step
!       end do
!    end do

!    call get_coordination_number(ncoord, mol, lattr, cutoff, cn, dcndr, dcndL)

!    if (any(abs(dcndL - numdL) > thr2)) then
!       call test_failed(error, "Derivative of coordination number does not match")
!    end if

! end subroutine test_numsigma






   subroutine test_overlap_diat_mol(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: ref(:, :)

      type(basis_type) :: bas
      type(tb_hamiltonian) :: h0
      type(erf_ncoord_type) :: ncoord
      type(erf_en_ncoord_type) :: ncoord_en
      type(adjacency_list) :: list
      real(wp), parameter :: cn_cutoff = 30.0_wp
      real(wp), allocatable :: lattr(:, :), cn(:), cn_en(:), rcov(:), en(:)
      real(wp), allocatable :: overlap(:, :), overlap_diat(:, :), dpint(:, :, :)
      real(wp), allocatable :: hamiltonian(:, :), selfenergy(:)
      real(wp) :: cutoff
      integer :: ii, jj

      call make_basis(bas, mol, 6)

      call check(error, bas%nao, size(ref, 1))
      if (allocated(error)) return

      call new_hamiltonian(h0, mol, bas, ceh_h0spec(mol))

      allocate(cn(mol%nat), cn_en(mol%nat), rcov(mol%nid), en(mol%nid))
      ! test with the standard Pyykko radii and Pauling EN (not as in CEH parametrization)
      rcov(:) = get_covalent_rad(mol%num)
      en(:) = get_pauling_en(mol%num)
      call new_erf_ncoord(ncoord, mol, cn_cutoff, rcov)
      call new_erf_en_ncoord(ncoord_en, mol, cn_cutoff, rcov)
      call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
      call get_coordination_number(ncoord , mol, lattr, cn_cutoff, cn)
      call get_coordination_number(ncoord_en, mol, lattr, cn_cutoff, cn_en)

      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      allocate(selfenergy(bas%nsh))
      call get_scaled_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy)

      allocate(overlap(bas%nao, bas%nao), overlap_diat(bas%nao, bas%nao), &
      & dpint(3, bas%nao, bas%nao), hamiltonian(bas%nao, bas%nao))

      call get_hamiltonian(mol, lattr, list, bas, h0, selfenergy, &
      & overlap, overlap_diat, dpint, hamiltonian)

      do ii = 1, size(overlap_diat, 2)
         do jj = 1, size(overlap_diat, 1)
            call check(error, overlap_diat(jj, ii), ref(jj, ii), thr=thr2)
            if (allocated(error)) then
               print '(3es21.13)', overlap_diat(jj, ii), ref(jj, ii), &
               & overlap_diat(jj, ii) - ref(jj, ii)
               return
            end if
         end do
      end do

   end subroutine test_overlap_diat_mol

   subroutine test_q_gen(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      !> Reference CEH charges
      real(wp), intent(in) :: ref(:)

      type(context_type) :: ctx
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp), allocatable :: cn(:), cn_en(:)
      real(wp), parameter :: accuracy = 1e-8_wp
      real(wp), allocatable :: lattr(:, :)
      integer :: i
      allocate(cn(mol%nat))

      call new_ceh_calculator(calc, mol)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      ctx%verbosity = 0
      call ceh_guess(ctx, calc, mol, error, wfn, accuracy, .false.)

      do i = 1, mol%nat
         call check(error, wfn%qat(i,1), ref(i), thr=1e-6_wp)
         if (allocated(error)) then
            print '(3es21.13)',  wfn%qat(i,1), ref(i), &
            & wfn%qat(i,1) - ref(i)
            return
         end if
      enddo

   end subroutine test_q_gen


   subroutine test_scaled_selfenergy_h2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 2
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & -5.1041627058615E-01_wp, -5.1041627058615E-01_wp &
      &],shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "H2")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_h2

   subroutine test_scaled_selfenergy_lih(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 3
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & -3.7307764740843E-01_wp, -3.9446056938638E-01_wp, -3.3732801377653E-01_wp &
      &],shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "LiH")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_lih

   subroutine test_scaled_selfenergy_s2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 6
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & -5.7642093144110E-01_wp, -5.3180811904793E-01_wp, -2.9175444046022E-01_wp, & 
      & -5.7642093144110E-01_wp, -5.3180811904793E-01_wp, -2.9175444046022E-01_wp &
      &], shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "S2")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_s2

   subroutine test_scaled_selfenergy_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nsh = 7
      real(wp), parameter :: scaled_selfenergy(nsh) = reshape([&
      & -5.2140420559246E-01_wp, -4.8639401015524E-01_wp, -2.4597945091348E-01_wp, & 
      & -4.7036092200061E-01_wp, -4.7036092200061E-01_wp, -4.7036092200061E-01_wp, & 
      & -4.7036092200061E-01_wp], shape(scaled_selfenergy))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_scaled_selfenergy_mol(error, mol, scaled_selfenergy)

   end subroutine test_scaled_selfenergy_sih4


   subroutine test_scaled_selfenergy_numgrad_h2(error)
   
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol
   
      call get_structure(mol, "MB16-43", "H2")
      call test_scaled_selfenergy_numgrad_mol(error, mol)
   
   end subroutine test_scaled_selfenergy_numgrad_h2
   
   subroutine test_scaled_selfenergy_numgrad_lih(error)
   
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol
   
      call get_structure(mol, "MB16-43", "LiH")
      call test_scaled_selfenergy_numgrad_mol(error, mol)
   
   end subroutine test_scaled_selfenergy_numgrad_lih
   
   subroutine test_scaled_selfenergy_numgrad_s2(error)
   
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol
   
      call get_structure(mol, "MB16-43", "S2")
      call test_scaled_selfenergy_numgrad_mol(error, mol)
   
   end subroutine test_scaled_selfenergy_numgrad_s2

   subroutine test_scaled_selfenergy_numgrad_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_scaled_selfenergy_numgrad_mol(error, mol)

   end subroutine test_scaled_selfenergy_numgrad_sih4


   subroutine test_hamiltonian_h2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 2
      real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -5.1041627058615E-01_wp, -3.9044476391402E-01_wp, -3.9044476391402E-01_wp, &
      & -5.1041627058615E-01_wp],shape(hamiltonian))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "H2")
      call test_hamiltonian_mol(error, mol, hamiltonian)

   end subroutine test_hamiltonian_h2

   subroutine test_hamiltonian_lih(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 5
      real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -3.7307764740843E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -1.7831565052851E-01_wp, +0.0000000000000E+00_wp, &
      & -3.9446056938638E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -3.9446056938638E-01_wp, +0.0000000000000E+00_wp, -1.8525661956590E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -3.9446056938638E-01_wp, +0.0000000000000E+00_wp, -1.7831565052851E-01_wp, &
      & +0.0000000000000E+00_wp, -1.8525661956590E-01_wp, +0.0000000000000E+00_wp, &
      & -3.3732801377653E-01_wp],shape(hamiltonian))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "LiH")
      call test_hamiltonian_mol(error, mol, hamiltonian)

   end subroutine test_hamiltonian_lih

   subroutine test_hamiltonian_s2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 18
      real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -5.7642093144110E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -6.8813732574872E-02_wp, +0.0000000000000E+00_wp, +1.5327215355447E-01_wp, &
      & +0.0000000000000E+00_wp, -1.0104862230667E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -5.3180811904793E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -8.8345872545035E-02_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.1751194343267E-02_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -5.3180811904793E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -1.5327215355447E-01_wp, +0.0000000000000E+00_wp, +2.3803538679166E-01_wp, &
      & +0.0000000000000E+00_wp, -1.3983152997378E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -5.3180811904793E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -8.8345872545035E-02_wp, +0.0000000000000E+00_wp, +9.1751194343267E-02_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -1.0104862230667E-01_wp, +0.0000000000000E+00_wp, +1.3983152997378E-01_wp, &
      & +0.0000000000000E+00_wp, -9.0501662901705E-02_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -9.1751194343267E-02_wp, +0.0000000000000E+00_wp, +9.2973802728680E-02_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -2.9175444046022E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -9.1751194343267E-02_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.2973802728680E-02_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -3.4662194100260E-02_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -3.4662194100260E-02_wp, &
      & -6.8813732574872E-02_wp, +0.0000000000000E+00_wp, -1.5327215355447E-01_wp, &
      & +0.0000000000000E+00_wp, -1.0104862230667E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -5.7642093144110E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -8.8345872545035E-02_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -9.1751194343267E-02_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -5.3180811904793E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +1.5327215355447E-01_wp, +0.0000000000000E+00_wp, +2.3803538679166E-01_wp, &
      & +0.0000000000000E+00_wp, +1.3983152997378E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -5.3180811904793E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -8.8345872545035E-02_wp, +0.0000000000000E+00_wp, -9.1751194343267E-02_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -5.3180811904793E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -1.0104862230667E-01_wp, +0.0000000000000E+00_wp, -1.3983152997378E-01_wp, &
      & +0.0000000000000E+00_wp, -9.0501662901705E-02_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.1751194343267E-02_wp, +0.0000000000000E+00_wp, +9.2973802728680E-02_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +9.1751194343267E-02_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.2973802728680E-02_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -2.9175444046022E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -3.4662194100260E-02_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -2.9175444046022E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -3.4662194100260E-02_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -2.9175444046022E-01_wp],&
      & shape(hamiltonian))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "S2")
      call test_hamiltonian_mol(error, mol, hamiltonian)

   end subroutine test_hamiltonian_s2

   subroutine test_hamiltonian_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 13
      real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -5.2140420559246E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -2.1916633126570E-01_wp, -2.1916633126570E-01_wp, -2.1916633126570E-01_wp, &
      & -2.1916633126570E-01_wp, +0.0000000000000E+00_wp, -4.8639401015524E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -1.6024153492284E-01_wp, +1.6024153492284E-01_wp, &
      & +1.6024153492284E-01_wp, -1.6024153492284E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -4.8639401015524E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +1.6024153492284E-01_wp, &
      & +1.6024153492284E-01_wp, -1.6024153492284E-01_wp, -1.6024153492284E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -4.8639401015524E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -1.6024153492284E-01_wp, +1.6024153492284E-01_wp, -1.6024153492284E-01_wp, &
      & +1.6024153492284E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -2.4597945091348E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -2.0759708971706E-17_wp, -2.0759708971706E-17_wp, &
      & -2.0759708971706E-17_wp, -2.0759708971706E-17_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -2.4597945091348E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +1.2688631551999E-01_wp, &
      & -1.2688631551999E-01_wp, -1.2688631551999E-01_wp, +1.2688631551999E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -2.4597945091348E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +1.2688631551999E-01_wp, -1.2688631551999E-01_wp, +1.2688631551999E-01_wp, &
      & -1.2688631551999E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -2.4597945091348E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -2.4597945091348E-01_wp, -1.2688631551999E-01_wp, &
      & -1.2688631551999E-01_wp, +1.2688631551999E-01_wp, +1.2688631551999E-01_wp, &
      & -2.1916633126570E-01_wp, -1.6024153492284E-01_wp, +1.6024153492284E-01_wp, &
      & -1.6024153492284E-01_wp, -2.0759708971706E-17_wp, +1.2688631551999E-01_wp, &
      & +1.2688631551999E-01_wp, +0.0000000000000E+00_wp, -1.2688631551999E-01_wp, &
      & -4.7036092200061E-01_wp, -3.3290465798408E-02_wp, -3.3290465798408E-02_wp, &
      & -3.3290465798408E-02_wp, -2.1916633126570E-01_wp, +1.6024153492284E-01_wp, &
      & +1.6024153492284E-01_wp, +1.6024153492284E-01_wp, -2.0759708971706E-17_wp, &
      & -1.2688631551999E-01_wp, -1.2688631551999E-01_wp, +0.0000000000000E+00_wp, &
      & -1.2688631551999E-01_wp, -3.3290465798408E-02_wp, -4.7036092200061E-01_wp, &
      & -3.3290465798408E-02_wp, -3.3290465798408E-02_wp, -2.1916633126570E-01_wp, &
      & +1.6024153492284E-01_wp, -1.6024153492284E-01_wp, -1.6024153492284E-01_wp, &
      & -2.0759708971706E-17_wp, -1.2688631551999E-01_wp, +1.2688631551999E-01_wp, &
      & +0.0000000000000E+00_wp, +1.2688631551999E-01_wp, -3.3290465798408E-02_wp, &
      & -3.3290465798408E-02_wp, -4.7036092200061E-01_wp, -3.3290465798408E-02_wp, &
      & -2.1916633126570E-01_wp, -1.6024153492284E-01_wp, -1.6024153492284E-01_wp, &
      & +1.6024153492284E-01_wp, -2.0759708971706E-17_wp, +1.2688631551999E-01_wp, &
      & -1.2688631551999E-01_wp, +0.0000000000000E+00_wp, +1.2688631551999E-01_wp, &
      & -3.3290465798408E-02_wp, -3.3290465798408E-02_wp, -3.3290465798408E-02_wp, &
      & -4.7036092200061E-01_wp],shape(hamiltonian))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_hamiltonian_mol(error, mol, hamiltonian)

   end subroutine test_hamiltonian_sih4


   subroutine test_hamiltonian_numgrad_h2(error)
   
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol
   
      call get_structure(mol, "MB16-43", "H2")
      call test_hamiltonian_numgrad(error, mol)
   
   end subroutine test_hamiltonian_numgrad_h2
   
   subroutine test_hamiltonian_numgrad_lih(error)
   
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol
   
      call get_structure(mol, "MB16-43", "LiH")
      call test_hamiltonian_numgrad(error, mol)
   
   end subroutine test_hamiltonian_numgrad_lih
   
   subroutine test_hamiltonian_numgrad_s2(error)
   
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol
   
      call get_structure(mol, "MB16-43", "S2")
      call test_hamiltonian_numgrad(error, mol)
   
   end subroutine test_hamiltonian_numgrad_s2

   subroutine test_hamiltonian_numgrad_pcl(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "PCl")
      call test_hamiltonian_numgrad(error, mol)

   end subroutine test_hamiltonian_numgrad_pcl



   subroutine test_hamiltonian_numgrad_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_hamiltonian_numgrad(error, mol)

   end subroutine test_hamiltonian_numgrad_sih4

   subroutine test_hamiltonian_numgrad_cecl3(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "CeCl3")
      call test_hamiltonian_numgrad(error, mol)

   end subroutine test_hamiltonian_numgrad_cecl3

   subroutine test_overlap_diat_h2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 2
      real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & +9.9999999988149E-01_wp, +1.1600752205819E+00_wp, +1.1600752205819E+00_wp, &
      & +9.9999999988149E-01_wp],shape(overlap_diat))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "H2")
      call test_overlap_diat_mol(error, mol, overlap_diat)

   end subroutine test_overlap_diat_h2

   subroutine test_overlap_diat_lih(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 5
      real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & +1.0000000000060E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +7.6131450009737E-01_wp, +0.0000000000000E+00_wp, &
      & +9.9999999992569E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.9999999992569E-01_wp, +0.0000000000000E+00_wp, +6.9022133546973E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.9999999992569E-01_wp, +0.0000000000000E+00_wp, +7.6131450009737E-01_wp, &
      & +0.0000000000000E+00_wp, +6.9022133546973E-01_wp, +0.0000000000000E+00_wp, &
      & +9.9999999988149E-01_wp],shape(overlap_diat))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "LiH")
      call test_overlap_diat_mol(error, mol, overlap_diat)

   end subroutine test_overlap_diat_lih

   subroutine test_overlap_diat_s2(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 18
      real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & +9.9999999986933E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +1.8104495601005E-01_wp, +0.0000000000000E+00_wp, -3.7708043848244E-01_wp, &
      & +0.0000000000000E+00_wp, +2.5954275281405E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +9.9999999999806E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +2.0567486368782E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -2.2945817477562E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999999806E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +3.7708043848244E-01_wp, +0.0000000000000E+00_wp, -5.5416166393395E-01_wp, &
      & +0.0000000000000E+00_wp, +3.4970103521297E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.9999999999806E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +2.0567486368782E-01_wp, +0.0000000000000E+00_wp, -2.2945817477562E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +2.5954275281405E-01_wp, +0.0000000000000E+00_wp, -3.4970103521297E-01_wp, &
      & +0.0000000000000E+00_wp, +2.7344681334107E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +2.2945817477562E-01_wp, +0.0000000000000E+00_wp, -2.8091627562657E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +2.2945817477562E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -2.8091627562657E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +1.0473030236384E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +1.0473030236384E-01_wp, &
      & +1.8104495601005E-01_wp, +0.0000000000000E+00_wp, +3.7708043848244E-01_wp, &
      & +0.0000000000000E+00_wp, +2.5954275281405E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.9999999986933E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +2.0567486368782E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +2.2945817477562E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +9.9999999999806E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -3.7708043848244E-01_wp, +0.0000000000000E+00_wp, -5.5416166393395E-01_wp, &
      & +0.0000000000000E+00_wp, -3.4970103521297E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999999806E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +2.0567486368782E-01_wp, +0.0000000000000E+00_wp, +2.2945817477562E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.9999999999806E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +2.5954275281405E-01_wp, +0.0000000000000E+00_wp, +3.4970103521297E-01_wp, &
      & +0.0000000000000E+00_wp, +2.7344681334107E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -2.2945817477562E-01_wp, +0.0000000000000E+00_wp, -2.8091627562657E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, -2.2945817477562E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -2.8091627562657E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +1.0473030236384E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +1.0473030236384E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999983021E-01_wp],&
      & shape(overlap_diat))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "S2")
      call test_overlap_diat_mol(error, mol, overlap_diat)

   end subroutine test_overlap_diat_s2

   subroutine test_overlap_diat_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 13
      real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & +9.9999999986933E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +6.7026425174062E-01_wp, +6.7026425174062E-01_wp, +6.7026425174062E-01_wp, &
      & +6.7026425174062E-01_wp, +0.0000000000000E+00_wp, +9.9999999999806E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +4.5664074590110E-01_wp, -4.5664074590110E-01_wp, &
      & -4.5664074590110E-01_wp, +4.5664074590110E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +9.9999999999806E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -4.5664074590110E-01_wp, &
      & -4.5664074590110E-01_wp, +4.5664074590110E-01_wp, +4.5664074590110E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.9999999999806E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +4.5664074590110E-01_wp, -4.5664074590110E-01_wp, +4.5664074590110E-01_wp, &
      & -4.5664074590110E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +6.4623102437481E-17_wp, +6.4623102437481E-17_wp, &
      & +6.4623102437481E-17_wp, +6.4623102437481E-17_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, -3.9498566077871E-01_wp, &
      & +3.9498566077871E-01_wp, +3.9498566077871E-01_wp, -3.9498566077871E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +9.9999999983021E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & -3.9498566077871E-01_wp, +3.9498566077871E-01_wp, -3.9498566077871E-01_wp, &
      & +3.9498566077871E-01_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, +0.0000000000000E+00_wp, &
      & +0.0000000000000E+00_wp, +9.9999999983021E-01_wp, +3.9498566077871E-01_wp, &
      & +3.9498566077871E-01_wp, -3.9498566077871E-01_wp, -3.9498566077871E-01_wp, &
      & +6.7026425174062E-01_wp, +4.5664074590110E-01_wp, -4.5664074590110E-01_wp, &
      & +4.5664074590110E-01_wp, +6.4623102437481E-17_wp, -3.9498566077871E-01_wp, &
      & -3.9498566077871E-01_wp, +0.0000000000000E+00_wp, +3.9498566077871E-01_wp, &
      & +9.9999999988149E-01_wp, +1.0733458859608E-01_wp, +1.0733458859608E-01_wp, &
      & +1.0733458859608E-01_wp, +6.7026425174062E-01_wp, -4.5664074590110E-01_wp, &
      & -4.5664074590110E-01_wp, -4.5664074590110E-01_wp, +6.4623102437481E-17_wp, &
      & +3.9498566077871E-01_wp, +3.9498566077871E-01_wp, +0.0000000000000E+00_wp, &
      & +3.9498566077871E-01_wp, +1.0733458859608E-01_wp, +9.9999999988149E-01_wp, &
      & +1.0733458859608E-01_wp, +1.0733458859608E-01_wp, +6.7026425174062E-01_wp, &
      & -4.5664074590110E-01_wp, +4.5664074590110E-01_wp, +4.5664074590110E-01_wp, &
      & +6.4623102437481E-17_wp, +3.9498566077871E-01_wp, -3.9498566077871E-01_wp, &
      & +0.0000000000000E+00_wp, -3.9498566077871E-01_wp, +1.0733458859608E-01_wp, &
      & +1.0733458859608E-01_wp, +9.9999999988149E-01_wp, +1.0733458859608E-01_wp, &
      & +6.7026425174062E-01_wp, +4.5664074590110E-01_wp, +4.5664074590110E-01_wp, &
      & -4.5664074590110E-01_wp, +6.4623102437481E-17_wp, -3.9498566077871E-01_wp, &
      & +3.9498566077871E-01_wp, +0.0000000000000E+00_wp, -3.9498566077871E-01_wp, &
      & +1.0733458859608E-01_wp, +1.0733458859608E-01_wp, +1.0733458859608E-01_wp, &
      & +9.9999999988149E-01_wp],shape(overlap_diat))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_overlap_diat_mol(error, mol, overlap_diat)

   end subroutine test_overlap_diat_sih4


   subroutine test_q_h2(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(2) = reshape([ & 
      & 0.0000000000000_wp, 0.00000000000000_wp &
      &], shape(charges))

      call get_structure(mol, "MB16-43", "H2")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_h2


   subroutine test_q_lih(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(2) = reshape([ &
      & 0.383963042442916_wp, -0.383963042442917_wp &
      &], shape(charges))

      call get_structure(mol, "MB16-43", "LiH")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_lih

   subroutine test_q_mb01(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(16) = reshape([ &
      & 0.504694175287870_wp, -0.049320603215405_wp, -0.442887082747601_wp, &     
      &-0.040897715947291_wp, -0.225355610004319_wp,  0.079338661074347_wp, &
      &-0.012184763492912_wp, -0.354104643499813_wp, -0.225243243808309_wp, &     
      & 0.087276255050057_wp,  0.085471100209451_wp, -0.038779744241030_wp, &
      & 0.038616630415300_wp,  0.127389947162501_wp, -0.039323802393347_wp, &
      & 0.505310440150523_wp], shape(charges))

      call get_structure(mol, "MB16-43", "01")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_mb01

   subroutine test_q_mb02(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(16) = reshape([ &
      &-0.059404506571487_wp, -0.118702191502332_wp, -0.195235491894287_wp, &
      &-0.210021399475394_wp,  0.496594666771036_wp,  0.158190196939209_wp, &
      &-0.075318639484230_wp, -0.045052163348822_wp,  0.286946319499636_wp, &
      & 0.157881990181518_wp, -0.029236175252485_wp,  0.281560060657950_wp, &
      &-0.312312030061058_wp, -0.064265894380938_wp,  0.066269279024472_wp, &
      &-0.337894021272329_wp], shape(charges))

      call get_structure(mol, "MB16-43", "02")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_mb02

   subroutine test_q_mb03(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(16) = reshape([ &
      & 0.047060647546981_wp, -0.452721646051320_wp,  0.000750024797478_wp, &
      & 0.303998197523052_wp,  0.479998103945730_wp,  0.099926859160474_wp, &
      &-0.229477809903280_wp,  0.016858971022812_wp,  0.008564268752175_wp, &
      & 0.004606322535074_wp, -0.354502824466586_wp, -0.187631548264506_wp, &     
      &-0.300801554641486_wp,  0.041952785802619_wp,  0.482401552816662_wp, &     
      & 0.039017649441556_wp], shape(charges))

      call get_structure(mol, "MB16-43", "03")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_mb03

   subroutine test_q_mb04(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: charges(16) = reshape([ &
      & 0.017867676479831_wp, -0.103259042275922_wp, -0.035078199415400_wp, &
      &-0.183433026064262_wp,  0.263017445728753_wp,  0.012000420842077_wp, &
      & 0.040327964255115_wp, -0.119339642424892_wp, -0.010898440812181_wp, &
      &-0.056108839369092_wp, -0.137112266394061_wp,  0.334034245022961_wp, &
      & 0.147772637227760_wp, -0.242482916232017_wp,  0.033049020618869_wp, &
      & 0.039642962813086_wp], shape(charges))

      call get_structure(mol, "MB16-43", "04")
      call test_q_gen(error, mol, charges)

   end subroutine test_q_mb04


   subroutine test_q_ef_chrg_mb01(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp), parameter :: accuracy = 1e-8_wp
      class(container_type), allocatable :: cont      
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: ref(16) = reshape([ &
      &-6.90111445726635_wp, -0.99979611590118_wp,   4.2837802019751_wp, &
      &-0.99980021634224_wp,  7.00023420054907_wp,   0.9454454365410_wp, &
      &-0.98804624196187_wp,  5.84908677934913_wp,   2.2839742642425_wp, & 
      & 0.95737100853797_wp,  0.99908251026128_wp, -10.9475535476315_wp, &
      &-4.99687086866152_wp,  2.83244134467160_wp,   4.8453679895474_wp, &
      &-2.16360228791069_wp], shape(ref))

      real(wp) :: efield(3)
      integer :: i

      efield = 0.0_wp
      efield(3) = 0.2_wp

      call get_structure(mol, "MB16-43", "01")
      mol%charge = 2.0_wp
      call new_ceh_calculator(calc, mol)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      cont = electric_field(efield)
      call calc%push_back(cont)
      ctx%verbosity = 0
      call ceh_guess(ctx, calc, mol, error, wfn, accuracy, .false.)
      do i = 1, mol%nat
         call check(error, wfn%qat(i,1), ref(i), thr=5e-6_wp, message="Calculated charge&
         & does not match reference")
         if (allocated(error)) return
      enddo

   end subroutine test_q_ef_chrg_mb01

   subroutine test_d_mb01(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      real(wp) :: dipole(3), tmp(3)
      real(wp), parameter :: accuracy = 1e-8_wp
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: ref(3) = reshape([ &
      & 0.438025381623586_wp, -0.735884148841272_wp, -2.76541717331434_wp &
      &], shape(ref))

      call get_structure(mol, "MB16-43", "01")

      call new_ceh_calculator(calc, mol)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
      ctx%verbosity = 0
      call ceh_guess(ctx, calc, mol, error, wfn, accuracy, .false.)
      tmp = 0.0_wp
      dipole = 0.0_wp
      call gemv(mol%xyz, wfn%qat(:, 1), tmp)
      dipole(:) = tmp + sum(wfn%dpat(:, :, 1), 2)

      if (any(abs(dipole - ref) > 1e-5_wp)) then
         call test_failed(error, "Numerical dipole moment does not match")
         print '(3es21.14)', dipole
         print '("---")'
         print '(3es21.14)', ref
         print '("---")'
         print '(3es21.14)', dipole - ref
      end if

   end subroutine test_d_mb01

   subroutine test_d_field_mb04(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
      class(container_type), allocatable :: cont
      real(wp) :: energy, efield(3), dipole(3), tmp(3)
      real(wp), parameter :: accuracy = 1e-8_wp
      ! calculated with GP3 standalone (full matrix diagonalization)
      real(wp), parameter :: ref(3) = reshape([ & 
         & -20.835631606789_wp,  97.0349021135889_wp, -7.70521258527074_wp &
         ], shape(ref))

      call get_structure(mol, "MB16-43", "04")
      energy = 0.0_wp
      efield = 0.0_wp
      efield(2) = 0.2_wp

      call new_ceh_calculator(calc, mol)
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)

      cont = electric_field(efield)
      call calc%push_back(cont)

      ctx%verbosity = 0
      call ceh_guess(ctx, calc, mol, error, wfn, accuracy, .false.)
      tmp = 0.0_wp
      dipole = 0.0_wp
      call gemv(mol%xyz, wfn%qat(:, 1), tmp)
      dipole(:) = tmp + sum(wfn%dpat(:, :, 1), 2)
      
      if (any(abs(dipole - ref) > 1e-5_wp)) then
         call test_failed(error, "Numerical dipole moment does not match")
         print '(3es21.14)', dipole
         print '("---")'
         print '(3es21.14)', ref
         print '("---")'
         print '(3es21.14)', dipole - ref
      end if

   end subroutine test_d_field_mb04

   subroutine test_d_hcn(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(context_type) :: ctx
      type(structure_type) :: mol1,mol2
      type(xtb_calculator) :: calc1,calc2
      type(wavefunction_type) :: wfn1,wfn2
      class(container_type), allocatable :: cont1,cont2
      real(wp), parameter :: accuracy = 1e-8_wp
      real(wp) :: efield(3), dip1(3), dip2(3), tmp(3)
      integer, parameter :: num(3) = reshape([ &
         7, &
         6, &
         1], shape(num))
      integer, parameter :: nat = 3
      real(wp) :: xyz(3, nat) = reshape([ &
      & -0.09604091224796_wp,  0.0_wp, 0.0_wp, &
      &  2.09604091224796_wp,  0.0_wp, 0.0_wp, &
      &  4.10859879422050_wp,  0.0_wp, 0.0_wp], &
      & shape(xyz))

      ctx%verbosity = 0
      call new(mol1, num, xyz)
      efield = 0.0_wp
      efield(1) = -0.1_wp
      call new_ceh_calculator(calc1, mol1)
      call new_wavefunction(wfn1, mol1%nat, calc1%bas%nsh, calc1%bas%nao, 1, kt)
      cont1 = electric_field(efield)
      call calc1%push_back(cont1)
      call ceh_guess(ctx, calc1, mol1, error, wfn1, accuracy, .false.)
      tmp = 0.0_wp
      dip1 = 0.0_wp
      call gemv(mol1%xyz, wfn1%qat(:, 1), tmp)
      dip1(:) = tmp + sum(wfn1%dpat(:, :, 1), 2)

      xyz(1, :) = xyz(1, :) - 1.0_wp
      call new(mol2, num, xyz)
      call new_ceh_calculator(calc2, mol2)
      call new_wavefunction(wfn2, mol2%nat, calc2%bas%nsh, calc2%bas%nao, 1, kt)
      cont2 = electric_field(efield)
      call calc2%push_back(cont2)
      call ceh_guess(ctx, calc2, mol2, error, wfn2, accuracy, .false.)
      tmp = 0.0_wp
      dip2 = 0.0_wp
      call gemv(mol2%xyz, wfn2%qat(:, 1), tmp)
      dip2(:) = tmp + sum(wfn2%dpat(:, :, 1), 2)

      if (any(abs(dip1 - dip2) > 1e-7_wp)) then
         call test_failed(error, "Numerical dipole moment does not match")
         print '(3es21.14)', dip1
         print '("---")'
         print '(3es21.14)', dip2
         print '("---")'
         print '(3es21.14)', dip1 - dip2
      end if

   end subroutine test_d_hcn
end module test_ceh

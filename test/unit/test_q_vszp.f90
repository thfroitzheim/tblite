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

module test_qvszp

   
   use iso_fortran_env, only: output_unit

   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_qvszp, only : qvszp_basis_type, qvszp_cgto_type, new_basis, &
   & add_qvszp_cgtos
   use tblite_basis_type, only : get_cutoff, maxg
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_overlap, only : overlap_cgto, msao, &
   & get_overlap, get_overlap_gradient
   use tblite_context_type, only : context_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_ceh_singlepoint, only : ceh_singlepoint
   use tblite_ceh_ceh, only : new_ceh_calculator
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction

   implicit none
   private

   public :: collect_qvszp

   ! CEH temperature
   real(wp), parameter :: kt = 4000.0_wp * 3.166808578545117e-06_wp
   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains

   !> Collect all exported unit tests
   subroutine collect_qvszp(testsuite)

      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("norm-qvszp-mb01", test_qvszp_norm_cgto_mb01), &
         new_unittest("norm-qvszp-accl3", test_qvszp_norm_cgto_accl3), &
         new_unittest("overlap-qvszp-lih", test_qvszp_overlap_lih), &
         new_unittest("overlap-qvszp-sih4", test_qvszp_overlap_sih4), &
         new_unittest("overlap-qvszp-accl3", test_qvszp_overlap_accl3), &
         new_unittest("coefficient-grad-qvszp-lih", test_qvszp_coefficient_grad_lih), &
         new_unittest("coefficient-grad-qvszp-sih4", test_qvszp_coefficient_grad_sih4), &
         new_unittest("coefficient-grad-qvszp-accl3", test_qvszp_coefficient_grad_accl3) &
         !new_unittest("overlap-diat-grad-qvszp-lih", test_qvszp_overlap_diat_grad_lih) &
         & ]

   end subroutine collect_qvszp

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
            write (iunit, '(1x,f15.12)', advance='no') matrix(j, k)
          end do
          write (iunit, '(a)')
        end do
      end do
  
   end subroutine write_2d_matrix

   subroutine test_qvszp_norm_cgto(error, mol)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      type(qvszp_cgto_type), allocatable :: cgto(:, :)
      type(qvszp_basis_type) :: bas
      type(context_type) :: ctx
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
            
      integer, allocatable :: nsh_id(:)
      real(wp), parameter :: accuracy = 1e-8_wp
      real(wp) :: vec(3), r2

      real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
      real(wp), allocatable :: overlap(:,:)

      integer :: iat, isp, ish, iao, l

      ! Setup q-vSZP CGTOs and basis 
      call add_qvszp_cgtos(cgto, mol, nsh_id, .false.)
      call new_basis(bas, mol, nsh_id, cgto, 1.0_wp)
      
      ! Obtain CEH charges for charge adaptation
      call new_ceh_calculator(calc, mol, error)
      if (allocated(error)) return
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .false.)
      ctx%verbosity = 0
      call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .false.)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if

      ! Get q-vSZP specific coordination number for scaling
      allocate(cn(mol%nat))
      call bas%ncoord%get_cn(mol, cn)

      ! Scale the basis set 
      call bas%scale_basis(mol, wfn%qat(:,1), cn, .true.)

      ! Check for normalization
      vec = 0.0
      r2 = 0.0
      do iat = 1, mol%nat
         isp = mol%id(iat)
         do ish = 1, nsh_id(isp)
            l = bas%cgto(ish,iat)%ang
            allocate(overlap(msao(l), (msao(l))), source=0.0_wp)   
            call overlap_cgto(bas%cgto(ish,iat), bas%cgto(ish,iat), r2, vec, 100.0_wp, overlap)
            do iao = 1, msao(l)
               call check(error, overlap(iao, iao), 1.0_wp, thr=thr)
               if (allocated(error)) then
                  print*, (overlap(iao, iao) - 1.0_wp)
               end if
            end do
            deallocate(overlap)
         end do
      end do

   end subroutine test_qvszp_norm_cgto

   subroutine test_qvszp_overlap(error, mol, ref)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      real(wp), intent(in) :: ref(:, :)

      type(qvszp_cgto_type), allocatable :: cgto(:, :)
      type(qvszp_basis_type) :: bas
      type(adjacency_list) :: list
      type(context_type) :: ctx
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
            
      integer, allocatable :: nsh_id(:)
      real(wp), parameter :: accuracy = 1e-8_wp

      real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
      real(wp), allocatable :: overlap(:,:), lattr(:, :)
      real(wp) :: cutoff

      integer :: ii, jj

      ! Setup q-vSZP CGTOs and basis 
      call add_qvszp_cgtos(cgto, mol, nsh_id, .false.)
      call new_basis(bas, mol, nsh_id, cgto, 1.0_wp)
      
      ! Obtain CEH charges for charge adaptation
      call new_ceh_calculator(calc, mol, error)
      if (allocated(error)) return
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .false.)
      ctx%verbosity = 0
      call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .false.)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if

      ! Get q-vSZP specific coordination number for scaling
      allocate(cn(mol%nat))
      call bas%ncoord%get_cn(mol, cn)

      ! Scale the basis set 
      call bas%scale_basis(mol, wfn%qat(:,1), cn, .true.)

      ! Get overlap matrix
      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      allocate(overlap(bas%nao, bas%nao))
      call get_overlap(mol, lattr, list, bas, overlap)

      !print '(*(6x,"&", 3(es22.15e2, "_wp":, ","), "&", /))', overlap
      !call write_2d_matrix(overlap, "overlap")
      ! Check overlap matrix
      do ii = 1, size(overlap, 2)
         do jj = 1, size(overlap, 1)
            call check(error, overlap(jj, ii), ref(jj, ii), thr=thr)
            if (allocated(error)) return
         end do
      end do

   end subroutine test_qvszp_overlap

   subroutine test_qvszp_coefficient_grad(mol, error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Molecular structure data
      type(structure_type), intent(inout) :: mol

      type(qvszp_cgto_type), allocatable :: cgtor(:, :), cgtol(:, :)
      type(qvszp_basis_type) :: basr, basl
      type(context_type) :: ctx
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
            
      integer, allocatable :: nsh_id(:)
      real(wp), parameter :: accuracy = 1e-8_wp, cn_cutoff = 30.0_wp, step = 1.0e-6_wp
      real(wp) :: vec(3), r2

      real(wp), allocatable :: lattr(:, :), cn(:), dcndr(:, :, :), dcndL(:, :, :)
      real(wp), allocatable :: ql(:), qr(:)
      real(wp), allocatable :: numdr(:, :, :, :, :), dqatnumdr(:, :, :)

      integer :: ic, iat, izp, jat, jzp, jsh, jpr

      allocate(ql(mol%nat), qr(mol%nat), dqatnumdr(3, mol%nat, mol%nat))
      allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), source=0.0_wp)

      ! Setup right q-vSZP CGTOs and basis 
      call add_qvszp_cgtos(cgtor, mol, nsh_id, .false.)
      call new_basis(basr, mol, nsh_id, cgtor, 1.0_wp)
      ! Setup left q-vSZP CGTOs and basis 
      call add_qvszp_cgtos(cgtol, mol, nsh_id, .false.)
      call new_basis(basl, mol, nsh_id, cgtol, 1.0_wp)

      allocate(numdr(3, mol%nat, mol%nat, maxval(basr%nsh_id), maxg), source=0.0_wp)
      
      ! Setup the CEH 
      call new_ceh_calculator(calc, mol, error)
      if (allocated(error)) return
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .false.)
      ctx%verbosity = 0

      do ic = 1, 3
         do iat = 1, mol%nat
            izp = mol%id(iat)

            ! Right hand
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            
            ! Obtain CEH charges for charge adaptation
            call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .false.)
            if (ctx%failed()) then
               call ctx%get_error(error)
               return
            end if
            qr = wfn%qat(:,1)

            ! Obtain CN for charge adaptation
            call basr%ncoord%get_cn(mol, cn)

            ! Scale the basis set 
            call basr%scale_basis(mol, wfn%qat(:,1), cn, .true.)

            ! Left hand
            mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step

            ! Obtain CEH charges for charge adaptation
            call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .false.)
            if (ctx%failed()) then
               call ctx%get_error(error)
               return
            end if
            ql = wfn%qat(:,1)

            ! Obtain CN for charge adaptation
            call basl%ncoord%get_cn(mol, cn)

            ! Scale the basis set 
            call basl%scale_basis(mol, wfn%qat(:,1), cn, .true.)

            ! Geometry reset 
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

            ! Numerical gradient of the hamiltonian matrix
            do jat = 1, mol%nat
               jzp = mol%id(jat)
               do jsh = 1, basr%nsh_id(jzp)
                  ! Upper triangular matrix and diagonal
                  numdr(ic, iat, jat, jsh, :) = & 
                     & + 0.5_wp*(basr%cgto(jsh,jat)%coeff - basl%cgto(jsh,jat)%coeff)/step
               end do 
            end do
            ! CEH gradient (intermediate solution)
            dqatnumdr(ic, iat, :) = 0.5_wp*(qr - ql)/step

         end do
      end do

      ! Obtain CEH charges and derivatives for charge adaptation
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .true.)
      call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .true.)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
      wfn%dqatdr(:, :, :, 1) = dqatnumdr
      wfn%dqatdL(:, :, :, 1) = 0.0_wp 

      ! Obtain CN for charge adaptation
      call basr%ncoord%get_cn(mol, cn, dcndr, dcndL)

      ! Scale the basis set 
      call basr%scale_basis(mol, wfn%qat(:,1), cn, .true., dqatdr=wfn%dqatdr(:, :, :, 1), &
      & dqatdL=wfn%dqatdL(:, :, :, 1), dcndr=dcndr(:, :, :), dcndL=dcndL(:, :, :))

      ! Check results
      num: do ic = 1, 3
         do iat = 1, mol%nat
            do jat = 1, mol%nat
               jzp = mol%id(jat)
               do jsh = 1, basr%nsh_id(jzp)
                  do jpr = 1, basr%cgto(jsh,jat)%nprim
                     call check(error, numdr(ic, iat, jat, jsh, jpr), basr%cgto(jsh,jat)%dcoeffdr(ic,iat,jpr), thr=thr)
                     if (allocated(error)) then 
                        call test_failed(error, "Coefficient derivative does not match")
                        exit num
                     end if
                  end do
               end do
            end do 
         end do
      end do num

   end subroutine test_qvszp_coefficient_grad


   subroutine test_qvszp_overlap_grad(mol, ksig, kpi, kdel, error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Scaling factors for the diatomic frame overlap
      real(wp), allocatable, intent(in) :: ksig(:,:), kpi(:,:), kdel(:,:)
      !> Molecular structure data
      type(structure_type), intent(inout) :: mol

      type(qvszp_cgto_type), allocatable :: cgto(:, :)
      type(qvszp_basis_type) :: bas
      type(adjacency_list) :: list
      type(context_type) :: ctx
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn
            
      integer, allocatable :: nsh_id(:)
      real(wp), parameter :: accuracy = 1e-8_wp, cn_cutoff = 30.0_wp, step = 1.0e-6_wp
      real(wp) :: vec(3), r2

      real(wp), allocatable :: lattr(:, :), cn(:), dcndr(:, :, :), dcndL(:, :, :)
      real(wp), allocatable :: ql(:), qr(:)
      real(wp), allocatable :: overlapr(:,:), overlapl(:,:) 
      real(wp), allocatable :: overlap_diatr(:,:), overlap_diatl(:,:) 
      real(wp), allocatable :: doverlap(:, :, :), doverlap_diat(:, :, :)
      real(wp), allocatable :: numoverlapdr(:, :, :), numdiatdr(:, :, :), dqatnumdr(:, :, :)

      integer :: ic, iat, is, ish, izp, iao, ii, jat, js, jsh, jzp, jao, jj

      real(wp) :: cutoff

      allocate(ql(mol%nat), qr(mol%nat), dqatnumdr(3, mol%nat, mol%nat))
      allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), source=0.0_wp)

      ! Setup q-vSZP CGTOs and basis 
      call add_qvszp_cgtos(cgto, mol, nsh_id, .false.)
      call new_basis(bas, mol, nsh_id, cgto, 1.0_wp)

      allocate(overlapr(bas%nao, bas%nao), overlapl(bas%nao, bas%nao), &
      & overlap_diatr(bas%nao, bas%nao), overlap_diatl(bas%nao, bas%nao), &
      & doverlap(3, bas%nao, bas%nao), doverlap_diat(3, bas%nao, bas%nao), &
      & numoverlapdr(3, bas%nao, bas%nao), numdiatdr(3, bas%nao, bas%nao), source=0.0_wp)
      
      ! Setup the CEH 
      call new_ceh_calculator(calc, mol, error)
      if (allocated(error)) return
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .false.)
      ctx%verbosity = 0

      ! Obtain CEH charges for charge adaptation
      call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .false.)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
      qr = wfn%qat(:,1)

      ! Obtain CN for charge adaptation
      call bas%ncoord%get_cn(mol, cn)

      ! Scale the basis set 
      call bas%scale_basis(mol, wfn%qat(:,1), cn, .true.)


      do ic = 1, 3
         do iat = 1, mol%nat
            izp = mol%id(iat)
            is = bas%ish_at(iat)

            ! Right hand
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            
            ! Obtain CEH charges for charge adaptation
            call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .false.)
            if (ctx%failed()) then
               call ctx%get_error(error)
               return
            end if
            qr = wfn%qat(:,1)

            ! Obtain CN for charge adaptation
            call bas%ncoord%get_cn(mol, cn)

            ! Scale the basis set 
            call bas%scale_basis(mol, wfn%qat(:,1), cn, .true.)

            ! Cutoffs and adjacency list
            cutoff = get_cutoff(bas)
            call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
            call new_adjacency_list(list, mol, lattr, cutoff)

            ! Get overlap matrix
            call get_overlap(mol, lattr, list, bas, &
            & ksig, kpi, kdel, overlapr, overlap_diatr)

            ! Left hand
            mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step

            ! Obtain CEH charges for charge adaptation
            call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .false.)
            if (ctx%failed()) then
               call ctx%get_error(error)
               return
            end if
            ql = wfn%qat(:,1)

            ! Obtain CN for charge adaptation
            call bas%ncoord%get_cn(mol, cn)

            ! Scale the basis set 
            call bas%scale_basis(mol, wfn%qat(:,1), cn, .true.)

            ! Cutoffs and adjacency list
            cutoff = get_cutoff(bas)
            call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
            call new_adjacency_list(list, mol, lattr, cutoff)

            ! Get overlap matrix
            call get_overlap(mol, lattr, list, bas, &
            & ksig, kpi, kdel, overlapl, overlap_diatl)

            ! Geometry reset 
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

            ! Numerical gradient of the overlap matrices
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
                           numoverlapdr(ic, jj+jao, ii+iao) = & 
                              & + 0.5_wp*(overlapr(jj+jao, ii+iao) - overlapl(jj+jao, ii+iao))/step
                           numdiatdr(ic, jj+jao, ii+iao) = & 
                              & + 0.5_wp*(overlap_diatr(jj+jao, ii+iao) - overlap_diatl(jj+jao, ii+iao))/step

                           ! Lower triangular matrix
                           if(ii+iao /= jj+jao) then
                              numoverlapdr(ic, ii+iao, jj+jao) = &
                                 & - 0.5_wp*(overlapr(ii+iao, jj+jao) - overlapl(ii+iao, jj+jao))/step
                              numdiatdr(ic, ii+iao, jj+jao) = &
                                 & - 0.5_wp*(overlap_diatr(ii+iao, jj+jao) - overlap_diatl(ii+iao, jj+jao))/step
                           end if
                        end do
                     end do 
                  end do
               end do
            end do
            ! CEH gradient (intermediate solution)
            dqatnumdr(ic, iat, :) = 0.5_wp*(qr - ql)/step

         end do
      end do

      ! Obtain CEH charges and derivatives for charge adaptation
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt, .true.)
      call ceh_singlepoint(ctx, calc, mol, wfn, accuracy, .true.)
      if (ctx%failed()) then
         call ctx%get_error(error)
         return
      end if
      wfn%dqatdr(:, :, :, 1) = dqatnumdr
      wfn%dqatdL(:, :, :, 1) = 0.0_wp
       
      ! Obtain CN for charge adaptation
      call bas%ncoord%get_cn(mol, cn, dcndr, dcndL)

      ! Scale the basis set 
      call bas%scale_basis(mol, wfn%qat(:,1), cn, .true., dqatdr=wfn%dqatdr(:, :, :, 1), &
      & dqatdL=wfn%dqatdL(:, :, :, 1), dcndr=dcndr(:, :, :), dcndL=dcndL(:, :, :))

      ! Cutoffs and adjacency list
      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      ! Get overlap matrix gradient
      call get_overlap_gradient(mol, lattr, list, bas, ksig, kpi, kdel, &
      & overlapl, overlap_diatl, doverlap, doverlap_diat)
      
      ! Check results
      num: do ic = 1, 3
         write(*,*) "CURRENT DIMENSION", ic
         call write_2d_matrix(numoverlapdr(ic, :, :), "numerical overlap")
         call write_2d_matrix(doverlap(ic, :, :), "analytical overlap")
         call write_2d_matrix(doverlap(ic, :, :)-numoverlapdr(ic, :, :), "diff")
         ! Check overlap matrix derivative
         do ii = 1, size(numoverlapdr,2)
            do jj = 1, size(numoverlapdr,3)
               call check(error, numoverlapdr(ic, ii, jj), doverlap(ic, ii, jj), thr=thr2)
               if (allocated(error)) then 
                  call test_failed(error, "Overlap derivative does not match")
                  exit num
               end if
            end do           
         end do

         ! call write_2d_matrix(numdiatdr(ic, :, :), "numerical diat")
         ! call write_2d_matrix(doverlap_diat(ic, :, :), "analytical diat")
         ! call write_2d_matrix(doverlap_diat(ic, :, :)-numdiatdr(ic, :, :), "diff")

         ! ! Check diatomic frame overlap matrix derivative
         ! do ii = 1, size(numdiatdr,2)
         !    do jj = 1, size(numdiatdr,3)
         !       call check(error, numdiatdr(ic, ii, jj), doverlap_diat(ic, ii, jj), thr=thr2)
         !       if (allocated(error)) then 
         !          call test_failed(error, "Overlap derivative does not match")
         !          !exit num
         !       end if
         !    end do           
         ! end do
      end do num
   
   end subroutine test_qvszp_overlap_grad

   subroutine test_qvszp_coefficient_grad_lih(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol

      real(wp), allocatable :: ksig(:,:), kpi(:,:), kdel(:,:) 

      call get_structure(mol, "MB16-43", "LiH")
      call test_qvszp_coefficient_grad(mol, error)
   
   end subroutine test_qvszp_coefficient_grad_lih
   
   subroutine test_qvszp_coefficient_grad_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol

      real(wp), allocatable :: ksig(:,:), kpi(:,:), kdel(:,:) 

      call get_structure(mol, "MB16-43", "SiH4")
      call test_qvszp_coefficient_grad(mol, error)
   
   end subroutine test_qvszp_coefficient_grad_sih4
   
   subroutine test_qvszp_coefficient_grad_accl3(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol

      real(wp), allocatable :: ksig(:,:), kpi(:,:), kdel(:,:) 

      call get_structure(mol, "f-block", "AcCl3")
      call test_qvszp_coefficient_grad(mol, error)
   
   end subroutine test_qvszp_coefficient_grad_accl3

   subroutine test_qvszp_overlap_diat_grad_lih(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      type(structure_type) :: mol

      real(wp), allocatable :: ksig(:,:), kpi(:,:), kdel(:,:) 

      call get_structure(mol, "MB16-43", "LiH")

      ! Set diatomic frame scaling constants
      allocate(ksig(mol%nid, mol%nid), kpi(mol%nid, mol%nid), kdel(mol%nid, mol%nid))
      ksig = 0.5_wp
      kpi = 1.2_wp
      kdel = 1.5_wp
      
      call test_qvszp_overlap_grad(mol, ksig, kpi, kdel, error)
   
   end subroutine test_qvszp_overlap_diat_grad_lih
   
   subroutine test_qvszp_norm_cgto_mb01(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "01")
      call test_qvszp_norm_cgto(error, mol)

   end subroutine test_qvszp_norm_cgto_mb01

   subroutine test_qvszp_norm_cgto_accl3(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol

      call get_structure(mol, "f-block", "AcCl3")
      call test_qvszp_norm_cgto(error, mol)

   end subroutine test_qvszp_norm_cgto_accl3

   subroutine test_qvszp_overlap_lih(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 8
      real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 5.028293313832283E-01_wp, 0.000000000000000E+00_wp,&
      &-3.551095997005579E-01_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 5.675379345204321E-01_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 5.418471011098696E-01_wp,&
      & 0.000000000000000E+00_wp,-1.851146945604232E-05_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 5.675379345204321E-01_wp, 5.028293313832282E-01_wp,&
      & 0.000000000000000E+00_wp, 5.418471011098696E-01_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 5.675379345204322E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 9.999999999999999E-01_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      &-3.551095997005579E-01_wp, 0.000000000000000E+00_wp,-1.851146945604232E-05_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 9.999999999999999E-01_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 5.675379345204322E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 9.999999999999999E-01_wp],shape(overlap))
   
      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "LiH")
      call test_qvszp_overlap(error, mol, overlap)

   end subroutine test_qvszp_overlap_lih

   subroutine test_qvszp_overlap_sih4(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 25
      real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 4.260461489313604E-01_wp,-3.027344163816017E-01_wp, 3.027344163816017E-01_wp,&
      &-3.027344163816017E-01_wp, 4.260461489313604E-01_wp, 3.027344163816017E-01_wp,&
      & 3.027344163816017E-01_wp, 3.027344163816017E-01_wp, 4.260461489313604E-01_wp,&
      & 3.027344163816017E-01_wp,-3.027344163816017E-01_wp,-3.027344163816017E-01_wp,&
      & 4.260461489313604E-01_wp,-3.027344163816017E-01_wp,-3.027344163816017E-01_wp,&
      & 3.027344163816017E-01_wp, 0.000000000000000E+00_wp, 9.999999999999999E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 2.984622561250789E-01_wp, 2.663239724594116E-01_wp,&
      & 2.028550846193388E-01_wp,-2.028550846193388E-01_wp,-2.984622561250789E-01_wp,&
      & 2.663239724594116E-01_wp,-2.028550846193388E-01_wp,-2.028550846193388E-01_wp,&
      &-2.984622561250789E-01_wp, 2.663239724594116E-01_wp, 2.028550846193388E-01_wp,&
      & 2.028550846193388E-01_wp, 2.984622561250789E-01_wp, 2.663239724594116E-01_wp,&
      &-2.028550846193388E-01_wp, 2.028550846193388E-01_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 9.999999999999999E-01_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-2.984622561250789E-01_wp,&
      & 2.028550846193388E-01_wp, 2.663239724594116E-01_wp, 2.028550846193388E-01_wp,&
      &-2.984622561250789E-01_wp,-2.028550846193388E-01_wp, 2.663239724594116E-01_wp,&
      &-2.028550846193388E-01_wp, 2.984622561250789E-01_wp, 2.028550846193388E-01_wp,&
      & 2.663239724594116E-01_wp,-2.028550846193388E-01_wp, 2.984622561250789E-01_wp,&
      &-2.028550846193388E-01_wp, 2.663239724594116E-01_wp, 2.028550846193388E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 9.999999999999999E-01_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 2.984622561250789E-01_wp,-2.028550846193388E-01_wp, 2.028550846193388E-01_wp,&
      & 2.663239724594116E-01_wp,-2.984622561250789E-01_wp,-2.028550846193388E-01_wp,&
      &-2.028550846193388E-01_wp, 2.663239724594116E-01_wp, 2.984622561250789E-01_wp,&
      & 2.028550846193388E-01_wp,-2.028550846193388E-01_wp, 2.663239724594116E-01_wp,&
      &-2.984622561250789E-01_wp, 2.028550846193388E-01_wp, 2.028550846193388E-01_wp,&
      & 2.663239724594116E-01_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 9.999999999999998E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 2.059001252264858E-01_wp, 5.472570098843951E-02_wp,&
      & 1.577117175416436E-01_wp, 5.472570098843951E-02_wp, 2.059001252264858E-01_wp,&
      &-5.472570098843951E-02_wp, 1.577117175416436E-01_wp,-5.472570098843951E-02_wp,&
      &-2.059001252264858E-01_wp, 5.472570098843951E-02_wp, 1.577117175416436E-01_wp,&
      &-5.472570098843951E-02_wp,-2.059001252264858E-01_wp,-5.472570098843951E-02_wp,&
      & 1.577117175416436E-01_wp, 5.472570098843951E-02_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 9.999999999999998E-01_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-2.059001252264858E-01_wp,&
      &-5.472570098843951E-02_wp, 5.472570098843951E-02_wp, 1.577117175416436E-01_wp,&
      & 2.059001252264858E-01_wp,-5.472570098843951E-02_wp,-5.472570098843951E-02_wp,&
      & 1.577117175416436E-01_wp,-2.059001252264858E-01_wp, 5.472570098843951E-02_wp,&
      &-5.472570098843951E-02_wp, 1.577117175416436E-01_wp, 2.059001252264858E-01_wp,&
      & 5.472570098843951E-02_wp, 5.472570098843951E-02_wp, 1.577117175416436E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 9.999999999999998E-01_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp,-1.226508007742927E-01_wp,-2.453016015485854E-01_wp,&
      &-1.226508007742927E-01_wp, 0.000000000000000E+00_wp, 1.226508007742927E-01_wp,&
      &-2.453016015485854E-01_wp, 1.226508007742927E-01_wp, 0.000000000000000E+00_wp,&
      & 1.226508007742927E-01_wp, 2.453016015485854E-01_wp,-1.226508007742927E-01_wp,&
      & 0.000000000000000E+00_wp,-1.226508007742927E-01_wp, 2.453016015485854E-01_wp,&
      & 1.226508007742927E-01_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 9.999999999999998E-01_wp,&
      & 0.000000000000000E+00_wp,-2.059001252264858E-01_wp, 1.577117175416436E-01_wp,&
      & 5.472570098843951E-02_wp,-5.472570098843951E-02_wp, 2.059001252264858E-01_wp,&
      & 1.577117175416436E-01_wp,-5.472570098843951E-02_wp,-5.472570098843951E-02_wp,&
      & 2.059001252264858E-01_wp, 1.577117175416436E-01_wp, 5.472570098843951E-02_wp,&
      & 5.472570098843951E-02_wp,-2.059001252264858E-01_wp, 1.577117175416436E-01_wp,&
      &-5.472570098843951E-02_wp, 5.472570098843951E-02_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 9.999999999999998E-01_wp, 0.000000000000000E+00_wp,&
      &-2.124374185300831E-01_wp, 0.000000000000000E+00_wp, 2.124374185300831E-01_wp,&
      & 0.000000000000000E+00_wp, 2.124374185300831E-01_wp, 0.000000000000000E+00_wp,&
      &-2.124374185300831E-01_wp, 0.000000000000000E+00_wp, 2.124374185300831E-01_wp,&
      & 0.000000000000000E+00_wp, 2.124374185300831E-01_wp, 0.000000000000000E+00_wp,&
      &-2.124374185300831E-01_wp, 0.000000000000000E+00_wp,-2.124374185300831E-01_wp,&
      & 4.260461489313605E-01_wp, 2.984622561250790E-01_wp,-2.984622561250790E-01_wp,&
      & 2.984622561250790E-01_wp, 2.059001252264858E-01_wp,-2.059001252264858E-01_wp,&
      & 0.000000000000000E+00_wp,-2.059001252264858E-01_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 8.718465237414751E-02_wp, 1.712464585372617E-01_wp,&
      & 0.000000000000000E+00_wp, 1.712464585372617E-01_wp, 8.718465237414751E-02_wp,&
      & 1.712464585372617E-01_wp,-1.712464585372617E-01_wp, 0.000000000000000E+00_wp,&
      & 8.718465237414751E-02_wp, 0.000000000000000E+00_wp,-1.712464585372617E-01_wp,&
      & 1.712464585372617E-01_wp,-3.027344163816017E-01_wp, 2.663239724594116E-01_wp,&
      & 2.028550846193388E-01_wp,-2.028550846193388E-01_wp, 5.472570098843951E-02_wp,&
      &-5.472570098843951E-02_wp,-1.226508007742927E-01_wp, 1.577117175416436E-01_wp,&
      &-2.124374185300831E-01_wp, 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-1.712464585372618E-01_wp,&
      &-4.548596649725484E-02_wp, 0.000000000000000E+00_wp,-2.305326168717275E-01_wp,&
      &-1.712464585372618E-01_wp,-4.548596649725484E-02_wp, 2.305326168717275E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.850466503744726E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 3.027344163816017E-01_wp,&
      & 2.028550846193388E-01_wp, 2.663239724594116E-01_wp, 2.028550846193388E-01_wp,&
      & 1.577117175416436E-01_wp, 5.472570098843951E-02_wp,-2.453016015485854E-01_wp,&
      & 5.472570098843951E-02_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.850466503744726E-01_wp,&
      & 0.000000000000000E+00_wp, 1.712464585372618E-01_wp, 2.305326168717275E-01_wp,&
      &-4.548596649725484E-02_wp, 0.000000000000000E+00_wp, 1.712464585372618E-01_wp,&
      & 0.000000000000000E+00_wp,-4.548596649725484E-02_wp, 2.305326168717275E-01_wp,&
      &-3.027344163816017E-01_wp,-2.028550846193388E-01_wp, 2.028550846193388E-01_wp,&
      & 2.663239724594116E-01_wp, 5.472570098843951E-02_wp, 1.577117175416436E-01_wp,&
      &-1.226508007742927E-01_wp,-5.472570098843951E-02_wp, 2.124374185300831E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp,-1.712464585372618E-01_wp,-2.305326168717275E-01_wp,&
      & 0.000000000000000E+00_wp,-4.548596649725484E-02_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.850466503744726E-01_wp,&
      &-1.712464585372618E-01_wp, 0.000000000000000E+00_wp, 2.305326168717275E-01_wp,&
      &-4.548596649725484E-02_wp, 4.260461489313605E-01_wp,-2.984622561250790E-01_wp,&
      &-2.984622561250790E-01_wp,-2.984622561250790E-01_wp, 2.059001252264858E-01_wp,&
      & 2.059001252264858E-01_wp, 0.000000000000000E+00_wp, 2.059001252264858E-01_wp,&
      & 0.000000000000000E+00_wp, 8.718465237414751E-02_wp,-1.712464585372617E-01_wp,&
      & 0.000000000000000E+00_wp,-1.712464585372617E-01_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 8.718465237414751E-02_wp, 0.000000000000000E+00_wp,-1.712464585372617E-01_wp,&
      &-1.712464585372617E-01_wp, 8.718465237414751E-02_wp,-1.712464585372617E-01_wp,&
      &-1.712464585372617E-01_wp, 0.000000000000000E+00_wp, 3.027344163816017E-01_wp,&
      & 2.663239724594116E-01_wp,-2.028550846193388E-01_wp,-2.028550846193388E-01_wp,&
      &-5.472570098843951E-02_wp,-5.472570098843951E-02_wp, 1.226508007742927E-01_wp,&
      & 1.577117175416436E-01_wp, 2.124374185300831E-01_wp, 1.712464585372618E-01_wp,&
      &-4.548596649725484E-02_wp, 0.000000000000000E+00_wp,-2.305326168717275E-01_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.850466503744726E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.712464585372618E-01_wp,&
      &-4.548596649725484E-02_wp,-2.305326168717275E-01_wp, 0.000000000000000E+00_wp,&
      & 3.027344163816017E-01_wp,-2.028550846193388E-01_wp, 2.663239724594116E-01_wp,&
      &-2.028550846193388E-01_wp, 1.577117175416436E-01_wp,-5.472570098843951E-02_wp,&
      &-2.453016015485854E-01_wp,-5.472570098843951E-02_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.850466503744726E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.712464585372618E-01_wp,&
      & 0.000000000000000E+00_wp,-4.548596649725484E-02_wp,-2.305326168717275E-01_wp,&
      & 1.712464585372618E-01_wp,-2.305326168717275E-01_wp,-4.548596649725484E-02_wp,&
      & 0.000000000000000E+00_wp, 3.027344163816017E-01_wp,-2.028550846193388E-01_wp,&
      &-2.028550846193388E-01_wp, 2.663239724594116E-01_wp,-5.472570098843951E-02_wp,&
      & 1.577117175416436E-01_wp, 1.226508007742927E-01_wp,-5.472570098843951E-02_wp,&
      &-2.124374185300831E-01_wp, 1.712464585372618E-01_wp,-2.305326168717275E-01_wp,&
      & 0.000000000000000E+00_wp,-4.548596649725484E-02_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,&
      & 1.712464585372618E-01_wp, 0.000000000000000E+00_wp,-2.305326168717275E-01_wp,&
      &-4.548596649725484E-02_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.850466503744726E-01_wp, 4.260461489313605E-01_wp,&
      &-2.984622561250790E-01_wp, 2.984622561250790E-01_wp, 2.984622561250790E-01_wp,&
      &-2.059001252264858E-01_wp,-2.059001252264858E-01_wp, 0.000000000000000E+00_wp,&
      & 2.059001252264858E-01_wp, 0.000000000000000E+00_wp, 8.718465237414751E-02_wp,&
      &-1.712464585372617E-01_wp, 1.712464585372617E-01_wp, 0.000000000000000E+00_wp,&
      & 8.718465237414751E-02_wp, 0.000000000000000E+00_wp, 1.712464585372617E-01_wp,&
      & 1.712464585372617E-01_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 8.718465237414751E-02_wp,&
      &-1.712464585372617E-01_wp, 0.000000000000000E+00_wp, 1.712464585372617E-01_wp,&
      & 3.027344163816017E-01_wp, 2.663239724594116E-01_wp, 2.028550846193388E-01_wp,&
      & 2.028550846193388E-01_wp, 5.472570098843951E-02_wp, 5.472570098843951E-02_wp,&
      & 1.226508007742927E-01_wp, 1.577117175416436E-01_wp, 2.124374185300831E-01_wp,&
      & 1.712464585372618E-01_wp,-4.548596649725484E-02_wp, 2.305326168717275E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.850466503744726E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.712464585372618E-01_wp,-4.548596649725484E-02_wp, 0.000000000000000E+00_wp,&
      & 2.305326168717275E-01_wp,-3.027344163816017E-01_wp, 2.028550846193388E-01_wp,&
      & 2.663239724594116E-01_wp,-2.028550846193388E-01_wp, 1.577117175416436E-01_wp,&
      &-5.472570098843951E-02_wp, 2.453016015485854E-01_wp, 5.472570098843951E-02_wp,&
      & 0.000000000000000E+00_wp,-1.712464585372618E-01_wp, 2.305326168717275E-01_wp,&
      &-4.548596649725484E-02_wp, 0.000000000000000E+00_wp,-1.712464585372618E-01_wp,&
      & 0.000000000000000E+00_wp,-4.548596649725484E-02_wp,-2.305326168717275E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.850466503744726E-01_wp, 0.000000000000000E+00_wp,-3.027344163816017E-01_wp,&
      & 2.028550846193388E-01_wp,-2.028550846193388E-01_wp, 2.663239724594116E-01_wp,&
      &-5.472570098843951E-02_wp, 1.577117175416436E-01_wp,-1.226508007742927E-01_wp,&
      & 5.472570098843951E-02_wp, 2.124374185300831E-01_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.850466503744726E-01_wp,&
      &-1.712464585372618E-01_wp, 0.000000000000000E+00_wp,-2.305326168717275E-01_wp,&
      &-4.548596649725484E-02_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,-1.712464585372618E-01_wp,&
      & 2.305326168717275E-01_wp, 0.000000000000000E+00_wp,-4.548596649725484E-02_wp,&
      & 4.260461489313605E-01_wp, 2.984622561250790E-01_wp, 2.984622561250790E-01_wp,&
      &-2.984622561250790E-01_wp,-2.059001252264858E-01_wp, 2.059001252264858E-01_wp,&
      & 0.000000000000000E+00_wp,-2.059001252264858E-01_wp, 0.000000000000000E+00_wp,&
      & 8.718465237414751E-02_wp, 0.000000000000000E+00_wp, 1.712464585372617E-01_wp,&
      &-1.712464585372617E-01_wp, 8.718465237414751E-02_wp, 1.712464585372617E-01_wp,&
      & 1.712464585372617E-01_wp, 0.000000000000000E+00_wp, 8.718465237414751E-02_wp,&
      & 1.712464585372617E-01_wp, 0.000000000000000E+00_wp,-1.712464585372617E-01_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp,-3.027344163816017E-01_wp, 2.663239724594116E-01_wp,&
      &-2.028550846193388E-01_wp, 2.028550846193388E-01_wp,-5.472570098843951E-02_wp,&
      & 5.472570098843951E-02_wp,-1.226508007742927E-01_wp, 1.577117175416436E-01_wp,&
      &-2.124374185300831E-01_wp, 0.000000000000000E+00_wp, 1.850466503744726E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-1.712464585372618E-01_wp,&
      &-4.548596649725484E-02_wp,-2.305326168717275E-01_wp, 0.000000000000000E+00_wp,&
      &-1.712464585372618E-01_wp,-4.548596649725484E-02_wp, 0.000000000000000E+00_wp,&
      & 2.305326168717275E-01_wp, 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-3.027344163816017E-01_wp,&
      &-2.028550846193388E-01_wp, 2.663239724594116E-01_wp, 2.028550846193388E-01_wp,&
      & 1.577117175416436E-01_wp, 5.472570098843951E-02_wp, 2.453016015485854E-01_wp,&
      &-5.472570098843951E-02_wp, 0.000000000000000E+00_wp,-1.712464585372618E-01_wp,&
      & 0.000000000000000E+00_wp,-4.548596649725484E-02_wp, 2.305326168717275E-01_wp,&
      &-1.712464585372618E-01_wp,-2.305326168717275E-01_wp,-4.548596649725484E-02_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.850466503744726E-01_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 3.027344163816017E-01_wp, 2.028550846193388E-01_wp, 2.028550846193388E-01_wp,&
      & 2.663239724594116E-01_wp, 5.472570098843951E-02_wp, 1.577117175416436E-01_wp,&
      & 1.226508007742927E-01_wp, 5.472570098843951E-02_wp,-2.124374185300831E-01_wp,&
      & 1.712464585372618E-01_wp, 0.000000000000000E+00_wp, 2.305326168717275E-01_wp,&
      &-4.548596649725484E-02_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.850466503744726E-01_wp, 1.712464585372618E-01_wp,&
      & 2.305326168717275E-01_wp, 0.000000000000000E+00_wp,-4.548596649725484E-02_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp],shape(overlap))

      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "SiH4")
      call test_qvszp_overlap(error, mol, overlap)

   end subroutine test_qvszp_overlap_sih4


   subroutine test_qvszp_overlap_accl3(error)

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: nao = 43
      real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 9.999999999999996E-01_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 2.546244146896991E-01_wp, 8.904521860262070E-02_wp,&
      & 1.194320942043015E-01_wp, 6.183091478742130E-02_wp, 6.047256261536093E-03_wp,&
      & 1.168083121507140E-02_wp, 5.319101370119733E-03_wp, 8.110895686925459E-03_wp,&
      &-2.254913440674616E-03_wp, 2.545368264525664E-01_wp,-1.479529014292982E-01_wp,&
      &-1.507479321872774E-02_wp, 6.248399057037402E-02_wp,-1.017180391830452E-02_wp,&
      & 2.454034054646125E-03_wp,-8.048565621540937E-03_wp,-1.036396307531422E-03_wp,&
      &-9.894771926742629E-03_wp, 2.546290040058253E-01_wp, 8.702694713955239E-03_wp,&
      &-1.537100523094241E-02_wp,-1.603241986148280E-01_wp,-1.532338905733210E-03_wp,&
      &-1.469122536653884E-04_wp,-8.023297592600592E-03_wp, 2.706470824241631E-03_wp,&
      & 1.407306246330118E-02_wp, 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp,-1.387778780781446E-17_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-1.973102049062196E-01_wp,&
      & 1.121432576074864E-01_wp,-6.045790687668998E-02_wp,-3.129952391125229E-02_wp,&
      & 9.226408301780580E-03_wp, 1.782165554648692E-02_wp,-4.754718969258193E-03_wp,&
      & 2.286344247452108E-03_wp,-1.146804083903603E-02_wp, 3.278091162485747E-01_wp,&
      & 3.270230911069708E-02_wp,-1.268353560535391E-02_wp, 5.257239072303660E-02_wp,&
      & 1.233135000684818E-02_wp,-2.975042883210660E-03_wp, 1.413549772308417E-02_wp,&
      & 4.817287462261999E-04_wp, 2.260288857247389E-02_wp,-1.928389975975782E-02_wp,&
      & 1.567901565919992E-01_wp, 7.604487167628129E-04_wp, 7.931708386725814E-03_wp,&
      &-1.954549296314300E-02_wp,-1.873914712647955E-03_wp,-8.323523586443188E-04_wp,&
      & 7.459203931459591E-05_wp,-6.708106897412503E-04_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      &-2.646427438736851E-01_wp,-6.045790687668996E-02_wp, 7.612967691506115E-02_wp,&
      &-4.198055490211861E-02_wp, 2.286344247452107E-03_wp, 1.524869736605941E-02_wp,&
      & 1.878768563866048E-02_wp, 1.058833839992669E-02_wp,-8.525367787670529E-04_wp,&
      & 3.340018745778095E-02_wp,-1.268353560535391E-02_wp, 1.558936851198954E-01_wp,&
      & 5.356555441007389E-03_wp, 4.817287462261997E-04_wp,-1.811990762199503E-02_wp,&
      &-1.736985804602402E-03_wp, 7.652463223438904E-03_wp, 4.686087259199180E-04_wp,&
      & 3.405990142396841E-02_wp, 7.604487167628129E-04_wp, 1.558775731802586E-01_wp,&
      &-1.400926208604901E-02_wp, 7.459203931459594E-05_wp, 1.065824624332955E-03_wp,&
      &-1.768572441644346E-03_wp,-1.963500781960325E-02_wp,-6.850563048498959E-04_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp,-1.387778780781446E-17_wp, 0.000000000000000E+00_wp,&
      & 1.387778780781446E-17_wp,-1.370075862319917E-01_wp,-3.129952391125229E-02_wp,&
      &-4.198055490211861E-02_wp, 1.354853202229285E-01_wp, 1.201607135749007E-02_wp,&
      & 2.286344247452108E-03_wp,-3.301565519630095E-03_wp, 1.611658198895281E-02_wp,&
      & 7.080409322663008E-03_wp,-1.384415008471236E-01_wp, 5.257239072303659E-02_wp,&
      & 5.356555441007388E-03_wp, 1.349834423346564E-01_wp,-2.000041929617723E-02_wp,&
      & 4.817287462262003E-04_wp,-5.969753197836450E-03_wp,-2.037825431370820E-03_wp,&
      & 5.661029063890270E-03_wp, 3.552549959260523E-01_wp, 7.931708386725814E-03_wp,&
      &-1.400926208604901E-02_wp, 1.109989558061869E-02_wp, 1.836690509563453E-03_wp,&
      & 7.459203931459594E-05_wp, 1.533389705729007E-02_wp,-3.244027322347748E-03_wp,&
      &-2.664860363771760E-02_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 9.999999999999998E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 9.336173778708011E-02_wp,&
      & 1.585601315035428E-02_wp, 1.131744594455876E-01_wp,-4.009246557822630E-02_wp,&
      & 2.017654828951673E-03_wp, 3.228777988126517E-02_wp, 6.236274278382165E-02_wp,&
      &-4.582692336889597E-03_wp,-1.322267895469086E-02_wp,-1.566460435169708E-01_wp,&
      & 1.660461878936972E-01_wp, 2.396835040534949E-02_wp, 6.449243520200129E-02_wp,&
      &-2.309479251972241E-03_wp,-1.948093548852921E-02_wp, 2.674150463780350E-02_wp,&
      &-7.477946722906489E-04_wp, 9.714846767178610E-02_wp,-2.366042499400333E-02_wp,&
      & 1.755955375970396E-01_wp, 3.691329429922923E-03_wp, 2.885657677088460E-02_wp,&
      &-9.995084309531216E-02_wp,-1.186547464654780E-02_wp, 4.024635858893757E-03_wp,&
      &-3.366551588210808E-03_wp,-2.091604547977351E-02_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 9.999999999999998E-01_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.803367765267521E-01_wp, 3.062734657555167E-02_wp, 1.199231145415588E-01_wp,&
      & 1.131744594455876E-01_wp, 3.228777988126517E-02_wp, 4.766883667501617E-02_wp,&
      & 3.014243757627399E-02_wp, 6.487439978089043E-02_wp, 2.660368501404787E-02_wp,&
      & 3.779218793477304E-02_wp,-4.006005257356781E-02_wp, 1.580570462019545E-01_wp,&
      & 2.396835040534949E-02_wp,-1.948093548852921E-02_wp,-7.835668817848246E-02_wp,&
      &-2.537132714397605E-02_wp, 4.282141048512239E-02_wp,-1.251459174206803E-02_wp,&
      &-2.268431836811374E-03_wp, 1.683513749174360E-02_wp,-9.291192814470115E-03_wp,&
      & 3.691329429922922E-03_wp,-1.186547464654780E-02_wp, 2.267202257691399E-02_wp,&
      & 1.521968071603357E-03_wp, 6.455916960609892E-03_wp,-2.661248062415271E-03_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 9.999999999999997E-01_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 8.211997737530544E-02_wp, 1.311945676470954E-01_wp,&
      &-5.328890077053486E-02_wp, 9.109843582911245E-02_wp, 6.236274278382165E-02_wp,&
      & 3.014243757627405E-02_wp,-6.084476365669887E-02_wp, 2.093020286222212E-02_wp,&
      &-2.325394870315192E-02_wp,-1.239481188122672E-01_wp, 9.154332634880821E-02_wp,&
      & 3.824122009652011E-02_wp,-3.866090009118883E-02_wp, 2.674150463780352E-02_wp,&
      &-2.537132714397604E-02_wp, 4.308469549368481E-02_wp, 1.071490826274628E-02_wp,&
      & 2.601319210379564E-02_wp,-1.238855387564933E-01_wp,-5.374298642707931E-03_wp,&
      & 3.899864345132841E-02_wp, 9.900727893249393E-02_wp, 4.024635858893756E-03_wp,&
      & 1.521968071603357E-03_wp, 4.293684846532067E-02_wp,-2.803824785510204E-02_wp,&
      &-3.696241844564563E-02_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 9.999999999999998E-01_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.252216350012496E-01_wp,&
      & 1.131744594455876E-01_wp, 8.327180271578394E-02_wp,-5.377410714426065E-02_wp,&
      &-4.582692336889596E-03_wp, 6.487439978089043E-02_wp, 2.093020286222212E-02_wp,&
      &-7.121644583707668E-04_wp,-5.394287007897097E-02_wp,-1.596052995066541E-02_wp,&
      & 2.396835040534949E-02_wp,-6.675120858771078E-02_wp, 6.571078467879550E-03_wp,&
      &-7.477946722906489E-04_wp, 4.282141048512239E-02_wp, 1.071490826274629E-02_wp,&
      & 4.953626045965858E-03_wp, 1.451154029554059E-02_wp, 4.178987409106262E-02_wp,&
      & 3.691329429922922E-03_wp, 1.711657240793602E-01_wp,-5.096749996079889E-02_wp,&
      &-3.366551588210808E-03_wp, 6.455916960609892E-03_wp,-2.803824785510203E-02_wp,&
      &-9.591078218594690E-02_wp, 2.485883202263236E-02_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 9.999999999999996E-01_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      &-3.481291816916520E-02_wp, 6.722007672201609E-02_wp,-4.220072685990615E-02_wp,&
      &-9.037132978064609E-02_wp,-1.322267895469085E-02_wp, 2.660368501404788E-02_wp,&
      &-2.325394870315194E-02_wp,-5.394287007897097E-02_wp,-2.851259805292793E-02_wp,&
      &-1.523797436792693E-01_wp, 6.499308133004139E-02_wp, 2.331556552072368E-02_wp,&
      &-1.658347532571246E-01_wp, 9.714846767178610E-02_wp,-1.251459174206803E-02_wp,&
      & 2.601319210379562E-02_wp, 1.451154029554059E-02_wp,-7.675290382092213E-03_wp,&
      & 2.172982997449486E-01_wp, 2.883921128113612E-02_wp,-3.390131872627126E-02_wp,&
      &-1.759154509479025E-01_wp,-2.091604547977351E-02_wp,-2.661248062415271E-03_wp,&
      &-3.696241844564566E-02_wp, 2.485883202263235E-02_wp, 8.986553028404413E-02_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,-9.574257339115482E-17_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp,-6.862885944333728E-03_wp,-2.352698810721638E-02_wp,&
      &-1.004434735863699E-02_wp, 3.781143000996426E-02_wp, 2.439371290099738E-03_wp,&
      &-2.356763940357238E-02_wp,-9.573576200081020E-03_wp, 3.006495495986367E-02_wp,&
      & 1.272603818558988E-02_wp,-3.275114472081729E-02_wp,-1.080899641447738E-02_wp,&
      & 6.049946259626841E-03_wp,-9.722842322371636E-02_wp, 7.285918079573005E-02_wp,&
      &-2.226924160852676E-03_wp, 7.860618424361789E-03_wp, 1.195187549207985E-02_wp,&
      &-3.655550738162817E-02_wp,-1.459634274178628E-02_wp, 9.855306067937843E-02_wp,&
      & 2.749421233824322E-03_wp, 1.777688848623613E-02_wp,-7.983685125583712E-02_wp,&
      &-1.078538960739731E-02_wp, 3.485514967844198E-03_wp,-3.634650375220898E-03_wp,&
      &-2.067213199823942E-02_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.000000000000001E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-7.013249410404437E-02_wp,&
      &-2.942554931395643E-02_wp,-6.752544056763578E-02_wp, 1.469535563625872E-02_wp,&
      &-7.929227533829629E-03_wp,-4.796279121035552E-02_wp,-4.567622608156514E-02_wp,&
      &-1.101216010416653E-02_wp, 1.936340914448410E-02_wp,-1.484974722030525E-02_wp,&
      & 2.092008883615450E-02_wp,-5.616858926330325E-02_wp, 2.842922138796015E-03_wp,&
      & 3.228646134519802E-03_wp, 4.814976599275308E-02_wp, 1.460196370138135E-02_wp,&
      &-9.891343763986659E-04_wp, 1.795140989110575E-02_wp,-2.287485624110323E-03_wp,&
      & 1.547578577064782E-02_wp,-8.469260304161204E-03_wp, 3.640903566789722E-03_wp,&
      &-1.510885403830548E-02_wp, 2.497148719603813E-02_wp, 2.247471923823931E-03_wp,&
      & 8.103144802471892E-03_wp,-3.942121353900569E-03_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-1.665334536937735E-16_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      &-6.802976922578063E-02_wp,-4.453680850300276E-02_wp,-1.376138619501135E-02_wp,&
      &-6.265193940464103E-02_wp,-4.096922349440277E-02_wp,-4.262713084311717E-02_wp,&
      & 1.597426822799089E-02_wp,-3.899407999281288E-02_wp,-9.913117013673638E-03_wp,&
      &-6.203146574730595E-02_wp, 4.327699079137380E-02_wp, 2.943686192571209E-02_wp,&
      &-2.886620266640593E-02_wp, 2.963744056419991E-02_wp,-2.696830351838712E-02_wp,&
      & 3.899632680726240E-02_wp, 1.390006426015647E-02_wp, 1.781922459229880E-02_wp,&
      & 3.644875522761238E-03_wp,-2.481186422430601E-02_wp,-1.765906564399253E-03_wp,&
      &-4.346579276711811E-03_wp, 1.989155217895336E-02_wp, 5.949286951402332E-03_wp,&
      &-2.282711389178952E-03_wp, 2.139810523794881E-03_wp, 4.149032133885583E-03_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-1.387778780781446E-17_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,-5.551115123125783E-17_wp,&
      & 0.000000000000000E+00_wp, 1.106450525808549E-02_wp,-4.047119256872936E-02_wp,&
      & 5.763707305879297E-02_wp,-2.810224847927740E-02_wp,-3.845808049292056E-02_wp,&
      & 8.974546488196547E-03_wp, 4.978089179485133E-02_wp, 6.231714940740243E-03_wp,&
      & 1.434032871363807E-02_wp,-1.595284220086065E-02_wp, 1.791315680943756E-02_wp,&
      &-5.959318406889177E-02_wp,-7.565147491895574E-03_wp, 9.952273555086244E-03_wp,&
      & 4.337691668692271E-02_wp, 1.897528876075850E-02_wp,-1.831909227229836E-02_wp,&
      & 9.681220535810810E-03_wp,-1.626917205522779E-02_wp,-1.074097956085215E-03_wp,&
      &-5.945723090416467E-02_wp, 1.978742213800153E-02_wp, 1.531204046169161E-03_wp,&
      &-2.542208638003138E-03_wp, 1.933726609315561E-02_wp, 4.683348962545750E-02_wp,&
      &-1.406263986718217E-02_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.942890293094024E-16_wp,-4.723827882077183E-02_wp,&
      &-6.265193940464103E-02_wp,-9.555584348416628E-03_wp, 2.186761373287551E-03_wp,&
      &-1.084628766033663E-02_wp,-3.899407999281288E-02_wp, 1.109215781707620E-02_wp,&
      &-1.354674420833593E-02_wp, 4.032125164772303E-02_wp, 2.619734715154162E-02_wp,&
      &-2.886620266640594E-02_wp,-1.243187923466680E-02_wp,-1.288305595534449E-02_wp,&
      & 2.838918112505760E-03_wp, 1.390006426015647E-02_wp,-1.646906612148192E-02_wp,&
      & 7.468054342921154E-05_wp,-2.331090475596503E-02_wp,-6.714721892983815E-02_wp,&
      &-4.346579276711811E-03_wp, 3.253217124943841E-02_wp, 5.502645539805805E-02_wp,&
      & 5.022744764462082E-03_wp, 2.139810523794881E-03_wp, 4.205293718418014E-02_wp,&
      &-3.335492577637602E-02_wp,-3.598740316775215E-02_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 2.615114966916810E-02_wp,-3.929888146419137E-02_wp, 2.517902614643967E-02_wp,&
      & 6.691788327147559E-02_wp, 1.936340914448410E-02_wp,-1.401759762715070E-02_wp,&
      & 1.703184579190498E-02_wp, 5.004969665344108E-02_wp, 3.677954743243223E-02_wp,&
      &-1.444531010377486E-02_wp, 1.197635916854938E-02_wp,-5.463882165553090E-02_wp,&
      &-1.706282353260835E-02_wp, 1.795140989110575E-02_wp, 3.296613799925560E-02_wp,&
      & 1.420427503992098E-02_wp,-3.380965573848587E-02_wp, 2.237134098318997E-03_wp,&
      & 2.100836045574692E-02_wp, 3.093775544053226E-03_wp, 7.778202904884629E-02_wp,&
      &-2.555517665681595E-02_wp,-3.942121353900568E-03_wp, 6.112772769082304E-03_wp,&
      &-2.064087301454750E-02_wp,-6.163884273211702E-02_wp, 2.066651179470221E-02_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 8.326672684688674E-17_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 2.687395269863170E-02_wp,-1.368659923965709E-02_wp,&
      & 3.933204165042306E-02_wp, 4.324296291121656E-03_wp, 9.067390599316298E-03_wp,&
      & 1.555359763288201E-03_wp, 3.748857784386517E-02_wp, 1.293692249384755E-02_wp,&
      & 1.011079258997193E-04_wp, 8.394360253427904E-02_wp,-8.003771584250008E-02_wp,&
      &-1.550645904138616E-02_wp,-5.913517130087542E-03_wp,-2.817063608821062E-02_wp,&
      & 1.736265410576016E-02_wp,-2.014734551457818E-02_wp,-3.036490348359322E-03_wp,&
      &-6.931805241115378E-02_wp, 8.892806632054284E-02_wp, 2.038431321445513E-02_wp,&
      &-1.675082026709908E-02_wp,-7.460636655623722E-02_wp,-2.116930580241157E-02_wp,&
      &-2.802806838903087E-03_wp,-2.123546368463537E-02_wp, 1.842507427669205E-02_wp,&
      & 7.067773510318534E-02_wp, 2.546244146896991E-01_wp,-1.973102049062196E-01_wp,&
      &-2.646427438736851E-01_wp,-1.370075862319917E-01_wp, 9.336173778708012E-02_wp,&
      & 1.803367765267521E-01_wp, 8.211997737530551E-02_wp, 1.252216350012496E-01_wp,&
      &-3.481291816916521E-02_wp,-6.862885944333721E-03_wp,-7.013249410404437E-02_wp,&
      &-6.802976922578063E-02_wp, 1.106450525808548E-02_wp,-4.723827882077185E-02_wp,&
      & 2.615114966916811E-02_wp, 2.687395269863169E-02_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.321302646831464E-03_wp,&
      &-4.840842595230028E-03_wp,-2.747224633155594E-03_wp, 1.349004251090167E-05_wp,&
      &-1.821272091314209E-06_wp, 3.708990204362797E-04_wp,-6.714094920657612E-05_wp,&
      &-1.033589391621077E-06_wp,-3.267752137599800E-04_wp, 1.347613073731518E-03_wp,&
      &-1.672739882063343E-03_wp,-2.806610975486558E-03_wp,-4.625269223552698E-03_wp,&
      & 2.131661632649750E-04_wp, 1.293491177497988E-04_wp,-6.710463509118858E-05_wp,&
      & 3.576614032086451E-04_wp, 2.561653298769958E-04_wp, 8.904521860262070E-02_wp,&
      & 1.121432576074864E-01_wp,-6.045790687668996E-02_wp,-3.129952391125229E-02_wp,&
      & 1.585601315035428E-02_wp, 3.062734657555167E-02_wp, 1.311945676470954E-01_wp,&
      & 1.131744594455876E-01_wp, 6.722007672201609E-02_wp,-2.352698810721638E-02_wp,&
      &-2.942554931395643E-02_wp,-4.453680850300275E-02_wp,-4.047119256872936E-02_wp,&
      &-6.265193940464103E-02_wp,-3.929888146419138E-02_wp,-1.368659923965709E-02_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 4.840842595230028E-03_wp,-1.227789019756612E-02_wp,-8.161704368013259E-03_wp,&
      & 4.007744308824233E-05_wp,-9.601871392670786E-06_wp, 1.955405076968439E-03_wp,&
      &-1.606921694867899E-04_wp,-6.093185202772807E-06_wp,-1.519164074610300E-03_wp,&
      & 1.672739882063343E-03_wp, 4.577598973848464E-04_wp,-2.824599855246675E-03_wp,&
      &-4.654914732904501E-03_wp, 3.423087656910362E-05_wp, 2.077127821882286E-05_wp,&
      &-5.223671899238223E-05_wp, 7.137127420796188E-04_wp, 6.526348981436349E-04_wp,&
      & 1.194320942043015E-01_wp,-6.045790687668998E-02_wp, 7.612967691506115E-02_wp,&
      &-4.198055490211861E-02_wp, 1.131744594455876E-01_wp, 1.199231145415587E-01_wp,&
      &-5.328890077053486E-02_wp, 8.327180271578392E-02_wp,-4.220072685990612E-02_wp,&
      &-1.004434735863699E-02_wp,-6.752544056763576E-02_wp,-1.376138619501130E-02_wp,&
      & 5.763707305879297E-02_wp,-9.555584348416642E-03_wp, 2.517902614643968E-02_wp,&
      & 3.933204165042305E-02_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 2.747224633155594E-03_wp,-8.161704368013261E-03_wp,&
      &-2.528121064115085E-03_wp, 2.274433359894836E-05_wp,-6.093185202772808E-06_wp,&
      & 8.336356213639305E-04_wp,-4.914846471812638E-04_wp,-2.323103829461939E-06_wp,&
      &-1.093247904369124E-03_wp, 2.806610975486558E-03_wp,-2.824599855246675E-03_wp,&
      &-2.598041086701424E-03_wp,-7.810260829800153E-03_wp, 7.137127420796187E-04_wp,&
      & 2.916231196675856E-04_wp,-4.987385488180778E-04_wp, 8.063629346907818E-04_wp,&
      & 8.576804930573188E-04_wp, 6.183091478742130E-02_wp,-3.129952391125229E-02_wp,&
      &-4.198055490211861E-02_wp, 1.354853202229285E-01_wp,-4.009246557822630E-02_wp,&
      & 1.131744594455876E-01_wp, 9.109843582911245E-02_wp,-5.377410714426065E-02_wp,&
      &-9.037132978064609E-02_wp, 3.781143000996425E-02_wp, 1.469535563625872E-02_wp,&
      &-6.265193940464105E-02_wp,-2.810224847927740E-02_wp, 2.186761373287554E-03_wp,&
      & 6.691788327147559E-02_wp, 4.324296291121660E-03_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-1.349004251090167E-05_wp,&
      & 4.007744308824233E-05_wp, 2.274433359894836E-05_wp, 2.103612812962821E-03_wp,&
      &-4.072014513167806E-04_wp,-6.093185202772807E-06_wp, 4.478030745477695E-07_wp,&
      &-2.310907317698089E-04_wp, 6.503149886306556E-06_wp, 4.625269223552698E-03_wp,&
      &-4.654914732904501E-03_wp,-7.810260829800153E-03_wp,-1.073001727925391E-02_wp,&
      & 1.034734669126724E-03_wp, 7.137127420796189E-04_wp,-1.444390076936583E-04_wp,&
      & 1.736132264333461E-03_wp, 1.022307758041926E-03_wp, 6.047256261536093E-03_wp,&
      & 9.226408301780580E-03_wp, 2.286344247452107E-03_wp, 1.201607135749007E-02_wp,&
      & 2.017654828951673E-03_wp, 3.228777988126517E-02_wp, 6.236274278382165E-02_wp,&
      &-4.582692336889596E-03_wp,-1.322267895469085E-02_wp, 2.439371290099738E-03_wp,&
      &-7.929227533829630E-03_wp,-4.096922349440278E-02_wp,-3.845808049292056E-02_wp,&
      &-1.084628766033663E-02_wp, 1.936340914448410E-02_wp, 9.067390599316298E-03_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      &-1.821272091314209E-06_wp, 9.601871392670786E-06_wp, 6.093185202772808E-06_wp,&
      & 4.072014513167806E-04_wp,-6.144729534573771E-07_wp,-2.538703927610979E-08_wp,&
      & 2.722175039254806E-09_wp,-3.626020197309276E-07_wp, 2.325732235599854E-08_wp,&
      & 2.131661632649750E-04_wp,-3.423087656910362E-05_wp,-7.137127420796187E-04_wp,&
      &-1.034734669126724E-03_wp, 1.149574814148520E-06_wp, 7.293465151063494E-07_wp,&
      &-3.166311608540729E-07_wp, 2.892653938852755E-06_wp, 2.166188392389130E-06_wp,&
      & 1.168083121507140E-02_wp, 1.782165554648693E-02_wp, 1.524869736605941E-02_wp,&
      & 2.286344247452108E-03_wp, 3.228777988126517E-02_wp, 4.766883667501617E-02_wp,&
      & 3.014243757627404E-02_wp, 6.487439978089045E-02_wp, 2.660368501404789E-02_wp,&
      &-2.356763940357238E-02_wp,-4.796279121035552E-02_wp,-4.262713084311716E-02_wp,&
      & 8.974546488196553E-03_wp,-3.899407999281289E-02_wp,-1.401759762715071E-02_wp,&
      & 1.555359763288199E-03_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 3.708990204362797E-04_wp,-1.955405076968439E-03_wp,&
      &-8.336356213639305E-04_wp, 6.093185202772807E-06_wp,-2.538703927610979E-08_wp,&
      & 4.555430545465183E-06_wp,-1.182538937280366E-06_wp,-1.320006357409028E-08_wp,&
      &-4.373639282681602E-06_wp, 1.293491177497988E-04_wp,-2.077127821882286E-05_wp,&
      &-2.916231196675856E-04_wp,-7.137127420796189E-04_wp, 7.293465151063494E-07_wp,&
      & 3.901855150951473E-07_wp,-4.204276800721357E-07_wp, 1.618024824302476E-06_wp,&
      & 1.446248786416078E-06_wp, 5.319101370119733E-03_wp,-4.754718969258193E-03_wp,&
      & 1.878768563866048E-02_wp,-3.301565519630095E-03_wp, 6.236274278382165E-02_wp,&
      & 3.014243757627402E-02_wp,-6.084476365669887E-02_wp, 2.093020286222214E-02_wp,&
      &-2.325394870315194E-02_wp,-9.573576200081018E-03_wp,-4.567622608156515E-02_wp,&
      & 1.597426822799090E-02_wp, 4.978089179485133E-02_wp, 1.109215781707619E-02_wp,&
      & 1.703184579190498E-02_wp, 3.748857784386515E-02_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-6.714094920657612E-05_wp,&
      & 1.606921694867899E-04_wp, 4.914846471812638E-04_wp,-4.478030745477695E-07_wp,&
      & 2.722175039254806E-09_wp,-1.182538937280366E-06_wp,-2.868236394280212E-07_wp,&
      & 3.295397489360130E-09_wp, 4.884164944858348E-07_wp,-6.710463509118858E-05_wp,&
      & 5.223671899238223E-05_wp, 4.987385488180778E-04_wp, 1.444390076936583E-04_wp,&
      &-3.166311608540729E-07_wp,-4.204276800721357E-07_wp,-3.164540880755009E-07_wp,&
      &-1.162518590140049E-06_wp,-3.805009412713238E-07_wp, 8.110895686925459E-03_wp,&
      & 2.286344247452108E-03_wp, 1.058833839992669E-02_wp, 1.611658198895281E-02_wp,&
      &-4.582692336889597E-03_wp, 6.487439978089043E-02_wp, 2.093020286222212E-02_wp,&
      &-7.121644583707671E-04_wp,-5.394287007897097E-02_wp, 3.006495495986367E-02_wp,&
      &-1.101216010416653E-02_wp,-3.899407999281288E-02_wp, 6.231714940740240E-03_wp,&
      &-1.354674420833593E-02_wp, 5.004969665344108E-02_wp, 1.293692249384755E-02_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      &-1.033589391621077E-06_wp, 6.093185202772807E-06_wp, 2.323103829461939E-06_wp,&
      & 2.310907317698089E-04_wp,-3.626020197309276E-07_wp,-1.320006357409028E-08_wp,&
      & 3.295397489360130E-09_wp,-1.813175951776987E-07_wp, 1.420942667549927E-08_wp,&
      & 3.576614032086451E-04_wp,-7.137127420796188E-04_wp,-8.063629346907818E-04_wp,&
      &-1.736132264333461E-03_wp, 2.892653938852755E-06_wp, 1.618024824302476E-06_wp,&
      &-1.162518590140049E-06_wp, 4.279000292432886E-06_wp, 3.270087771174437E-06_wp,&
      &-2.254913440674616E-03_wp,-1.146804083903603E-02_wp,-8.525367787670529E-04_wp,&
      & 7.080409322663008E-03_wp,-1.322267895469086E-02_wp, 2.660368501404787E-02_wp,&
      &-2.325394870315192E-02_wp,-5.394287007897095E-02_wp,-2.851259805292793E-02_wp,&
      & 1.272603818558988E-02_wp, 1.936340914448410E-02_wp,-9.913117013673664E-03_wp,&
      & 1.434032871363807E-02_wp, 4.032125164772303E-02_wp, 3.677954743243225E-02_wp,&
      & 1.011079258997202E-04_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp,-3.267752137599800E-04_wp, 1.519164074610300E-03_wp,&
      & 1.093247904369124E-03_wp,-6.503149886306556E-06_wp, 2.325732235599854E-08_wp,&
      &-4.373639282681602E-06_wp, 4.884164944858348E-07_wp, 1.420942667549927E-08_wp,&
      & 3.558259083686315E-06_wp, 2.561653298769958E-04_wp,-6.526348981436349E-04_wp,&
      &-8.576804930573188E-04_wp,-1.022307758041926E-03_wp, 2.166188392389130E-06_wp,&
      & 1.446248786416078E-06_wp,-3.805009412713238E-07_wp, 3.270087771174437E-06_wp,&
      & 1.950141181023217E-06_wp, 2.545368264525664E-01_wp, 3.278091162485747E-01_wp,&
      & 3.340018745778095E-02_wp,-1.384415008471236E-01_wp,-1.566460435169708E-01_wp,&
      & 3.779218793477303E-02_wp,-1.239481188122672E-01_wp,-1.596052995066541E-02_wp,&
      &-1.523797436792692E-01_wp,-3.275114472081729E-02_wp,-1.484974722030525E-02_wp,&
      &-6.203146574730595E-02_wp,-1.595284220086066E-02_wp, 2.619734715154162E-02_wp,&
      &-1.444531010377486E-02_wp, 8.394360253427902E-02_wp, 1.321302646831464E-03_wp,&
      & 4.840842595230028E-03_wp, 2.747224633155594E-03_wp,-1.349004251090167E-05_wp,&
      &-1.821272091314209E-06_wp, 3.708990204362797E-04_wp,-6.714094920657612E-05_wp,&
      &-1.033589391621077E-06_wp,-3.267752137599800E-04_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.328190143562808E-03_wp,&
      & 3.216156559734708E-03_wp,-6.041883010632351E-06_wp,-4.573897085653133E-03_wp,&
      &-4.089499315072396E-04_wp,-5.402018447513957E-07_wp,-2.509008997861063E-04_wp,&
      & 7.682547778633844E-07_wp, 1.470191147200401E-04_wp,-1.479529014292982E-01_wp,&
      & 3.270230911069708E-02_wp,-1.268353560535391E-02_wp, 5.257239072303659E-02_wp,&
      & 1.660461878936972E-01_wp,-4.006005257356780E-02_wp, 9.154332634880823E-02_wp,&
      & 2.396835040534949E-02_wp, 6.499308133004139E-02_wp,-1.080899641447741E-02_wp,&
      & 2.092008883615449E-02_wp, 4.327699079137379E-02_wp, 1.791315680943756E-02_wp,&
      &-2.886620266640593E-02_wp, 1.197635916854938E-02_wp,-8.003771584250011E-02_wp,&
      &-4.840842595230028E-03_wp,-1.227789019756612E-02_wp,-8.161704368013259E-03_wp,&
      & 4.007744308824233E-05_wp, 9.601871392670786E-06_wp,-1.955405076968439E-03_wp,&
      & 1.606921694867899E-04_wp, 6.093185202772807E-06_wp, 1.519164074610300E-03_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      &-3.216156559734708E-03_wp,-4.201361488852514E-03_wp, 1.186321088327884E-05_wp,&
      & 8.980826935912253E-03_wp, 1.207714942072668E-03_wp, 1.595329377454439E-06_wp,&
      & 8.209367842672483E-04_wp,-2.992657852805844E-06_wp,-8.436287474274741E-04_wp,&
      &-1.507479321872774E-02_wp,-1.268353560535391E-02_wp, 1.558936851198954E-01_wp,&
      & 5.356555441007388E-03_wp, 2.396835040534950E-02_wp, 1.580570462019544E-01_wp,&
      & 3.824122009652010E-02_wp,-6.675120858771076E-02_wp, 2.331556552072368E-02_wp,&
      & 6.049946259626848E-03_wp,-5.616858926330325E-02_wp, 2.943686192571210E-02_wp,&
      &-5.959318406889177E-02_wp,-1.243187923466680E-02_wp,-5.463882165553090E-02_wp,&
      &-1.550645904138615E-02_wp,-2.747224633155594E-03_wp,-8.161704368013261E-03_wp,&
      &-2.528121064115085E-03_wp, 2.274433359894836E-05_wp, 6.093185202772808E-06_wp,&
      &-8.336356213639305E-04_wp, 4.914846471812638E-04_wp, 2.323103829461939E-06_wp,&
      & 1.093247904369124E-03_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 6.041883010632352E-06_wp, 1.186321088327884E-05_wp,&
      & 2.113525572614489E-03_wp,-1.687141302909504E-05_wp,-2.992657852805846E-06_wp,&
      & 2.709268809615223E-04_wp,-2.423779426335597E-06_wp,-3.853020362159290E-04_wp,&
      & 1.075872311698144E-06_wp, 6.248399057037403E-02_wp, 5.257239072303660E-02_wp,&
      & 5.356555441007389E-03_wp, 1.349834423346564E-01_wp, 6.449243520200129E-02_wp,&
      & 2.396835040534949E-02_wp,-3.866090009118883E-02_wp, 6.571078467879551E-03_wp,&
      &-1.658347532571246E-01_wp,-9.722842322371636E-02_wp, 2.842922138796015E-03_wp,&
      &-2.886620266640593E-02_wp,-7.565147491895572E-03_wp,-1.288305595534449E-02_wp,&
      &-1.706282353260835E-02_wp,-5.913517130087542E-03_wp, 1.349004251090167E-05_wp,&
      & 4.007744308824233E-05_wp, 2.274433359894836E-05_wp, 2.103612812962821E-03_wp,&
      & 4.072014513167806E-04_wp, 6.093185202772807E-06_wp,-4.478030745477695E-07_wp,&
      & 2.310907317698089E-04_wp,-6.503149886306556E-06_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 4.573897085653133E-03_wp,&
      & 8.980826935912253E-03_wp,-1.687141302909504E-05_wp,-1.065864695998768E-02_wp,&
      &-1.994606086619037E-03_wp,-2.992657852805845E-06_wp,-1.167505466641572E-03_wp,&
      & 3.747074000850726E-06_wp, 4.291618080490168E-04_wp,-1.017180391830452E-02_wp,&
      & 1.233135000684818E-02_wp, 4.817287462261997E-04_wp,-2.000041929617723E-02_wp,&
      &-2.309479251972242E-03_wp,-1.948093548852921E-02_wp, 2.674150463780352E-02_wp,&
      &-7.477946722906489E-04_wp, 9.714846767178610E-02_wp, 7.285918079573005E-02_wp,&
      & 3.228646134519803E-03_wp, 2.963744056419990E-02_wp, 9.952273555086245E-03_wp,&
      & 2.838918112505758E-03_wp, 1.795140989110575E-02_wp,-2.817063608821062E-02_wp,&
      &-1.821272091314209E-06_wp,-9.601871392670786E-06_wp,-6.093185202772808E-06_wp,&
      &-4.072014513167806E-04_wp,-6.144729534573771E-07_wp,-2.538703927610979E-08_wp,&
      & 2.722175039254806E-09_wp,-3.626020197309276E-07_wp, 2.325732235599854E-08_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      &-4.089499315072396E-04_wp,-1.207714942072668E-03_wp, 2.992657852805846E-06_wp,&
      & 1.994606086619037E-03_wp, 5.727149785103583E-06_wp, 7.907856893085177E-09_wp,&
      & 3.558796095575019E-06_wp,-1.179230025999622E-08_wp,-2.358855833689911E-06_wp,&
      & 2.454034054646125E-03_wp,-2.975042883210660E-03_wp,-1.811990762199503E-02_wp,&
      & 4.817287462262000E-04_wp,-1.948093548852921E-02_wp,-7.835668817848246E-02_wp,&
      &-2.537132714397605E-02_wp, 4.282141048512239E-02_wp,-1.251459174206802E-02_wp,&
      &-2.226924160852675E-03_wp, 4.814976599275308E-02_wp,-2.696830351838712E-02_wp,&
      & 4.337691668692271E-02_wp, 1.390006426015647E-02_wp, 3.296613799925561E-02_wp,&
      & 1.736265410576017E-02_wp, 3.708990204362797E-04_wp, 1.955405076968439E-03_wp,&
      & 8.336356213639305E-04_wp,-6.093185202772807E-06_wp,-2.538703927610979E-08_wp,&
      & 4.555430545465183E-06_wp,-1.182538937280366E-06_wp,-1.320006357409028E-08_wp,&
      &-4.373639282681602E-06_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp,-5.402018447513957E-07_wp,-1.595329377454439E-06_wp,&
      &-2.709268809615223E-04_wp, 2.992657852805845E-06_wp, 7.907856893085177E-09_wp,&
      &-2.593382717908868E-07_wp, 5.625907935246967E-09_wp, 4.042412949936638E-07_wp,&
      &-3.649930827548552E-09_wp,-8.048565621540937E-03_wp, 1.413549772308417E-02_wp,&
      &-1.736985804602402E-03_wp,-5.969753197836450E-03_wp, 2.674150463780350E-02_wp,&
      &-2.537132714397605E-02_wp, 4.308469549368484E-02_wp, 1.071490826274629E-02_wp,&
      & 2.601319210379564E-02_wp, 7.860618424361802E-03_wp, 1.460196370138135E-02_wp,&
      & 3.899632680726239E-02_wp, 1.897528876075850E-02_wp,-1.646906612148191E-02_wp,&
      & 1.420427503992098E-02_wp,-2.014734551457818E-02_wp,-6.714094920657612E-05_wp,&
      &-1.606921694867899E-04_wp,-4.914846471812638E-04_wp, 4.478030745477695E-07_wp,&
      & 2.722175039254806E-09_wp,-1.182538937280366E-06_wp,-2.868236394280212E-07_wp,&
      & 3.295397489360130E-09_wp, 4.884164944858348E-07_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-2.509008997861063E-04_wp,&
      &-8.209367842672483E-04_wp, 2.423779426335597E-06_wp, 1.167505466641572E-03_wp,&
      & 3.558796095575019E-06_wp, 5.625907935246967E-09_wp, 2.208314275925427E-06_wp,&
      &-8.000955000555548E-09_wp,-1.279401244822832E-06_wp,-1.036396307531422E-03_wp,&
      & 4.817287462261997E-04_wp, 7.652463223438904E-03_wp,-2.037825431370820E-03_wp,&
      &-7.477946722906489E-04_wp, 4.282141048512239E-02_wp, 1.071490826274628E-02_wp,&
      & 4.953626045965858E-03_wp, 1.451154029554059E-02_wp, 1.195187549207984E-02_wp,&
      &-9.891343763986659E-04_wp, 1.390006426015647E-02_wp,-1.831909227229836E-02_wp,&
      & 7.468054342921150E-05_wp,-3.380965573848588E-02_wp,-3.036490348359322E-03_wp,&
      &-1.033589391621077E-06_wp,-6.093185202772807E-06_wp,-2.323103829461939E-06_wp,&
      &-2.310907317698089E-04_wp,-3.626020197309276E-07_wp,-1.320006357409028E-08_wp,&
      & 3.295397489360130E-09_wp,-1.813175951776987E-07_wp, 1.420942667549927E-08_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 7.682547778633844E-07_wp, 2.992657852805844E-06_wp, 3.853020362159290E-04_wp,&
      &-3.747074000850726E-06_wp,-1.179230025999622E-08_wp, 4.042412949936638E-07_wp,&
      &-8.000955000555548E-09_wp,-5.499909544027137E-07_wp, 3.671915150243327E-09_wp,&
      &-9.894771926742629E-03_wp, 2.260288857247389E-02_wp, 4.686087259199180E-04_wp,&
      & 5.661029063890270E-03_wp, 9.714846767178610E-02_wp,-1.251459174206803E-02_wp,&
      & 2.601319210379562E-02_wp, 1.451154029554059E-02_wp,-7.675290382092234E-03_wp,&
      &-3.655550738162818E-02_wp, 1.795140989110575E-02_wp, 1.781922459229878E-02_wp,&
      & 9.681220535810813E-03_wp,-2.331090475596502E-02_wp, 2.237134098318999E-03_wp,&
      &-6.931805241115380E-02_wp,-3.267752137599800E-04_wp,-1.519164074610300E-03_wp,&
      &-1.093247904369124E-03_wp, 6.503149886306556E-06_wp, 2.325732235599854E-08_wp,&
      &-4.373639282681602E-06_wp, 4.884164944858348E-07_wp, 1.420942667549927E-08_wp,&
      & 3.558259083686315E-06_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 1.470191147200401E-04_wp, 8.436287474274741E-04_wp,&
      &-1.075872311698144E-06_wp,-4.291618080490168E-04_wp,-2.358855833689911E-06_wp,&
      &-3.649930827548552E-09_wp,-1.279401244822832E-06_wp, 3.671915150243327E-09_wp,&
      & 1.374932036159839E-08_wp, 2.546290040058253E-01_wp,-1.928389975975781E-02_wp,&
      & 3.405990142396841E-02_wp, 3.552549959260523E-01_wp,-2.366042499400333E-02_wp,&
      &-2.268431836811374E-03_wp,-1.238855387564933E-01_wp, 4.178987409106263E-02_wp,&
      & 2.172982997449486E-01_wp,-1.459634274178628E-02_wp,-2.287485624110323E-03_wp,&
      & 3.644875522761238E-03_wp,-1.626917205522780E-02_wp,-6.714721892983819E-02_wp,&
      & 2.100836045574692E-02_wp, 8.892806632054290E-02_wp, 1.347613073731518E-03_wp,&
      & 1.672739882063343E-03_wp, 2.806610975486558E-03_wp, 4.625269223552698E-03_wp,&
      & 2.131661632649750E-04_wp, 1.293491177497988E-04_wp,-6.710463509118858E-05_wp,&
      & 3.576614032086451E-04_wp, 2.561653298769958E-04_wp, 1.328190143562808E-03_wp,&
      &-3.216156559734708E-03_wp, 6.041883010632351E-06_wp, 4.573897085653133E-03_wp,&
      &-4.089499315072396E-04_wp,-5.402018447513957E-07_wp,-2.509008997861063E-04_wp,&
      & 7.682547778633844E-07_wp, 1.470191147200401E-04_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 8.702694713955241E-03_wp,&
      & 1.567901565919992E-01_wp, 7.604487167628128E-04_wp, 7.931708386725814E-03_wp,&
      & 1.755955375970396E-01_wp, 1.683513749174361E-02_wp,-5.374298642707933E-03_wp,&
      & 3.691329429922923E-03_wp, 2.883921128113612E-02_wp, 9.855306067937845E-02_wp,&
      & 1.547578577064782E-02_wp,-2.481186422430601E-02_wp,-1.074097956085215E-03_wp,&
      &-4.346579276711809E-03_wp, 3.093775544053226E-03_wp, 2.038431321445513E-02_wp,&
      &-1.672739882063343E-03_wp, 4.577598973848464E-04_wp,-2.824599855246675E-03_wp,&
      &-4.654914732904501E-03_wp,-3.423087656910362E-05_wp,-2.077127821882286E-05_wp,&
      & 5.223671899238223E-05_wp,-7.137127420796188E-04_wp,-6.526348981436349E-04_wp,&
      & 3.216156559734708E-03_wp,-4.201361488852514E-03_wp, 1.186321088327884E-05_wp,&
      & 8.980826935912253E-03_wp,-1.207714942072668E-03_wp,-1.595329377454439E-06_wp,&
      &-8.209367842672483E-04_wp, 2.992657852805844E-06_wp, 8.436287474274741E-04_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      &-1.537100523094241E-02_wp, 7.604487167628129E-04_wp, 1.558775731802586E-01_wp,&
      &-1.400926208604901E-02_wp, 3.691329429922923E-03_wp,-9.291192814470115E-03_wp,&
      & 3.899864345132840E-02_wp, 1.711657240793602E-01_wp,-3.390131872627125E-02_wp,&
      & 2.749421233824322E-03_wp,-8.469260304161211E-03_wp,-1.765906564399253E-03_wp,&
      &-5.945723090416469E-02_wp, 3.253217124943841E-02_wp, 7.778202904884632E-02_wp,&
      &-1.675082026709908E-02_wp,-2.806610975486558E-03_wp,-2.824599855246675E-03_wp,&
      &-2.598041086701424E-03_wp,-7.810260829800153E-03_wp,-7.137127420796187E-04_wp,&
      &-2.916231196675856E-04_wp, 4.987385488180778E-04_wp,-8.063629346907818E-04_wp,&
      &-8.576804930573188E-04_wp,-6.041883010632352E-06_wp, 1.186321088327884E-05_wp,&
      & 2.113525572614489E-03_wp,-1.687141302909504E-05_wp, 2.992657852805846E-06_wp,&
      &-2.709268809615223E-04_wp, 2.423779426335597E-06_wp, 3.853020362159290E-04_wp,&
      &-1.075872311698144E-06_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp,-1.603241986148280E-01_wp, 7.931708386725814E-03_wp,&
      &-1.400926208604901E-02_wp, 1.109989558061869E-02_wp, 2.885657677088460E-02_wp,&
      & 3.691329429922923E-03_wp, 9.900727893249391E-02_wp,-5.096749996079889E-02_wp,&
      &-1.759154509479025E-01_wp, 1.777688848623613E-02_wp, 3.640903566789722E-03_wp,&
      &-4.346579276711810E-03_wp, 1.978742213800154E-02_wp, 5.502645539805806E-02_wp,&
      &-2.555517665681596E-02_wp,-7.460636655623726E-02_wp,-4.625269223552698E-03_wp,&
      &-4.654914732904501E-03_wp,-7.810260829800153E-03_wp,-1.073001727925391E-02_wp,&
      &-1.034734669126724E-03_wp,-7.137127420796189E-04_wp, 1.444390076936583E-04_wp,&
      &-1.736132264333461E-03_wp,-1.022307758041926E-03_wp,-4.573897085653133E-03_wp,&
      & 8.980826935912253E-03_wp,-1.687141302909504E-05_wp,-1.065864695998768E-02_wp,&
      & 1.994606086619037E-03_wp, 2.992657852805845E-06_wp, 1.167505466641572E-03_wp,&
      &-3.747074000850726E-06_wp,-4.291618080490168E-04_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,-1.532338905733210E-03_wp,&
      &-1.954549296314300E-02_wp, 7.459203931459594E-05_wp, 1.836690509563453E-03_wp,&
      &-9.995084309531216E-02_wp,-1.186547464654780E-02_wp, 4.024635858893756E-03_wp,&
      &-3.366551588210808E-03_wp,-2.091604547977351E-02_wp,-7.983685125583710E-02_wp,&
      &-1.510885403830548E-02_wp, 1.989155217895337E-02_wp, 1.531204046169162E-03_wp,&
      & 5.022744764462081E-03_wp,-3.942121353900568E-03_wp,-2.116930580241157E-02_wp,&
      & 2.131661632649750E-04_wp, 3.423087656910362E-05_wp, 7.137127420796187E-04_wp,&
      & 1.034734669126724E-03_wp, 1.149574814148520E-06_wp, 7.293465151063494E-07_wp,&
      &-3.166311608540729E-07_wp, 2.892653938852755E-06_wp, 2.166188392389130E-06_wp,&
      &-4.089499315072396E-04_wp, 1.207714942072668E-03_wp,-2.992657852805846E-06_wp,&
      &-1.994606086619037E-03_wp, 5.727149785103583E-06_wp, 7.907856893085177E-09_wp,&
      & 3.558796095575019E-06_wp,-1.179230025999622E-08_wp,-2.358855833689911E-06_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      &-1.469122536653885E-04_wp,-1.873914712647955E-03_wp, 1.065824624332955E-03_wp,&
      & 7.459203931459594E-05_wp,-1.186547464654780E-02_wp, 2.267202257691398E-02_wp,&
      & 1.521968071603357E-03_wp, 6.455916960609892E-03_wp,-2.661248062415271E-03_wp,&
      &-1.078538960739731E-02_wp, 2.497148719603814E-02_wp, 5.949286951402333E-03_wp,&
      &-2.542208638003138E-03_wp, 2.139810523794881E-03_wp, 6.112772769082304E-03_wp,&
      &-2.802806838903087E-03_wp, 1.293491177497988E-04_wp, 2.077127821882286E-05_wp,&
      & 2.916231196675856E-04_wp, 7.137127420796189E-04_wp, 7.293465151063494E-07_wp,&
      & 3.901855150951473E-07_wp,-4.204276800721357E-07_wp, 1.618024824302476E-06_wp,&
      & 1.446248786416078E-06_wp,-5.402018447513957E-07_wp, 1.595329377454439E-06_wp,&
      & 2.709268809615223E-04_wp,-2.992657852805845E-06_wp, 7.907856893085177E-09_wp,&
      &-2.593382717908868E-07_wp, 5.625907935246967E-09_wp, 4.042412949936638E-07_wp,&
      &-3.649930827548552E-09_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp,-8.023297592600592E-03_wp,-8.323523586443188E-04_wp,&
      &-1.768572441644346E-03_wp, 1.533389705729007E-02_wp, 4.024635858893757E-03_wp,&
      & 1.521968071603357E-03_wp, 4.293684846532067E-02_wp,-2.803824785510203E-02_wp,&
      &-3.696241844564566E-02_wp, 3.485514967844202E-03_wp, 2.247471923823931E-03_wp,&
      &-2.282711389178954E-03_wp, 1.933726609315561E-02_wp, 4.205293718418014E-02_wp,&
      &-2.064087301454750E-02_wp,-2.123546368463536E-02_wp,-6.710463509118858E-05_wp,&
      &-5.223671899238223E-05_wp,-4.987385488180778E-04_wp,-1.444390076936583E-04_wp,&
      &-3.166311608540729E-07_wp,-4.204276800721357E-07_wp,-3.164540880755009E-07_wp,&
      &-1.162518590140049E-06_wp,-3.805009412713238E-07_wp,-2.509008997861063E-04_wp,&
      & 8.209367842672483E-04_wp,-2.423779426335597E-06_wp,-1.167505466641572E-03_wp,&
      & 3.558796095575019E-06_wp, 5.625907935246967E-09_wp, 2.208314275925427E-06_wp,&
      &-8.000955000555548E-09_wp,-1.279401244822832E-06_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 1.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 2.706470824241631E-03_wp,&
      & 7.459203931459591E-05_wp,-1.963500781960325E-02_wp,-3.244027322347748E-03_wp,&
      &-3.366551588210809E-03_wp, 6.455916960609892E-03_wp,-2.803824785510204E-02_wp,&
      &-9.591078218594690E-02_wp, 2.485883202263235E-02_wp,-3.634650375220897E-03_wp,&
      & 8.103144802471892E-03_wp, 2.139810523794881E-03_wp, 4.683348962545748E-02_wp,&
      &-3.335492577637602E-02_wp,-6.163884273211701E-02_wp, 1.842507427669205E-02_wp,&
      & 3.576614032086451E-04_wp, 7.137127420796188E-04_wp, 8.063629346907818E-04_wp,&
      & 1.736132264333461E-03_wp, 2.892653938852755E-06_wp, 1.618024824302476E-06_wp,&
      &-1.162518590140049E-06_wp, 4.279000292432886E-06_wp, 3.270087771174437E-06_wp,&
      & 7.682547778633844E-07_wp,-2.992657852805844E-06_wp,-3.853020362159290E-04_wp,&
      & 3.747074000850726E-06_wp,-1.179230025999622E-08_wp, 4.042412949936638E-07_wp,&
      &-8.000955000555548E-09_wp,-5.499909544027137E-07_wp, 3.671915150243327E-09_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 1.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.407306246330118E-02_wp,-6.708106897412503E-04_wp,-6.850563048498959E-04_wp,&
      &-2.664860363771760E-02_wp,-2.091604547977351E-02_wp,-2.661248062415271E-03_wp,&
      &-3.696241844564562E-02_wp, 2.485883202263235E-02_wp, 8.986553028404413E-02_wp,&
      &-2.067213199823943E-02_wp,-3.942121353900569E-03_wp, 4.149032133885584E-03_wp,&
      &-1.406263986718216E-02_wp,-3.598740316775212E-02_wp, 2.066651179470221E-02_wp,&
      & 7.067773510318530E-02_wp, 2.561653298769958E-04_wp, 6.526348981436349E-04_wp,&
      & 8.576804930573188E-04_wp, 1.022307758041926E-03_wp, 2.166188392389130E-06_wp,&
      & 1.446248786416078E-06_wp,-3.805009412713238E-07_wp, 3.270087771174437E-06_wp,&
      & 1.950141181023217E-06_wp, 1.470191147200401E-04_wp,-8.436287474274741E-04_wp,&
      & 1.075872311698144E-06_wp, 4.291618080490168E-04_wp,-2.358855833689911E-06_wp,&
      &-3.649930827548552E-09_wp,-1.279401244822832E-06_wp, 3.671915150243327E-09_wp,&
      & 1.374932036159839E-08_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp,&
      & 1.000000000000000E+00_wp],shape(overlap))

      type(structure_type) :: mol

      call get_structure(mol, "f-block", "AcCl3")
      call test_qvszp_overlap(error, mol, overlap)

   end subroutine test_qvszp_overlap_accl3

end module test_qvszp
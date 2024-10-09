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

!> @file tblite/ceh/coupled_perturbed.f90
!> Implementation of the coupled-perturbed equations for the
!> density matrix and Mulliken charge derivatives
module tblite_ceh_coupled_perturbed
   use iso_fortran_env, only: output_unit
   
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
   use tblite_wavefunction_fermi, only : get_fermi_filling_gradient
   use tblite_blas, only : gemm
   use tblite_lapack, only : getrf, getrs
   implicit none
   private

   public :: get_density_matrix_gradient

contains


   !> Determine the density matrix gradient based on the coupled-perturbed formalism
   subroutine get_density_matrix_gradient(mol,bas,wfn,list,dh0dr,dh0dL,doverlap,ddensitydr,ddensitydL)
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
      real(wp), intent(in) :: dh0dL(:, :, :, :)
      !> Derivative of the electronic energy w.r.t. coordinate displacements
      real(wp), intent(in) :: doverlap(:, :, :)
      !> Derivative of the electronic energy w.r.t. coordinate displacements
      real(wp), intent(out) :: ddensitydr(:, :, :)
      !> Derivative of the electronic energy w.r.t. the lattice vector
      real(wp), intent(out) :: ddensitydL(:, :, :, :)

      integer :: ic, iao, jao, homo, spin
      real(wp) :: denom, e_fermi
      real(wp), allocatable :: dh0dr_mo(:,:), doverlap_mo(:,:), ddensitydr_mo(:,:), tmp(:,:)
      real(wp), allocatable :: u_mo(:,:), demo(:), a(:, :), b(:, :)
      real(wp), allocatable :: focc(:), dfocc(:)
      integer, allocatable :: ipiv(:)
      integer :: info

      allocate(dh0dr_mo(bas%nao,bas%nao), doverlap_mo(bas%nao,bas%nao), ddensitydr_mo(bas%nao,bas%nao), &
      & tmp(bas%nao,bas%nao), u_mo(bas%nao,bas%nao), demo(bas%nao), a(bas%nao, bas%nao), b(bas%nao, bas%nao), &
      & focc(bas%nao), dfocc(bas%nao), source = 0.0_wp)

      allocate(ipiv(bas%nao))

      ddensitydr(:, :, :) = 0.0_wp
      ddensitydL(:, :, :, :) = 0.0_wp

      do ic = 1, 3

         ddensitydr_mo = 0.0_wp
         demo = 0.0_wp
         tmp = 0.0_wp
         u_mo = 0.0_wp

         ! Transform the derivatives of Fock and overlap matrix into the MO basis: X_mo^(1) = C^T x X^(1) x C
         ! Fock matrix derivative
         call gemm(amat=wfn%coeff(:,:,1),bmat=dh0dr(ic,:,:),cmat=tmp,transa='T',transb='N')
         call gemm(amat=tmp,bmat=wfn%coeff(:,:,1),cmat=dh0dr_mo(:,:),transa='N',transb='N')
         ! Overlap derivative
         call gemm(amat=wfn%coeff(:,:,1),bmat=doverlap(ic,:,:),cmat=tmp,transa='T',transb='N')
         call gemm(amat=tmp,bmat=wfn%coeff(:,:,1),cmat=doverlap_mo(:,:),transa='N',transb='N')

         !write(*,*) "focc", wfn%focc, wfn%nspin

         do spin = 1, wfn%nspin
            ! Determine the eigen value derivative E^(1)_pp = ( (F_mo^{1}_pp - S_mo^(1)_pp * E^(0)_pp) ) * n^(0)_p
            do iao = 1, bas%nao
               demo(iao) = (dh0dr_mo(iao,iao) - doverlap_mo(iao,iao) * wfn%emo(iao,spin)) * wfn%focc(iao,spin)
            end do 

            ! Determine the occupation number derivative based on the Fermi-Smearing
            call get_fermi_filling_gradient(wfn%nel(spin), wfn%kt, wfn%emo(:,spin), demo, &
            & homo, focc, dfocc, e_fermi)
         
            !write(*,*) "dfocc", dfocc

            ! Solve the Coupled Perturbed equations u^(1)_qp as a linear equations of the form: A_pq * u^(1)_qp = B_qp
            ! A_pq = (n^(0)_q - n^(0)_p) * (E^(0)_pp + E^(0)_qq)
            ! B_qp = (n^(0)_q * F_mo^(1)_qp - S_mo^(1)_qp * (n^(0)_q - n^(0)_p) * E^(0)_pp)
            do iao = 1, bas%nao
               do jao = 1, bas%nao
                  ! Setup A_pq
                  a(iao, jao) = (wfn%focc(jao,spin) - wfn%focc(iao,spin)) * (wfn%emo(iao,spin) + wfn%emo(jao,spin))
                  ! Setup B_qp
                  b(iao, jao) = wfn%focc(jao,spin) * dh0dr_mo(iao,jao) &
                     & - doverlap_mo(iao,jao) * (wfn%focc(jao,spin) - wfn%focc(iao,spin)) * wfn%emo(iao,spin)
               end do 
            end do 

            !call write_2d_matrix(a, "a")
            !call write_2d_matrix(b, "b")

            ! LU decomoposition of the A matrix
            call getrf(a, ipiv, info)
            if (info == 0) then
               ! Generate inverse of the A matrix given its LU decomposition
               call getrs(a, b, ipiv, info) !, trans="t")
            endif

            !call write_2d_matrix(a, "after")
            !call write_2d_matrix(b, "solution")
            
            ! Determine the unitary transformation matrix derivative
            u_mo = b
            ! Find all the non-independent pairs
            do iao = 1, bas%nao
               do jao = 1, bas%nao 
                  !write(*,*) "cut", (wfn%focc(jao,spin) - wfn%focc(iao,spin)) * (wfn%emo(iao,spin) + wfn%emo(jao,spin))
                  if((wfn%focc(jao,spin) - wfn%focc(iao,spin)) * (wfn%emo(iao,spin) + wfn%emo(jao,spin)) == 0.0_wp) then 
                     !write(*,*) "non-independent pair", iao, jao
                     u_mo(iao, jao) = -0.5_wp * doverlap_mo(iao, jao)
                  end if
               end do 
            end do 

            !call write_2d_matrix(u_mo, "u_mo")

            ! ! Determine the unitary transformation matrix derivative from the unperturbed to the perturbed basis matrix
            ! do iao = 1, homo
            !    do jao = iao+1, homo
            !       ! Off-diagonal occ-occ terms: u^(1)_qp = -0.5 S_mo^(1)_qp
            !       u_mo(iao,jao) = -0.5_wp * doverlap_mo(iao,jao)
            !       !u_mo(jao,iao) = -0.5_wp * doverlap_mo(jao,iao)
            !    end do 
            !    do jao = homo + 1, bas%nao
            !       !denom = wfn%emo(jao,1)-wfn%emo(iao,1)
            !       ! occ-virt terms: u^(1)_qp = (F_mo^(1)_qp - S_mo^(1)_qp * E^(0)_pp) / (E^(0)_pp - E^(0)_qq)
            !       denom = wfn%emo(iao,1)-wfn%emo(jao,1)
            !       u_mo(iao,jao) = u_mo(iao,jao) + (dh0dr_mo(iao,jao) - doverlap_mo(iao,jao) * wfn%emo(iao,1))/denom 
            !       u_mo(jao,iao) = u_mo(jao,iao) - (dh0dr_mo(iao,jao) - doverlap_mo(iao,jao) * wfn%emo(iao,1))/denom 
            !       !denom = wfn%emo(jao,1)-wfn%emo(iao,1)
            !       !u_mo(jao,iao) = u_mo(jao,iao) + (dh0dr_mo(jao,iao) - doverlap_mo(jao,iao) * wfn%emo(jao,1))/denom
            !       !u_mo(jao,iao) = (dh0dr_mo(jao,iao) - doverlap_mo(jao,iao) * wfn%emo(iao,1))/denom

            !    end do 
            ! end do 
            ! !call write_2d_matrix(u_mo, "u_mo")
            ! do iao = 1, bas%nao
            !    ! Diagonal terms: u^(1)_pp = -0.5 S_mo^(1)_pp
            !    u_mo(iao,iao) = -0.5_wp * doverlap_mo(iao,iao)
            ! end do 

            ! Calculate the first derivative of the density matrix in the MO basis
            do iao = 1, bas%nao
               do jao = 1, bas%nao
                  ! Unitary transformation matrix derivative terms: P_mo^(1)_qp = u^(1)_qp*n^(0)_p + n^(0)_q*u^(1)_qp
                  ddensitydr_mo(iao,jao) = ddensitydr_mo(iao,jao) &
                     & + u_mo(iao,jao) + u_mo(jao,iao) !*wfn%focc(jao,spin) + wfn%focc(iao,spin)*
               end do
               ! Occupation number derivative terms: P_mo^(1)_pp = n^(1)_p
               ddensitydr_mo(iao,iao) = ddensitydr_mo(iao,iao) + dfocc(iao)
            end do
         
         end do

         ! Transform the first derivative of the density matrix back to the AO basis: P^(1) = C x P_mo^(1) x C^T
         call gemm(amat=wfn%coeff(:,:,1),bmat=ddensitydr_mo(:,:),cmat=tmp,transa='N',transb='N')
         call gemm(amat=tmp,bmat=wfn%coeff(:,:,1),cmat=ddensitydr(ic,:,:),transa='N',transb='T')

      end do

   end subroutine get_density_matrix_gradient


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

end module tblite_ceh_coupled_perturbed

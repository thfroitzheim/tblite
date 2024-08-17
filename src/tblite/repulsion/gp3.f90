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

!> @file tblite/repulsion/gp3.f90
!> Provides a screened Coulomb repulsion interaction as used in GP3-xTB

!> Classical repulsion interaction as used in the GP3-xTB Hamiltonian
module tblite_repulsion_gp3
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_container, only : container_cache
   use tblite_repulsion_type, only : repulsion_type, view, taint
   use tblite_repulsion_cache, only : repulsion_cache
   use tblite_ncoord, only : ncoord_type, new_ncoord
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_repulsion_gp3

   !> Container to evaluate classical repulsion interactions for the GP3 Hamiltonian
   type, public, extends(repulsion_type) :: gp3_repulsion
      !> Exponent for the damping of the repulsion (arithmetic average)
      real(wp), allocatable :: alpha(:, :)
      !> Effective nuclear charge
      real(wp), allocatable :: zeff(:, :)
      !> Linear CN dependence of effective nuclear charge
      real(wp), allocatable :: kcn(:)
      !> Linear atomic charge dependence of effective nuclear charge 
      real(wp), allocatable :: kq(:)
      !> Geometric mean of fitted covalent atomic radii
      real(wp), allocatable :: rcov(:, :) 
      !> Coordination number for modifying the self-energies
      class(ncoord_type), allocatable :: ncoord
   contains
      !> Update repulsion cache
      procedure :: update
      !> Evaluates the damping for classical repulsion
      procedure :: get_damping
      !> Evaluates the derivative of the damping for classical repulsion
      procedure :: get_damping_derivs
   end type gp3_repulsion

   character(len=*), parameter :: label = "screened Coulomb repulsion - GP3-xTB"

contains


subroutine new_repulsion_gp3(self, mol, alpha, zeff, kcn, kq, rexp,&
   & rcov_rep, rcov_cn, exp_cn, cutoff)
   !> Instance of the repulsion container
   type(gp3_repulsion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Exponent for the damping of the repulsion
   real(wp), intent(in) :: alpha(:)
   !> Effective nuclear charge
   real(wp), intent(in) :: zeff(:)
   !> Coefficient for CN dependence
   real(wp), intent(in) :: kcn(:)
   !> Coefficient for charge dependence
   real(wp), intent(in) :: kq(:)
   !> Exponent of the repulsion polynomial
   real(wp), intent(in) :: rexp
   !> Fitted covalent atomic radii for repulsion
   real(wp), intent(in) :: rcov_rep(:)
   !> Fitted covalent atomic radii for repulsion CN
   real(wp), intent(in) :: rcov_cn(:)
   !> Exponent of the repulsion CN
   real(wp), intent(in) :: exp_cn
   !> Real-space cutoff
   real(wp), intent(in), optional :: cutoff

   integer :: isp, izp, jsp, jzp

   self%label = label
   allocate(self%alpha(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%alpha(jsp, isp) = 0.5_wp * (alpha(isp) + alpha(jsp))
      end do
   end do

   allocate(self%zeff(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%zeff(jsp, isp) = zeff(isp)*zeff(jsp) 
      end do
   end do

   allocate(self%kcn(mol%nid))
   self%kcn = kcn

   allocate(self%kq(mol%nid))
   self%kq = kq
   
   allocate(self%rexp(mol%nid, mol%nid))
   self%rexp(:, :) = rexp

   allocate(self%rcov(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%rcov(jsp, isp) = sqrt(rcov_rep(isp)*rcov_rep(jsp))
      end do
   end do

   if (present(cutoff)) self%cutoff = cutoff

   ! Create the coordination number 
   call new_ncoord(self%ncoord, mol, cn_type="erf", &
      & kcn=exp_cn, rcov=rcov_cn)

end subroutine new_repulsion_gp3


!> Update repulsion cache
subroutine update(self, mol, cache, wfn)
   !> Instance of the classical repulsion
   class(gp3_repulsion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Optional auxiliary wavefunction data
   type(wavefunction_type), intent(in), optional :: wfn

   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   type(repulsion_cache), pointer :: ptr
   
   call taint(cache, ptr)
   
   if (.not.allocated(ptr%scaled_zeff)) allocate(ptr%scaled_zeff(mol%nat, mol%nat), source=0.0_wp)
   if (.not.allocated(ptr%scaled_dzeffdr)) then 
      allocate(ptr%scaled_dzeffdr(3, mol%nat, mol%nat, mol%nat), source=0.0_wp)
   end if
   if (.not.allocated(ptr%scaled_dzeffdL)) then 
      allocate(ptr%scaled_dzeffdL(3, 3, mol%nat, mol%nat), source=0.0_wp)
   end if 

   ! Calculate the coordination number required effective nuclear charge scaling
   allocate(cn(mol%nat))
   if (allocated(wfn%dqatdr)) then
      allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
   end if
   call self%ncoord%get_cn(mol, cn, dcndr, dcndL)

   ! Update the scaled effective nuclear charge
   if (allocated(wfn%dqatdr)) then
      call get_scaled_dzeff(self, mol, cn, dcndr, dcndL, &
         & wfn%qat(:,1), wfn%dqatdr(:,:,:,1), wfn%dqatdL(:,:,:,1), &
         & ptr%scaled_zeff, ptr%scaled_dzeffdr, ptr%scaled_dzeffdL)
   else
      call get_scaled_zeff(self, mol, cn, wfn%qat(:,1), ptr%scaled_zeff)
   end if 

end subroutine update


!> Scales effective nuclear charge for GP3 repulsion 
subroutine get_scaled_zeff(self, mol, cn, qat, scaled_zeff)
   !> Instance of the repulsion container
   type(gp3_repulsion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Atom-resolved charges
   real(wp), intent(in) :: qat(:)
   !> Scaled effective nuclear charge
   real(wp), intent(out) :: scaled_zeff(:, :)

   integer :: iat, izp, jat, jzp
   real(wp) :: scalei, scalej
   do iat = 1, mol%nat
      izp = mol%id(iat)
      scalei = 1.0_wp + self%kcn(izp) * cn(iat) - self%kq(izp) * qat(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         scalej = 1.0_wp + self%kcn(jzp) * cn(jat) - self%kq(jzp) * qat(jat)

         scaled_zeff(jat, iat) = self%zeff(jzp, izp) * scalei * scalej
      end do
   end do 

end subroutine get_scaled_zeff

!> Gradient for scaled effective nuclear charge for GP3 repulsion
subroutine get_scaled_dzeff(self, mol, cn, dcndr, dcndL, qat, dqatdr, dqatdL, &
   & scaled_zeff, scaled_dzeffdr, scaled_dzeffdL)
   !> Instance of the repulsion container
   type(gp3_repulsion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Derivative of CN w.r.t. coordinates
   real(wp), intent(in) :: dcndr(:, :, :)
   !> Derivative of CN w.r.t. lattice vectors
   real(wp), intent(in) :: dcndL(:, :, :)
   !> Atom-resolved charges
   real(wp), intent(in) :: qat(:)
   !> Derivative of charges w.r.t. coordinates
   real(wp), intent(in) :: dqatdr(:, :, :)
   !> Derivative of charges w.r.t. lattice vectors
   real(wp), intent(in) :: dqatdL(:, :, :)
   !> Scaled effective nuclear charge
   real(wp), intent(out) :: scaled_zeff(:, :)
   !> Derivative of scaled effective nuclear charge w.r.t. coordinates
   real(wp), intent(out) :: scaled_dzeffdr(:, :, :, :)
   !> Derivative of scaled effective nuclear charge w.r.t. lattice vectors
   real(wp), intent(out) :: scaled_dzeffdL(:, :, :, :)

   integer :: iat, izp, jat, jzp
   real(wp) :: scalei, scalej

   scaled_dzeffdr(:, :, :, :) = 0.0_wp 
   scaled_dzeffdL(:, :, :, :) = 0.0_wp 

   !$omp parallel do default(none) schedule(runtime) reduction(+: scaled_dzeffdr, scaled_dzeffdL) &
   !$omp shared(self, mol, cn, qat, dcndr, dcndL, dqatdr, dqatdL, scaled_zeff) &
   !$omp private(iat, jat, izp, jzp, scalei, scalej)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      scalei = 1.0_wp + self%kcn(izp) * cn(iat) - self%kq(izp) * qat(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         scalej = 1.0_wp + self%kcn(jzp) * cn(jat) - self%kq(jzp) * qat(jat)

         scaled_zeff(jat, iat) = self%zeff(jzp, izp) * scalei * scalej
         
         scaled_dzeffdr(:, :, jat, iat) = scaled_dzeffdr(:, :, jat, iat) + self%zeff(jzp, izp) &
         & * (scalei * (self%kcn(jzp) * dcndr(:, :, jat) - self%kq(jzp) * dqatdr(:, :, jat)) &
         & + (self%kcn(izp) * dcndr(:, :, iat) - self%kq(izp) * dqatdr(:, :, iat)) * scalej)
         
         scaled_dzeffdL(:, :, jat, iat) = scaled_dzeffdL(:, :, jat, iat) + self%zeff(jzp, izp) &
         & * (scalei * (self%kcn(jzp) * dcndL(:, :, jat) - self%kq(jzp) * dqatdL(:, :, jat)) &
         & + (self%kcn(izp) * dcndL(:, :, iat) - self%kq(izp) * dqatdL(:, :, iat)) * scalej)

      end do
   end do 

end subroutine get_scaled_dzeff


!> Evaluates the damping for GP3 repulsion
elemental function get_damping(self, izp, jzp, r) result(damp)
   !> Instance of the repulsion container
   class(gp3_repulsion), intent(in) :: self
   !> Atom i index
   integer, intent(in)  :: izp
   !> Atom j index
   integer, intent(in)  :: jzp
   !> Current distance
   real(wp), intent(in) :: r

   real(wp) :: damp, r1rel

   r1rel = (r - self%rcov(jzp, izp)) / self%rcov(jzp, izp)
   damp = 0.5_wp * (1.0_wp + erf(-self%alpha(jzp, izp) * r1rel))

end function get_damping

!> Evaluates the derivative of the damping for GP3 repulsion
elemental function get_damping_derivs(self, izp, jzp, r) result(ddamp)
   !> Instance of the repulsion container
   class(gp3_repulsion), intent(in) :: self
   !> Atom i index
   integer, intent(in)  :: izp
   !> Atom j index
   integer, intent(in)  :: jzp
   !> Current distance
   real(wp), intent(in) :: r

   real(wp) :: ddamp, r1rel, prefactor, expo

   r1rel = (r - self%rcov(jzp, izp)) / self%rcov(jzp, izp)
   prefactor = -self%alpha(jzp, izp) / (sqrt(pi) * self%rcov(jzp, izp))
   expo = self%alpha(jzp, izp) * r1rel
   ddamp = prefactor * exp(-expo**2)

end function get_damping_derivs

end module tblite_repulsion_gp3
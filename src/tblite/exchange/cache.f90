module tblite_exchange_cache
   use mctc_env, only : wp, dp
   use mctc_io, only : structure_type
   use tblite_coulomb_ewald, only : get_alpha
   use tblite_wignerseitz, only : wignerseitz_cell, new_wignerseitz_cell
   implicit none
   private

   public :: exchange_cache


   type :: exchange_cache
      real(dp), allocatable :: curr_D(:, :)
      real(dp), allocatable :: ref_D(:, :)
      real(dp), allocatable :: prev_F(:, :)
      real(dp), allocatable :: gamma_(:,:)
   contains
      procedure :: update
   end type exchange_cache


contains


subroutine update(self, mol)
   !> Instance of the electrostatic container
   class(exchange_cache), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol


end subroutine update

end module tblite_exchange_cache

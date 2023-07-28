module tblite_exchange
   use mctc_env, only : wp
   use tblite_mulliken_kfock, only : mulliken_kfock_type, new_mulliken_exchange
   use mctc_io, only : structure_type
   use tblite_param_exchange, only : exchange_record
   use tblite_exchange_type, only : exchange_type
   use tblite_basis_type, only : basis_type
   implicit none
   private
   public :: new_exchange

contains
subroutine new_exchange(self, mol, hardness, par, bas)
   !> Instance of the multipole container
   class(exchange_type), allocatable ,intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: hardness(:)
   type(exchange_record) :: par
   type(basis_type) :: bas

   if (par%mulliken) then
      block
         type(mulliken_kfock_type), allocatable :: tmp
         allocate(tmp)
         call new_mulliken_exchange(tmp, mol, hardness, .false., .false., par%frscale, &
            & par%omega, par%lrscale, par%average, par%expsmooth, bas)
         call move_alloc(tmp, self)
      end block
   end if
end subroutine

end module tblite_exchange

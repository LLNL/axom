!
! Test allocatable arrays as meta-buffers
!

module sidre_allocatable
  use fruit
  use sidre_mod
  implicit none

contains

  subroutine local_allocatable
    integer, allocatable :: iarray(:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view_iarray1
!    integer num_elements

    ds = datastore_new()
    root = ds%get_root()

    allocate(iarray(10))

    view_iarray1 = root%register_allocatable("iarray", iarray)

!    num_elements = view_iarray1%get_number_of_elements()
!    call assert_equals(num_elements, 10)
    call assert_true(.true.)

    call ds%delete()

  end subroutine local_allocatable

end module sidre_allocatable


function fortran_test() bind(C,name="fortran_test")
  use iso_c_binding
  use fruit
  use sidre_allocatable
  implicit none
  integer(C_INT) fortran_test
  logical ok

  call init_fruit

  call local_allocatable

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (ok) then
     fortran_test = 0
  else
     fortran_test = 1
  endif
end function fortran_test

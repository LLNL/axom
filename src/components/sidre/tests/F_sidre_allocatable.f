!
! Test allocatable arrays as meta-buffers
!

module sidre_allocatable
  use fruit
  use sidre_mod
  implicit none

contains

! Use Fortran to allocate, register with datastore then 
! query metadata using datastore API.
!----------------------------------------------------------------------

  subroutine local_allocatable_int
    integer, allocatable :: iarray(:)
    integer, pointer :: ipointer(:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view
    integer type
    integer num_elements
    integer i

    ds = datastore_new()
    root = ds%get_root()

    allocate(iarray(10))
    do i=1,10
       iarray(i) = i
    enddo

    view = root%create_allocatable_view("iarray", iarray)

!XXX    type = view%get_type_id()
!XXX    call assert_equals(type, ATK_C_INT_T)

    num_elements = view%get_number_of_elements()
    call assert_equals(num_elements, 10)

    ! get array via a pointer
    call view%get_value(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

  end subroutine local_allocatable_int

!----------------------------------------------------------------------
! Allocate array via the datastore
  subroutine ds_allocatable_int
    integer, allocatable :: iarray(:)
    integer, pointer :: ipointer(:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view
    integer num_elements
    integer type
    integer i

    ds = datastore_new()
    root = ds%get_root()

    view = root%create_allocatable_view("iarray", iarray)

!XXX    type = view%get_type_id()
!XXX    call assert_equals(type, ATK_C_INT_T)

    call view%declare(ATK_C_INT_T, 10)
    call view%allocate()
    
    call assert_true(allocated(iarray))

    ! Check size intrinsic
    call assert_equals(size(iarray), 10)

    do i=1,10
       iarray(i) = i
    enddo
    
    ! get array via a pointer
    call view%get_value(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

    ! deleting the datastore deallocates iarray
    call assert_false(allocated(iarray))

  end subroutine ds_allocatable_int

!----------------------------------------------------------------------
!
! register a static array with the datastore

  subroutine local_static_int_array
    integer :: iarray(10)
    integer, pointer :: ipointer(:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view
    integer type
    integer num_elements
    integer i

    ds = datastore_new()
    root = ds%get_root()

    do i=1,10
       iarray(i) = i
    enddo

    view = root%register_static("iarray", iarray)

!XXX    type = view%get_type_id()
!XXX    call assert_equals(type, ATK_C_INT_T)

    num_elements = view%get_number_of_elements()
    call assert_equals(num_elements, 10)

    ! get array via a pointer
    call view%get_value(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

  end subroutine local_static_int_array

!----------------------------------------------------------------------

  subroutine local_allocatable_double
    use iso_c_binding
    real(C_DOUBLE), allocatable :: iarray(:)
    real(C_DOUBLE), pointer :: ipointer(:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view
    integer num_elements
    integer type
    integer i

    ds = datastore_new()
    root = ds%get_root()

    allocate(iarray(10))
    do i=1,10
       iarray(i) = i + 0.5d0
    enddo

    view = root%create_allocatable_view("iarray", iarray)

!XXX    type = view%get_type_id()
!XXX    call assert_equals(type, ATK_C_DOUBLE_T)

    num_elements = view%get_number_of_elements()
    call assert_equals(num_elements, 10)

    ! get array via a pointer
    call view%get_value(ipointer)
    call assert_true(all(abs(iarray-ipointer).lt..0005))

    call ds%delete()

  end subroutine local_allocatable_double

end module sidre_allocatable


function fortran_test() bind(C,name="fortran_test")
  use iso_c_binding
  use fruit
  use sidre_allocatable
  implicit none
  integer(C_INT) fortran_test
  logical ok

  call init_fruit

  call local_allocatable_int
  call ds_allocatable_int
  call local_static_int_array
  call local_allocatable_double

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (ok) then
     fortran_test = 0
  else
     fortran_test = 1
  endif
end function fortran_test

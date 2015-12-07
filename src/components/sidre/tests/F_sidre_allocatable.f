!
! Test allocatable arrays as meta-buffers
!

module sidre_allocatable
  use iso_c_binding
  use fruit
  use sidre_mod
  implicit none

  ! Global variables to add to datastore.
  ! If they were local to a subroutine, then Fortran will attempt to
  ! free them before returning and we want the datastore to free them.
  integer, allocatable :: iarray(:)
  real(C_DOUBLE), allocatable :: darray(:)

contains

! Allocate array via Fortran
! Register with datastore then 
! Query metadata using datastore API.
!----------------------------------------------------------------------

  subroutine local_allocatable_int
    integer, pointer :: ipointer(:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view
    integer type
    integer num_elements
    integer i

    ds = datastore_new()
    root = ds%get_root()

    ! Allocate array via Fortran
    allocate(iarray(10))

    view = root%create_allocatable_view("iarray", iarray)

    type = view%get_type_id()
    call assert_equals(type, SIDRE_INT_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, 10)

    do i=1,10
       iarray(i) = i
    enddo

    ! get array via a pointer
    call view%get_value(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

  end subroutine local_allocatable_int

!----------------------------------------------------------------------

! Register with datastore, check type and length
! Allocate array via the datastore
! Check from Fortran with ALLOCATED and SIZE
! Check datastore metadata
  subroutine ds_allocatable_int
    integer, pointer :: ipointer(:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view
    integer num_elements
    integer type
    integer i

    ds = datastore_new()
    root = ds%get_root()

    ! Register with datastore, check type and length
    view = root%create_allocatable_view("iarray", iarray)

    type = view%get_type_id()
    call assert_equals(type, SIDRE_INT_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, 0)

    ! Allocate array via datastore
! To be consistent with actual Fortran code, the method that creates 
! an allocatable view should take the type and the allocate method
! should take only shape, length, etc.
    call view%allocate(SIDRE_INT_ID, 10)
    
    ! Check from Fortran with ALLOCATED and SIZE
    call assert_true(allocated(iarray))

    ! Check size intrinsic
    call assert_equals(size(iarray), 10)

! Check datastore metadata
    type = view%get_type_id()
    call assert_equals(type, SIDRE_INT_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, 10)

    ! get array via a pointer
    do i=1,10
       iarray(i) = i
    enddo
    
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

    view = root%create_array_view("iarray", iarray)

    type = view%get_type_id()
    call assert_equals(type, SIDRE_INT_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, 10)

    ! get array via a pointer
    call view%get_value(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

  end subroutine local_static_int_array

!----------------------------------------------------------------------
!--- check other types

  subroutine local_allocatable_double
    real(C_DOUBLE), pointer :: dpointer(:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view
    integer num_elements
    integer type
    integer i

    ds = datastore_new()
    root = ds%get_root()

    allocate(darray(10))

    view = root%create_allocatable_view("darray", darray)

    type = view%get_type_id()
    call assert_equals(type, SIDRE_DOUBLE_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, 10)

    do i=1,10
       darray(i) = i + 0.5d0
    enddo

    ! get array via a pointer
    call view%get_value(dpointer)
    call assert_true(all(abs(darray-dpointer).lt..0005))

    call ds%delete()

  end subroutine local_allocatable_double

!----------------------------------------------------------------------

end module sidre_allocatable


program fortran_test
  use fruit
  use sidre_allocatable
  implicit none
  logical ok

  call init_fruit

  call local_allocatable_int
  call ds_allocatable_int
  call local_static_int_array
  call local_allocatable_double

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif
end program fortran_test

!
! Test allocatable arrays as meta-buffers
!

module sidre_allocatable
  use iso_c_binding
  use fruit
  use sidre_mod
  implicit none

contains

! Allocate array via Fortran
! Register with datastore then 
! Query metadata using datastore API.
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
    integer rank
    integer(SIDRE_LENGTH) extents(7)

    ds = datastore_new()
    root = ds%get_root()

    allocate(iarray(10))

    do i=1,10
       iarray(i) = i
    enddo

    view = root%create_array_view("iarray", iarray)

    type = view%get_type_id()
    call assert_equals(type, SIDRE_INT_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, size(iarray))

    rank = view%get_num_dimensions()
    call assert_equals(rank, 1)

    rank = view%get_shape(7, extents)
    call assert_equals(rank, 1)
    call assert_true(extents(1) == size(iarray, 1))

    ! get array via a pointer
    call view%get_value(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

    deallocate(iarray)

  end subroutine local_allocatable_int

!----------------------------------------------------------------------

  subroutine local_allocatable_int_3d
    integer, allocatable :: iarray(:,:,:)
    integer, pointer :: ipointer(:,:,:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view
    integer type
    integer num_elements
    integer i, j, k
    integer rank
    integer(SIDRE_LENGTH) extents(7)

    ds = datastore_new()
    root = ds%get_root()

    allocate(iarray(2,3,4))

    do i=1,2
       do j=1,3
          do k=1,4
             iarray(i,j,k) = i*100 + j*1-0 + k
          enddo
       enddo
    enddo

    view = root%create_array_view("iarray", iarray)

    type = view%get_type_id()
    call assert_equals(type, SIDRE_INT_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, size(iarray))

    rank = view%get_num_dimensions()
    call assert_equals(rank, 3)

    rank = view%get_shape(7, extents)
    call assert_equals(rank, 3)
    call assert_true(extents(1) == size(iarray, 1))
    call assert_true(extents(2) == size(iarray, 2))
    call assert_true(extents(3) == size(iarray, 3))

    ! get array via a pointer
    call view%get_value(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

    deallocate(iarray)

  end subroutine local_allocatable_int_3d

!----------------------------------------------------------------------
!
! register a static (non-allocatable) array with the datastore as external view

  subroutine local_static_int
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

  end subroutine local_static_int

!----------------------------------------------------------------------
!--- check other types

  subroutine local_allocatable_double
    real(C_DOUBLE), allocatable :: darray(:)
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

    do i=1,10
       darray(i) = i + 0.5d0
    enddo

    view = root%create_array_view("darray", darray)

    type = view%get_type_id()
    call assert_equals(type, SIDRE_DOUBLE_ID)

    num_elements = view%get_num_elements()
    call assert_equals(num_elements, 10)

    ! get array via a pointer
    call view%get_value(dpointer)
    call assert_true(all(abs(darray-dpointer).lt..0005))

    call ds%delete()

    deallocate(darray)

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
  call local_allocatable_int_3d
  call local_static_int
  call local_allocatable_double

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif
end program fortran_test

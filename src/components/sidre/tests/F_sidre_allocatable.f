!
! Test allocatable arrays as meta-buffers
!

module sidre_allocatable
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
    type(dataview)  view_iarray1
    integer type
    integer num_elements
    integer i

    ds = datastore_new()
    root = ds%get_root()

    ! Allocate array via Fortran
    allocate(iarray(10))

    ! Register with datastore then 
    view_iarray1 = root%register_allocatable("iarray", iarray)

! Query metadata using datastore API.
!XXX    type = view_iarray1%get_type_id()
!XXX    call assert_equals(type, ATK_C_INT_T)

    num_elements = view_iarray1%get_number_of_elements()
    call assert_equals(num_elements, 10)

    do i=1,10
       iarray(i) = i
    enddo

    ! get array via a pointer
    call view_iarray1%get_value(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

  end subroutine local_allocatable_int

!----------------------------------------------------------------------
! Register with datastore, check type and length
! Allocate array via the datastore
! Check from Fortran with ALLOCATED and SIZE
! Check datastore metadata
  subroutine ds_allocatable_int
    integer, allocatable :: iarray(:)
    integer, pointer :: ipointer(:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view_iarray1
    integer num_elements
    integer type
    integer i

    ds = datastore_new()
    root = ds%get_root()

    ! Register with datastore, check type and length
    view_iarray1 = root%register_allocatable("iarray", iarray)

!XXX    type = view_iarray1%get_type_id()
!XXX    call assert_equals(type, ATK_C_INT_T)

    num_elements = view_iarray1%get_number_of_elements()
    call assert_equals(num_elements, 0)

    ! Allocate array via datastore
    call view_iarray1%declare(ATK_C_INT_T, 10)
    call view_iarray1%allocate()
    
    ! Check from Fortran with ALLOCATED and SIZE
    call assert_true(allocated(iarray))

    ! Check size intrinsic
    call assert_equals(size(iarray), 10)

! Check datastore metadata
!XXX    type = view_iarray1%get_type_id()
!XXX    call assert_equals(type, ATK_C_INT_T)

    num_elements = view_iarray1%get_number_of_elements()
    call assert_equals(num_elements, 10)

    ! get array via a pointer
    do i=1,10
       iarray(i) = i
    enddo
    
    call view_iarray1%get_value(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

  end subroutine ds_allocatable_int


!----------------------------------------------------------------------
! Register with datastore
! Allocate array via Fortran
! Check datastore metadata
  subroutine sync_allocatable_int
    integer, allocatable :: iarray(:)
    integer, pointer :: ipointer(:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view_iarray1
    integer num_elements
    integer type
    integer i

    ds = datastore_new()
    root = ds%get_root()

    ! Register with datastore
    view_iarray1 = root%register_allocatable("iarray", iarray)

!XXX    type = view_iarray1%get_type_id()
!XXX    call assert_equals(type, ATK_C_INT_T)

    ! Allocate array via Fortran
    allocate(iarray(10))
    call assert_true(allocated(iarray))
    call assert_equals(size(iarray), 10)

    ! Check datastore metadata
    call view_iarray1%sync
!XXX    type = view_iarray1%get_type_id()
!XXX    call assert_equals(type, ATK_C_INT_T)

    num_elements = view_iarray1%get_number_of_elements()
    call assert_equals(num_elements, 10)

    do i=1,10
       iarray(i) = i
    enddo

    ! get array via a pointer
    call view_iarray1%get_value(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

  end subroutine sync_allocatable_int

!----------------------------------------------------------------------
!
! register a static array with the datastore

  subroutine local_static_int_array
    integer :: iarray(10)
    integer, pointer :: ipointer(:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view_iarray1
    integer type
    integer num_elements
    integer i

    ds = datastore_new()
    root = ds%get_root()

    do i=1,10
       iarray(i) = i
    enddo

    view_iarray1 = root%register_static("iarray", iarray)

!XXX    type = view_iarray1%get_type_id()
!XXX    call assert_equals(type, ATK_C_INT_T)

    num_elements = view_iarray1%get_number_of_elements()
    call assert_equals(num_elements, 10)

    ! get array via a pointer
    call view_iarray1%get_value(ipointer)
    call assert_true(all(iarray.eq.ipointer))

    call ds%delete()

  end subroutine local_static_int_array

!--- check other types

  subroutine local_allocatable_double
    use iso_c_binding
    real(C_DOUBLE), allocatable :: iarray(:)
    real(C_DOUBLE), pointer :: ipointer(:)

    type(datastore) ds
    type(datagroup) root
    type(dataview)  view_iarray1
    integer num_elements
    integer type
    integer i

    ds = datastore_new()
    root = ds%get_root()

    allocate(iarray(10))
    do i=1,10
       iarray(i) = i + 0.5d0
    enddo

    view_iarray1 = root%register_allocatable("iarray", iarray)

!XXX    type = view_iarray1%get_type_id()
!XXX    call assert_equals(type, ATK_C_DOUBLE_T)

    num_elements = view_iarray1%get_number_of_elements()
    call assert_equals(num_elements, 10)

    ! get array via a pointer
    call view_iarray1%get_value(ipointer)
    call assert_true(all(abs(iarray-ipointer).lt..0005))

    call ds%delete()

  end subroutine local_allocatable_double

!----------------------------------------------------------------------

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
!  call sync_allocatable_int
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

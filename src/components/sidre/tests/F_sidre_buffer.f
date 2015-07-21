!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!

module sidre_buffer
  use iso_c_binding
  use fruit
  use sidre_mod
  implicit none

contains
!------------------------------------------------------------------------------

  subroutine create_buffers
    type(datastore) ds
    type(databuffer) dbuff_0, dbuff_1, dbuff_3

    ds = datastore_new()

    dbuff_0 = ds%create_buffer()
    dbuff_1 = ds%create_buffer()

    call assert_equals(dbuff_0%get_index(), 0)
    call assert_equals(dbuff_1%get_index(), 1)
    call ds%destroy_buffer(0)

    dbuff_3 = ds%create_buffer()
    call assert_equals(dbuff_3%get_index(), 0)

    call ds%print()
    call ds%delete()
  end subroutine create_buffers

!------------------------------------------------------------------------------

  subroutine alloc_buffer_for_int_array
    type(datastore) ds
    type(databuffer) dbuff
    type(C_PTR) data_ptr
    integer(C_INT), pointer :: data(:)
    integer i
    
    ds = datastore_new()
    dbuff = ds%create_buffer()

    call dbuff%declare(ATK_C_INT_T, 10_8)
    call dbuff%allocate()

    data_ptr = dbuff%get_data()
    call c_f_pointer(data_ptr, data, [ 10 ])

    do i = 1, 10
       data(i) = i * i
    enddo

    call dbuff%print()

!    call assert_equals(dbuff%getNode().schema().total_bytes(), &
!         dbuff%getSchema().total_bytes())

    call ds%print()
    call ds%delete()
  end subroutine alloc_buffer_for_int_array

!------------------------------------------------------------------------------

  subroutine init_buffer_for_int_array
    type(datastore) ds
    type(databuffer) dbuff
    type(C_PTR) data_ptr
    integer(C_INT), pointer :: data(:)
    integer i

    ds = datastore_new()
    dbuff = ds%create_buffer()

    call dbuff%allocate(ATK_C_INT_T, 10_8)
    data_ptr = dbuff%get_data()
    call c_f_pointer(data_ptr, data, [ 10 ])

    do i = 1, 10
       data(i) = i * i
    enddo

    call dbuff%print()

!  call assert_equals(dbuff%getNode().schema().total_bytes(),
!            dbuff%getSchema().total_bytes())

    call ds%print()
    call ds%delete()
  end subroutine init_buffer_for_int_array

!------------------------------------------------------------------------------

  subroutine realloc_buffer
    type(datastore) ds
    type(databuffer) dbuff
    type(C_PTR) data_ptr
    integer(C_LONG), pointer :: data(:)
    integer i

    ds = datastore_new()

    dbuff = ds%create_buffer()

    call dbuff%allocate(ATK_C_LONG_T, 5_8)

!    call assert_equals(dbuff%getNode().schema().total_bytes(), sizeof(long)*5)

    data_ptr = dbuff%get_data()
    call c_f_pointer(data_ptr, data, [ 5 ])

    data(:) = 5

    call dbuff%print()
  
    call dbuff%reallocate(ATK_C_LONG_T, 10_8)

    ! data buffer changes
    data_ptr = dbuff%get_data()
    call c_f_pointer(data_ptr, data, [ 10 ])

    do i = 1, 5
       call assert_equals(int(data(i)), 5)  ! XXX cast
    enddo

    do i = 6, 10
       data(i) = 10
    enddo

!  call assert_equals(dbuff%getNode().schema().total_bytes(), sizeof(long)*10)

    call dbuff%print()

    call ds%print()
    call ds%delete()
  end subroutine realloc_buffer

!----------------------------------------------------------------------
end module sidre_buffer
!----------------------------------------------------------------------

function fortran_test() bind(C,name="fortran_test")
  use fruit
  use sidre_buffer
  implicit none
  integer(C_INT) fortran_test
  logical ok

  call init_fruit

  call create_buffers
  call alloc_buffer_for_int_array
  call init_buffer_for_int_array
  call realloc_buffer

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (ok) then
     fortran_test = 0
  else
     fortran_test = 1
  endif
end function fortran_test

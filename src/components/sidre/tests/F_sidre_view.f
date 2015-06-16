!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!

module sidre_view
  use iso_c_binding
  use fruit
  use sidre_mod
  implicit none

contains
!------------------------------------------------------------------------------

  subroutine create_views()
    type(datastore) ds
    type(datagroup) root
    type(dataview) dv_0, dv_1
    type(databuffer) db_0, db_1

    ds = datastore_new()
    root = ds%get_root()

    dv_0 = root%create_view_and_buffer("field0")
    dv_1 = root%create_view_and_buffer("field1")

    db_0 = dv_0%get_buffer()
    db_1 = dv_1%get_buffer()

    call assert_equals(db_0%get_index(), 0)
    call assert_equals(db_1%get_index(), 1)
    call datastore_delete(ds)
  end subroutine create_views

!------------------------------------------------------------------------------

  subroutine int_buffer_from_view()
    type(datastore) ds
    type(datagroup) root
    type(dataview) dv
    type(C_PTR) data_ptr
    integer(C_INT), pointer :: data(:)
    integer i

    ds = datastore_new()
    root = ds%get_root()

    dv = root%create_view_and_buffer("u0")
    call assert_equals(dv%get_type_id(), ATK_INT32_T)  ! XXX NATIVE TYPE
    call dv%allocate(ATK_C_INT_T, 10_8)
    data_ptr = dv%get_data_in_buffer()
    call c_f_pointer(data_ptr, data, [ 10 ])

    do i = 1, 10
       data(i) = i * i
    enddo

!--#ifdef XXX
!--    dv->getNode().print_detailed()
!--#endif
!--
!--    !  EXPECT_EQ(ATK_dataview_get_total_bytes(dv), dv->getSchema().total_bytes())
!--    call assert_equals(dv%get_total_bytes(), sizeof(int) * 10)
    call datastore_delete(ds)
  end subroutine int_buffer_from_view

!------------------------------------------------------------------------------

  subroutine int_buffer_from_view_conduit_value()
    type(datastore) ds
    type(datagroup) root
    type(dataview) dv
    type(C_PTR) data_ptr
    integer(C_INT), pointer :: data(:)
    integer i

    ds = datastore_new()
    root = ds%get_root()

    dv = root%create_view_and_buffer("u0", ATK_C_INT_T, 10_8)
    data_ptr = dv%get_data_in_buffer()
    call c_f_pointer(data_ptr, data, [ 10 ])

    do i = 1, 10
       data(i) = i * i
    enddo

!--#ifdef XXX
!--    dv->getNode().print_detailed()
!--#endif

!--    EXPECT_EQ(ATK_dataview_get_total_bytes(dv), sizeof(int) * 10)
    call datastore_delete(ds)
  end subroutine int_buffer_from_view_conduit_value

!------------------------------------------------------------------------------

  subroutine int_array_multi_view()
    type(datastore) ds
    type(datagroup) root
    type(databuffer) dbuff
    type(dataview) dv_e, dv_o 
    type(C_PTR) data_ptr
    integer(C_INT), pointer :: data(:)
    integer i

    ds = datastore_new()
    root = ds%get_root()
    dbuff = ds%create_buffer()

    call dbuff%declare(ATK_C_INT_T, 10_8)
    call dbuff%allocate()
    data_ptr = dbuff%get_data()
    call c_f_pointer(data_ptr, data, [ 10 ])

    do i = 1, 10
       data(i) = i
    enddo

!--#ifdef XXX
!--    dbuff->getNode().print_detailed()
!--
!--    EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
!--              dbuff->getSchema().total_bytes())
!--#endif

    dv_e = root%create_view("even", dbuff)
    dv_o = root%create_view("odd", dbuff)

!--#ifdef XXX
!--  dv_e->apply(DataType::uint32(5,0,8))
!--
!--  dv_o->apply(DataType::uint32(5,4,8))
!--
!--  dv_e->getNode().print_detailed()
!--  dv_o->getNode().print_detailed()
!--
!--  uint32_array dv_e_ptr = dv_e->getNode().as_uint32_array()
!--  uint32_array dv_o_ptr = dv_o->getNode().as_uint32_array()
!--  for(int i=0  i<5  i++)
!--  {
!--    std::cout << "idx:" <<  i
!--              << " e:" << dv_e_ptr[i]
!--              << " o:" << dv_o_ptr[i]
!--              << " em:" << dv_e_ptr[i]  % 2
!--              << " om:" << dv_o_ptr[i]  % 2
!--              << std::endl
!--
!--    EXPECT_EQ(dv_e_ptr[i] % 2, 0u)
!--    EXPECT_EQ(dv_o_ptr[i] % 2, 1u)
!--  }
!--#endif
    call ds%print()
    call datastore_delete(ds)
  end subroutine int_array_multi_view

!------------------------------------------------------------------------------

  subroutine init_int_array_multi_view()
    type(datastore) ds
    type(datagroup) root
    type(databuffer) dbuff
    type(dataview) dv_e, dv_o 
    type(C_PTR) data_ptr
    integer, pointer :: data(:)
    integer i
    
    ds = datastore_new()
    root = ds%get_root()
    dbuff = ds%create_buffer()
    
    call dbuff%allocate(ATK_C_INT_T, 10_8)
    data_ptr = dbuff%get_data()
    call c_f_pointer(data_ptr, data, [ 10 ])

    do i = 1, 10
       data(i) = i
    enddo

!--#ifdef XXX
!--  dbuff->getNode().print_detailed()
!--
!--  EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
!--            dbuff->getSchema().total_bytes())
!--#endif

    dv_e = root%create_view("even", dbuff)
    dv_o = root%create_view("odd", dbuff)

!--#ifdef XXX
!--  ! uint32(num_elems, offset, stride)
!--  dv_e->apply(DataType::uint32(5,0,8))
!--
!--
!--  ! uint32(num_elems, offset, stride)
!--  dv_o->apply(DataType::uint32(5,4,8))
!--
!--
!--  dv_e->getNode().print_detailed()
!--  dv_o->getNode().print_detailed()
!--
!--  uint32_array dv_e_ptr = dv_e->getNode().as_uint32_array()
!--  uint32_array dv_o_ptr = dv_o->getNode().as_uint32_array()
!--  for(int i=0  i<5  i++)
!--  {
!--    std::cout << "idx:" <<  i
!--              << " e:" << dv_e_ptr[i]
!--              << " o:" << dv_o_ptr[i]
!--              << " em:" << dv_e_ptr[i]  % 2
!--              << " om:" << dv_o_ptr[i]  % 2
!--              << std::endl
!--
!--    EXPECT_EQ(dv_e_ptr[i] % 2, 0u)
!--    EXPECT_EQ(dv_o_ptr[i] % 2, 1u)
!--  }
!--#endif

    call ds%print()
    call datastore_delete(ds)
  end subroutine init_int_array_multi_view

!------------------------------------------------------------------------------

  subroutine int_array_multi_view_resize()
     !
     ! This example creates a 4 * 10 buffer of ints,
     ! and 4 views that point the 4 sections of 10 ints
     !
     ! We then create a new buffer to support 4*12 ints
     ! and 4 views that point into them
     !
     ! after this we use the old buffers to copy the values
     ! into the new views
     !
    type(datastore) ds
    type(datagroup) root, r_old
    type(dataview) base_old
    type(C_PTR) data_ptr
    integer(C_INT), pointer :: data(:)
    integer i

    ! create our main data store
    ds = datastore_new()

    ! get access to our root data Group
    root = ds%get_root()

    ! create a group to hold the "old" or data we want to copy
    r_old = root%create_group("r_old")
    ! create a view to hold the base buffer
    base_old = r_old%create_view_and_buffer("base_data")

    ! alloc our buffer
    ! we will create 4 sub views of this array
    call base_old%allocate(ATK_C_INT_T, 40_8)
    data_ptr = base_old%get_data_in_buffer()
    call c_f_pointer(data_ptr, data, [ 40 ])

    ! init the buff with values that align with the
    ! 4 subsections.
    do i = 1, 10
       data(i) = 1
    enddo
    do i = 11, 20
       data(i) = 2
    enddo
    do i = 21, 30
       data(i) = 3
    enddo
    do i = 31, 40
       data(i) = 4
    enddo

!--#ifdef XXX
!--  ! setup our 4 views
!--  ATK_databuffer * buff_old = ATK_dataview_get_buffer(base_old)
!--  buff_old->getNode().print()
!--  ATK_dataview * r0_old = ATK_dataview_create_view(r_old, "r0",buff_old)
!--  ATK_dataview * r1_old = ATK_dataview_create_view(r_old, "r1",buff_old)
!--  ATK_dataview * r2_old = ATK_dataview_create_view(r_old, "r2",buff_old)
!--  ATK_dataview * r3_old = ATK_dataview_create_view(r_old, "r3",buff_old)
!--
!--  ! each view is offset by 10 * the # of bytes in a uint32
!--  ! uint32(num_elems, offset)
!--  index_t offset =0
!--  r0_old->apply(r0_old, DataType::uint32(10,offset))
!--
!--  offset += sizeof(int) * 10
!--  r1_old->apply(r1_old, DataType::uint32(10,offset))
!--
!--  offset += sizeof(int) * 10
!--  r2_old->apply(r2_old, DataType::uint32(10,offset))
!--
!--  offset += sizeof(int) * 10
!--  r3_old->apply(r3_old, DataType::uint32(10,offset))
!--
!--  ! check that our views actually point to the expected data
!--  !
!--  uint32 * r0_ptr = r0_old->getNode().as_uint32_ptr()
!--  for(int i=0  i<10  i++)
!--  {
!--    EXPECT_EQ(r0_ptr[i], 1u)
!--    ! check pointer relation
!--    EXPECT_EQ(&r0_ptr[i], &data_ptr[i])
!--  }
!--
!--  uint32 * r3_ptr = r3_old->getNode().as_uint32_ptr()
!--  for(int i=0  i<10  i++)
!--  {
!--    EXPECT_EQ(r3_ptr[i], 4u)
!--    ! check pointer relation
!--    EXPECT_EQ(&r3_ptr[i], &data_ptr[i+30])
!--  }
!--
!--  ! create a group to hold the "old" or data we want to copy into
!--  ATK_datagroup * r_new = ATK_datagroup_create_group(root, "r_new")
!--  ! create a view to hold the base buffer
!--  ATK_dataview * base_new = ATK_datagroup_create_view_and_buffer(r_new, "base_data")
!--
!--  ! alloc our buffer
!--  ! create a buffer to hold larger subarrays
!--  base_new->allocate(base_new, DataType::uint32(4 * 12))
!--  int* base_new_data = (int *) ATK_databuffer_det_data(base_new)
!--  for (int i = 0 i < 4 * 12 ++i) 
!--  {
!--     base_new_data[i] = 0
!--  } 
!--
!--  ATK_databuffer * buff_new = ATK_dataview_get_buffer(base_new)
!--  buff_new->getNode().print()
!--
!--  ! create the 4 sub views of this array
!--  ATK_dataview * r0_new = ATK_datagroup_create_view(r_new, "r0",buff_new)
!--  ATK_dataview * r1_new = ATK_datagroup_create_view(r_new, "r1",buff_new)
!--  ATK_dataview * r2_new = ATK_datagroup_create_view(r_new, "r2",buff_new)
!--  ATK_dataview * r3_new = ATK_datagroup_create_view(r_new, "r3",buff_new)
!--
!--  ! apply views to r0,r1,r2,r3
!--  ! each view is offset by 12 * the # of bytes in a uint32
!--
!--  ! uint32(num_elems, offset)
!--  offset =0
!--  r0_new->apply(DataType::uint32(12,offset))
!--
!--  offset += sizeof(int) * 12
!--  r1_new->apply(DataType::uint32(12,offset))
!--
!--  offset += sizeof(int) * 12
!--  r2_new->apply(DataType::uint32(12,offset))
!--
!--  offset += sizeof(int) * 12
!--  r3_new->apply(DataType::uint32(12,offset))
!--
!--  ! update r2 as an example first
!--  buff_new->getNode().print()
!--  r2_new->getNode().print()
!--
!--  ! copy the subset of value
!--  r2_new->getNode().update(r2_old->getNode())
!--  r2_new->getNode().print()
!--  buff_new->getNode().print()
!--
!--
!--  ! check pointer values
!--  int * r2_new_ptr = (int *) ATK_dataview_get_data_in_buffer(r2_new)
!--
!--  for(int i=0  i<10  i++)
!--  {
!--    EXPECT_EQ(r2_new_ptr[i], 3)
!--  }
!--
!--  for(int i=10  i<12  i++)
!--  {
!--    EXPECT_EQ(r2_new_ptr[i], 0)     ! assumes zero-ed alloc
!--  }
!--
!--
!--  ! update the other views
!--  r0_new->getNode().update(r0_old->getNode())
!--  r1_new->getNode().update(r1_old->getNode())
!--  r3_new->getNode().update(r3_old->getNode())
!--
!--  buff_new->getNode().print()
!--#endif

    call ds%print()
    call datastore_delete(ds)

  end subroutine int_array_multi_view_resize

!------------------------------------------------------------------------------

  subroutine int_array_realloc()
    !
    ! info
    !
    type(datastore) ds
    type(datagroup) root
    type(dataview) a1, a2
    type(C_PTR) a1_ptr, a2_ptr
    real(C_FLOAT), pointer :: a1_data(:)
    integer(C_INT), pointer :: a2_data(:)
    integer i

    ! create our main data store
    ds = datastore_new()

    ! get access to our root data Group
    root = ds%get_root()

    ! create a view to hold the base buffer
    a1 = root%create_view_and_buffer("a1", ATK_C_FLOAT_T, 5_8)
    a2 = root%create_view_and_buffer("a2", ATK_C_INT_T, 5_8)

    a1_ptr = a1%get_data_in_buffer()
    a2_ptr = a2%get_data_in_buffer()
    call c_f_pointer(a1_ptr, a1_data, [ 5 ])
    call c_f_pointer(a2_ptr, a2_data, [ 5 ])

    do i = 1, 5
       a1_data(i) =  5.0
       a2_data(i) = -5
    enddo

!--  EXPECT_EQ(ATK_dataview_get_total_bytes(a1), sizeof(float)*5)
!--  EXPECT_EQ(ATK_dataview_get_total_bytes(a2), sizeof(int)*5)


    call a1%reallocate(ATK_C_FLOAT_T, 10_8)
    call a2%reallocate(ATK_C_INT_T, 15_8)

    a1_ptr = a1%get_data_in_buffer()
    a2_ptr = a2%get_data_in_buffer()
    call c_f_pointer(a1_ptr, a1_data, [ 10 ])
    call c_f_pointer(a2_ptr, a2_data, [ 15 ])

    do i = 1, 5
       call assert_equals(a1_data(i), 5.0)
       call assert_equals(a2_data(i), -5)
    enddo

    a1_data(6:10) = 10.0
    a2_data(6:10) = -10

    a2_data(11:15) = -15

!--  EXPECT_EQ(ATK_dataview_get_total_bytes(a1), sizeof(float)*10)
!--  EXPECT_EQ(ATK_dataview_get_total_bytes(a2), sizeof(int)*15)

    call ds%print()
    call datastore_delete(ds)

  end subroutine int_array_realloc

!------------------------------------------------------------------------------

  subroutine simple_opaque()
    type(datastore) ds
    type(datagroup) root
    type(dataview) opq_view
    integer(C_INT), target :: src_data
    integer(C_INT), pointer :: out_data
    type(C_PTR) src_ptr, opq_ptr

    ! create our main data store
    ds = datastore_new()

    ! get access to our root data Group
    root = ds%get_root()

    src_data = 42
   
    src_ptr = c_loc(src_data)

    opq_view = root%create_opaque_view("my_opaque", src_ptr)

    ! we shouldn't have any buffers
!XX    call assert_equals(ds%get_num_buffers(), 0)

    call assert_true(opq_view%is_opaque())

    opq_ptr = opq_view%get_opaque()
    call c_f_pointer(opq_ptr, out_data)

!XX    call assert_equals(opq_ptr, src_ptr)
    call assert_equals(out_data, 42)

    call ds%print()
    call datastore_delete(ds)
!--  free(src_data)
  end subroutine simple_opaque

!----------------------------------------------------------------------
end module sidre_view
!----------------------------------------------------------------------

program tester
  use fruit
  use sidre_view
  call init_fruit

  call create_views
  call int_buffer_from_view
  call int_buffer_from_view_conduit_value
!  call int_array_multi_view
!  call init_int_array_multi_view
  call int_array_multi_view_resize
  call int_array_realloc
  call simple_opaque

  call fruit_summary
  call fruit_finalize
end program tester

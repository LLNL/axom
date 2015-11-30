!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!

! API coverage tests
! Each test should be documented with the interface functions being tested

module sidre_group
  use iso_c_binding
  use fruit
  use sidre_mod
  implicit none

contains

  !------------------------------------------------------------------------------
  ! get_name()
  !------------------------------------------------------------------------------
  subroutine get_name
    type(datastore) ds
    type(datagroup) root, group
    character(30) name

    ds = datastore_new()
    root = ds%get_root()
    group = root%create_group("test")

    call group%get_name(name)
    call assert_true(name == "test" )
    
    call ds%delete()
  end subroutine get_name

  !------------------------------------------------------------------------------
  ! get_parent()
  !------------------------------------------------------------------------------
  subroutine get_parent
    type(datastore) ds
    type(datagroup) root, parent, child

    ds = datastore_new()
    root = ds%get_root()
    parent = root%create_group("parent")
    child = parent%create_group("child")

    call assert_true( child%get_parent() == parent )

    call ds%delete()
  end subroutine get_parent

!------------------------------------------------------------------------------
! Verify get_data_store()
!------------------------------------------------------------------------------
  subroutine get_datastore
    type(datastore) ds, const_ds
    type(datagroup) root, group

    ds = datastore_new()
    root = ds%get_root()
    group = root%create_group("parent")

    call assert_true( group%get_data_store() == ds )

    const_ds = group%get_data_store()
    call assert_true( const_ds == ds )

    call ds%delete()
  end subroutine get_datastore

  !------------------------------------------------------------------------------
  ! Verify has_group()
  !------------------------------------------------------------------------------
  subroutine has_group
    type(datastore) ds
    type(datagroup) root, parent, child

    ds = datastore_new()
    root = ds%get_root()

    parent = root%create_group("parent")
    child = parent%create_group("child")
    call assert_true( child%get_parent() == parent )

    call assert_true( parent%has_group("child") )

    call ds%delete()
  end subroutine has_group

!------------------------------------------------------------------------------
! Verify has_view()
!------------------------------------------------------------------------------
  subroutine has_view
    type(datastore) ds
    type(datagroup) root,parent
    type(dataview) view

    ds = datastore_new()
    root = ds%get_root()

    parent = root%create_group("parent")
    view = parent%create_view_empty("view")

    call assert_true( view%get_owning_group() == parent )

    call assert_true( parent%has_view("view") )

    call ds%delete()
  end subroutine has_view

!------------------------------------------------------------------------------
! Verify get_view_name(), get_view_index()
!------------------------------------------------------------------------------
  subroutine get_view_name_index
    type(datastore) ds
    type(datagroup) root, parent
    type(dataview) view1, view2
    integer idx1, idx2, idx3    ! IndexType
    character(len=30) name1, name2, name3
    character(len=30) tmpname

    ds = datastore_new()
    root = ds%get_root()

    parent = root%create_group("parent")
    view1 = parent%create_view_empty("view1")
    view2 = parent%create_view_empty("view2")

    call assert_true(parent%get_num_views() == 2)

    idx1 = parent%get_view_index("view1")
    idx2 = parent%get_view_index("view2")

    call parent%get_view_name(idx1, name1)
    call parent%get_view_name(idx2, name2)

    call assert_equals(name1, "view1")
    call view1%get_name(tmpname)
    call assert_equals(tmpname, name1)
    
    call assert_equals(name2, "view2")
    call view2%get_name(tmpname)
    call assert_equals(tmpname, name2)

    idx3 = parent%get_view_index("view3")
    call assert_equals(idx3, invalid_index)

    call parent%get_view_name(idx3, name3)
    call assert_true(name3 == " ")
    call assert_false(name_is_valid(name3))

    call ds%delete()
  end subroutine get_view_name_index

!------------------------------------------------------------------------------
! Verify get_group_name(), get_group_index()
!------------------------------------------------------------------------------
  subroutine get_group_name_index
    type(datastore) ds
    type(datagroup) root, parent, group1, group2
    integer idx1, idx2, idx3     ! IndexType
    character(len=30) name1, name2, name3
    character(len=30) tmpname

    ds = datastore_new()
    root = ds%get_root()

    parent = root%create_group("parent")
    group1 = parent%create_group("group1")
    group2 = parent%create_group("group2")

!--    call assert_equals(parent%get_num_groups(), 2)  ! size_t
    call assert_true(parent%get_num_groups() == 2)

    idx1 = parent%get_group_index("group1")
    idx2 = parent%get_group_index("group2")

    call parent%get_group_name(idx1, name1)
    call assert_equals(name1, "group1")
    call group1%get_name(tmpname)
    call assert_equals(tmpname, name1)

    call parent%get_group_name(idx2, name2)
    call assert_equals(name2, "group2")
    call group2%get_name(tmpname)
    call assert_equals(tmpname, name2)

    idx3 = parent%get_group_index("group3")
    call assert_equals(idx3, invalid_index)

    call parent%get_group_name(idx3, name3)
    call assert_true(name3 == " ")
    call assert_false(name_is_valid(name3))

    call ds%delete()
  end subroutine get_group_name_index

  !------------------------------------------------------------------------------
  ! create_view_empty()
  ! destroy_view()
  ! has_view()
  !------------------------------------------------------------------------------
  subroutine create_destroy_has_view
    type(datastore) ds
    type(datagroup) root,group
    type(dataview) view

    ds = datastore_new()
    root = ds%get_root()
    group = root%create_group("parent")

    view = group%create_view_empty("view")
    call assert_true( group%get_parent() == root )
    call assert_false( view%has_buffer() )

    call assert_true( group%has_view("view") )

    call group%destroy_view("view")

    call assert_false( group%has_view("view") )

    call ds%delete()
  end subroutine create_destroy_has_view

  !------------------------------------------------------------------------------
  ! create_group()
  ! destroy_group()
  ! has_group()
  !------------------------------------------------------------------------------
  subroutine create_destroy_has_group
    type(datastore) ds
    type(datagroup) root, group

    ds = datastore_new()
    root = ds%get_root()
    group = root%create_group("group")
    call assert_true( group%get_parent() == root )

    call assert_true( root%has_group("group") )

    call root%destroy_group("group")
    call assert_false( root%has_group("group") )

    call ds%delete()
  end subroutine create_destroy_has_group

  !------------------------------------------------------------------------------
  subroutine group_name_collisions
    type(datastore) ds
    type(datagroup) root, flds
    type(dataview) view

    ds = datastore_new()
    root = ds%get_root()
    flds = root%create_group("fields")
    view = flds%create_view_empty("a")

    call assert_true(flds%has_view("a"))

    call ds%delete()
  end subroutine group_name_collisions

  !------------------------------------------------------------------------------
  subroutine view_copy_move
    type(datastore) ds
    type(datagroup) root, flds, subgrp
    type(dataview) i0_view, f0_view, d0_view, tmpview

    ds = datastore_new()
    root = ds%get_root()
    flds = root%create_group("fields")

    i0_view = flds%create_view_and_allocate("i0", ATK_C_INT_T, 1)
    call i0_view%set_value(1)

    f0_view = flds%create_view_and_allocate("f0", ATK_C_FLOAT_T, 1)
    call f0_view%set_value(100.0)

    d0_view = flds%create_view_and_allocate("d0", ATK_C_DOUBLE_T, 1)
    call d0_view%set_value(3000.0d0)  ! XXX without d0, error in get_value_double

    call assert_true(flds%has_view("i0"))
    call assert_true(flds%has_view("f0"))
    call assert_true(flds%has_view("d0"))

    ! test moving a view from flds to sub
    subgrp = flds%create_group("sub")
    tmpview = subgrp%move_view(flds%get_view("d0"))
    call flds%print()
    call assert_false(flds%has_view("d0"))
    call assert_true(flds%has_group("sub"))
    call assert_true(subgrp%has_view("d0"))

    ! check the data value
    call assert_equals(tmpview%get_value_double(), 3000.0_C_DOUBLE)

    ! test copying a view from flds to sub
    tmpview = subgrp%copy_view(flds%get_view("i0"))

    call flds%print()

    call assert_true(flds%has_view("i0"))
    call assert_true(subgrp%has_view("i0"))

    ! we expect the actual data  pointers to be the same
!--    call assert_equals(flds%get_view("i0")%getDataBuffer(),
!--    call flds%get_group("sub")%get_view("i0")%getDataBuffer())
    
    call ds%delete()
  end subroutine view_copy_move

  !------------------------------------------------------------------------------
  subroutine groups_move_copy
    type(datastore) ds
    type(datagroup) root, flds, ga, gb, gc
    type(datagroup) subgrp, tmpgrp
    type(dataview) i0_view, f0_view, d0_view, tmpview

    ds = datastore_new()
    root = ds%get_root()
    flds = root%create_group("fields")

    ga = flds%create_group("a")
    gb = flds%create_group("b")
    gc = flds%create_group("c")

    i0_view = ga%create_view_and_allocate("i0", ATK_C_INT_T, 1_8)
    call i0_view%set_value(1)

    f0_view = gb%create_view_and_allocate("f0", ATK_C_FLOAT_T, 1_8)
    call f0_view%set_value(100.0)

    d0_view = gc%create_view_and_allocate("d0", ATK_C_DOUBLE_T, 1_8)
    call d0_view%set_value(3000.0d0)

    ! check that all sub groups exist
    call assert_true(flds%has_group("a"))
    call assert_true(flds%has_group("b"))
    call assert_true(flds%has_group("c"))

    !move "b" to a child of "sub"
    subgrp = flds%create_group("sub")
    tmpgrp = subgrp%move_group(gb)

    call flds%print()

    call assert_true(flds%has_group("a"))
    call assert_true(flds%has_group("sub"))
    call assert_true(flds%has_group("c"))
    
    call assert_true(tmpgrp == gb)
    
    call ds%delete()
  end subroutine groups_move_copy

  !------------------------------------------------------------------------------
  subroutine create_destroy_view_and_data
    type(datastore) ds
    type(datagroup) root, grp
    type(dataview) view1, view2
    type(databuffer) tmpbuf
!XX    type(databuffer) buffer1, tmpbuf
    integer bufferid1      ! IndexType
    character(len=30) view_name1, view_name2
!XX    logical buffvalid

    ds = datastore_new()
    root = ds%get_root()
    grp = root%create_group("grp")

    view_name1 = "viewBuffer1"
    view_name2 = "viewBuffer2"

    view1 = grp%create_view_and_allocate(view_name1, ATK_C_INT_T, 1_8)
    view2 = grp%create_view_and_allocate(view_name2, ATK_C_INT_T, 1_8)

    call assert_true(grp%has_view(view_name1))
    call assert_true(grp%get_view(view_name1) == view1)

    call assert_true(grp%has_view(view_name2))
    call assert_true(grp%get_view(view_name2) == view2)

    tmpbuf = view1%get_buffer()
    bufferid1 = tmpbuf%get_index()

    call grp%destroy_view_and_data(view_name1)


    call assert_false(grp%has_view(view_name1))
    call assert_true(ds%get_num_buffers() == 1)

!XX    buffer1 = ds%get_buffer(bufferId1)
!XX   buffvalid = .true.
!--    if( buffer1 == ATK_NULLPTR ) buffValid = .false.
!XX    call assert_false(buffValid)

    call ds%delete()
  end subroutine create_destroy_view_and_data


  !------------------------------------------------------------------------------
  subroutine create_destroy_alloc_view_and_data
    type(datastore) ds
    type(datagroup) root, grp
    type(dataview) view1
    integer i
    character(len=30) view_name1, view_name2
    integer(C_INT), pointer :: v1_vals(:)
!--    real(C_DOUBLE), pointer :: v2_vals(:)

    ds = datastore_new()
    root = ds%get_root()
    grp = root%create_group("grp")

    view_name1 = "viewBuffer1"
    view_name2 = "viewBuffer2"

    ! use create + alloc convenience methods
    ! this one is the DataType & method
    view1 = grp%create_view_and_allocate(view_name1, ATK_C_INT_T, 10)

!--    ! this one is the Schema & method
!--    Schema s
!--    s.set(DataType::c_double(10))
!--    DataView * const view2 = grp%create_view_and_buffer(view_name2,

    call assert_true(grp%has_view(view_name1))
    call assert_true(grp%get_view(view_name1) == view1)

!--    call assert_true(grp%has_view(view_name2))
!--    call assert_equals( grp%get_view(view_name2), view2 )


    call view1%get_value(v1_vals)
!--    double * v2_vals = view2%get_value()
  
    do i = 1, 10
       v1_vals(i) = i
!--       v2_vals(i) = i * 3.1415
    enddo

    call assert_true(view1%get_num_elements() == 10)
!--    call assert_equals(view2%get_num_elements(), 10)
!--    call assert_equals(view1%get_total_bytes(), 10 * sizeof(int))
!--    call assert_equals(view2%get_total_bytes(), 10 * sizeof(double))

    call grp%destroy_view_and_data(view_name1)
!--    call grp%destroy_view_and_buffer(view_name2)

    call ds%delete()
  end subroutine create_destroy_alloc_view_and_data

  !------------------------------------------------------------------------------
  subroutine create_view_of_buffer_with_schema
    type(datastore) ds
    type(datagroup) root
    type(dataview) base
    type(databuffer) base_buff
    integer(C_INT), pointer :: base_vals(:)
!--    integer i
!    integer(C_INT), pointer :: sub_a_vals(:)

    ds = datastore_new()
    root = ds%get_root()

    ! use create + alloc convenience methods
    ! this one is the DataType & method
    base =  root%create_view_and_allocate("base", ATK_C_INT_T, 10)
    call base%get_value(base_vals)

    base_vals(1:5) = 10
    base_vals(6:10) = 20

    base_buff = base%get_buffer()
    ! create two views into this buffer
    ! view for the first 5 values
!--    root%createView("sub_a", base_buff, DataType::c_int(5))
    ! view for the second 5 values
    !  (schema call path case)
!--    Schema s(DataType::c_int(5,5*sizeof(int)))
!--    root%createView("sub_b",base_buff,s)

!--    int * sub_a_vals = root%get_view("sub_a")%get_value(sub_a_vals)
!--    int * sub_b_vals = root%get_view("sub_b")%get_value(sub_b_vals)

!--    do i = 1, 5
!--       call assert_equals(sub_a_vals(i), 10)
!--       call assert_equals(sub_b_vals(i), 20)
!--    enddo

    call ds%delete()
  end subroutine create_view_of_buffer_with_schema

  !------------------------------------------------------------------------------
  subroutine save_restore_simple
    type(datastore) ds, ds2
    type(datagroup) root, root2, flds, ga
    type(dataview) i0_view

    ds = datastore_new()
    root = ds%get_root()
    flds = root%create_group("fields")

    ga = flds%create_group("a")

    i0_view = ga%create_view_and_allocate("i0", ATK_C_INT_T, 1)
    call i0_view%set_value(1)

    call assert_true(root%has_group("fields"))
    call assert_true(flds%has_group("a"))
    call assert_true(ga%has_view("i0"))

    call root%save("F_out_sidre_group_save_restore_simple","conduit")

    call ds%print()

    ds2 = datastore_new()
    root2 = ds2%get_root()

    call root2%load("F_out_sidre_group_save_restore_simple","conduit")

    call ds2%print()

    flds = root2%get_group("fields")
    ! check that all sub groups exist
    call assert_true(flds%has_group("a"))
    ga = flds%get_group("a")
    i0_view = ga%get_view("i0")
    call assert_equals(i0_view%get_value_int(), 1)

    call ds2%print()
    
    call ds%delete()
    call ds2%delete()
  end subroutine save_restore_simple

  !------------------------------------------------------------------------------
  subroutine save_restore_complex
    type(datastore) ds, ds2
    type(datagroup) root, flds, root2
    type(datagroup) ga, gb, gc
    type(dataview) i0_view, f0_view, d0_view
    
    ds = datastore_new()
    root = ds%get_root()
    flds = root%create_group("fields")

    ga = flds%create_group("a")
    gb = flds%create_group("b")
    gc = flds%create_group("c")

    i0_view = ga%create_view_and_allocate("i0", ATK_C_INT_T, 1)
    call i0_view%set_value(1)

    f0_view = gb%create_view_and_allocate("f0", ATK_C_FLOAT_T, 1)
    call f0_view%set_value(100.0)

    d0_view = gc%create_view_and_allocate("d0", ATK_C_DOUBLE_T, 1)
    call d0_view%set_value(3000.0d0)

    ! check that all sub groups exist
    call assert_true(flds%has_group("a"))
    call assert_true(flds%has_group("b"))
    call assert_true(flds%has_group("c"))

    call ds%print()

    call root%save("F_out_sidre_group_save_restore_complex","conduit")

    ds2 = datastore_new()
    root2 = ds2%get_root()

    call root2%load("F_out_sidre_group_save_restore_complex","conduit")

    flds = root2%get_group("fields")
    ! check that all sub groups exist
    call assert_true(flds%has_group("a"))
    call assert_true(flds%has_group("b"))
    call assert_true(flds%has_group("c"))
    
    ga = flds%get_group("a");
    gb = flds%get_group("b");
    gc = flds%get_group("c");

    i0_view = ga%get_view("i0");
    f0_view = gb%get_view("f0");
    d0_view = gc%get_view("d0");

    call assert_equals(i0_view%get_value_int(), 1)
    call assert_equals(f0_view%get_value_float(), 100.0)
    call assert_equals(d0_view%get_value_double(), 3000.0d0)

    call ds2%print()

    call ds%delete()
    call ds2%delete()
  end subroutine save_restore_complex

!----------------------------------------------------------------------
end module sidre_group
!----------------------------------------------------------------------

program fortran_test
  use fruit
  use sidre_group
  implicit none
  logical ok

  call init_fruit

  call get_name
  call get_parent
  call get_datastore
  call has_group
  call has_view
  call get_view_name_index
  call get_group_name_index
  call create_destroy_has_view
  call create_destroy_has_group
  call group_name_collisions
  call view_copy_move
  call groups_move_copy
  call create_destroy_view_and_data
  call create_destroy_alloc_view_and_data
  call create_view_of_buffer_with_schema
  call save_restore_simple
  call save_restore_complex

  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif
end program fortran_test

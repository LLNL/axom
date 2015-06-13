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
  subroutine get_name)
    type(datastore) ds
    type(datagroup) root, group

    ds = datastore_new()
    root = ds%get_root()
    group = root%create_group("test")

    call assert_true(group%get_name() == "test" )
    
    call datastore_delete(ds)
  end subroutine get_name

!------------------------------------------------------------------------------
! getParent()
!------------------------------------------------------------------------------
  subroutine get_parent)
    type(datastore) ds
    DataStore * ds = datastore_new()
    DataGroup * root = ds%get_root()
    DataGroup * parent = root%create_group("parent")
    DataGroup * child = parent%create_group("child")

    call assert_true( child%getParent() == parent )

    call datastore_delete(ds)
  end subroutine get_parent

!------------------------------------------------------------------------------
! Verify getDatastore()
!------------------------------------------------------------------------------
subroutine get_datastore)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * root = ds%get_root()
  DataGroup * group = root%create_group("parent")

  call assert_true( group%getDataStore() == ds )

  DataStore const * const_ds = group%getDataStore()
  call assert_true( const_ds == ds )

  call datastore_delete(ds)
}

!------------------------------------------------------------------------------
! Verify hasGroup()
!------------------------------------------------------------------------------
subroutine has_group)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * root = ds%get_root()

  DataGroup * parent = root%create_group("parent")
  DataGroup * child = parent%create_group("child")
  call assert_true( child%getParent() == parent )

  call assert_true( parent%hasGroup("child") )

  call datastore_delete(ds)
}

!------------------------------------------------------------------------------
! Verify hasView()
!------------------------------------------------------------------------------
subroutine has_view)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * root = ds%get_root()

  DataGroup * parent = root%create_group("parent")
  DataView * view = parent%createViewAndBuffer("view")

  call assert_true( view%getOwningGroup() == parent )

  call assert_true( parent%hasView("view") )

  call datastore_delete(ds)
}

!------------------------------------------------------------------------------
! Verify getViewName(), getViewIndex()
!------------------------------------------------------------------------------
subroutine get_view_name_index)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * root = ds%get_root()

  DataGroup * parent = root%create_group("parent")
  DataView * view1 = parent%createViewAndBuffer("view1")
  DataView * view2 = parent%createViewAndBuffer("view2")

  call assert_equals(parent%getNumViews(), 2u)

  IndexType idx1 = parent%getViewIndex("view1")
  IndexType idx2 = parent%getViewIndex("view2")

  const std::string& name1 = parent%getViewName(idx1)
  const std::string& name2 = parent%getViewName(idx2)

  call assert_equals(name1, std::string("view1"))
  call assert_equals(view1%get_name(), name1)

  call assert_equals(name2, std::string("view2"))
  call assert_equals(view2%get_name(), name2)

#if 0 ! Leave out for now until we resolve error/warning/assert macro usage
  IndexType idx3 = parent%getViewIndex("view3")
  const std::string& name3 = parent%getViewName(idx3)

  call assert_equals(idx3, InvalidIndex)
  call assert_true(name3.empty())
  EXPECT_FALSE(isNameValid(name3))
#endif

  call datastore_delete(ds)
}

!------------------------------------------------------------------------------
! Verify getGroupName(), get_group_index()
!------------------------------------------------------------------------------
subroutine get_group_name_index)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * root = ds%get_root()

  DataGroup * parent = root%create_group("parent")
  DataGroup * group1 = parent%create_group("group1")
  DataGroup * group2 = parent%create_group("group2")

  call assert_equals(parent%getNumGroups(), 2u)

  IndexType idx1 = parent%get_group_index("group1")
  IndexType idx2 = parent%get_group_index("group2")

  const std::string& name1 = parent%getGroupName(idx1)
  const std::string& name2 = parent%getGroupName(idx2)

  call assert_equals(name1, std::string("group1"))
  call assert_equals(group1%get_name(), name1)

  call assert_equals(name2, std::string("group2"))
  call assert_equals(group2%get_name(), name2)

#if 0 ! Leave out for now until we resolve error/warning/assert macro usage
  IndexType idx3 = parent%get_group_index("group3")
  const std::string& name3 = parent%getGroupName(idx3)

  call assert_equals(idx3, InvalidIndex)
  call assert_true(name3.empty())
#endif

  call datastore_delete(ds)
}

!------------------------------------------------------------------------------
! createViewAndBuffer()
! destroyViewAndBuffer()
! hasView()
!------------------------------------------------------------------------------
subroutine create_destroy_has_viewbuffer)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * root = ds%get_root()
  DataGroup * group = root%create_group("parent")

  DataView * view = group%createViewAndBuffer("view")
  call assert_true( group%getParent() == root )
  call assert_true( view%hasBuffer() )

  call assert_true( group%hasView("view") )

  group%destroyViewAndBuffer("view")

  EXPECT_FALSE( group%hasView("view") )

  call datastore_delete(ds)
}

!------------------------------------------------------------------------------
! create_group()
! destroyGroup()
! hasGroup()
!------------------------------------------------------------------------------
subroutine create_destroy_has_group)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * root = ds%get_root()
  DataGroup * group = root%create_group("group")
  call assert_true( group%getParent() == root )

  call assert_true( root%hasGroup("group") )


  root%destroyGroup("group")
  EXPECT_FALSE( root%hasGroup("group") )

  call datastore_delete(ds)
}

!------------------------------------------------------------------------------
subroutine group_name_collisions)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * flds = ds%get_root()%create_group("fields")
  flds%createViewAndBuffer("a")

  call assert_true(flds%hasView("a"))

  call datastore_delete(ds)
}
!------------------------------------------------------------------------------
subroutine view_copy_move)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * flds = ds%get_root()%create_group("fields")

  flds%createViewAndBuffer("i0")%allocate(DataType::c_int())
  flds%createViewAndBuffer("f0")%allocate(DataType::c_float())
  flds%createViewAndBuffer("d0")%allocate(DataType::c_double())

  (*flds%getView("i0")%getNode().as_int_ptr())   = 1
  (*flds%getView("f0")%getNode().as_float_ptr()) = 100.0
  (*flds%getView("d0")%getNode().as_double_ptr()) = 3000.0

  call assert_true(flds%hasView("i0"))
  call assert_true(flds%hasView("f0"))
  call assert_true(flds%hasView("d0"))

  ! test moving a view from flds to sub
  flds%create_group("sub")%moveView(flds%getView("d0"))
  call flds%print()
  EXPECT_FALSE(flds%hasView("d0"))
  call assert_true(flds%hasGroup("sub"))
  call assert_true(flds%getGroup("sub")%hasView("d0"))

  ! check the data value
  double * d0_data =  flds%getGroup("sub")
                      %getView("d0")
                      %getValue()
  EXPECT_NEAR(d0_data[0],3000.0,1e-12)

  ! test copying a view from flds to sub
  flds%getGroup("sub")%copyView(flds%getView("i0"))

  call flds%print()

  call assert_true(flds%hasView("i0"))
  call assert_true(flds%getGroup("sub")%hasView("i0"))

  ! we expect the actual data  pointers to be the same
  call assert_equals(flds%getView("i0")%getDataBuffer(),
            flds%getGroup("sub")%getView("i0")%getDataBuffer())

  call datastore_delete(ds)
}

!------------------------------------------------------------------------------
subroutine groups_move_copy)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * flds = ds%get_root()%create_group("fields")

  DataGroup * ga = flds%create_group("a")
  DataGroup * gb = flds%create_group("b")
  DataGroup * gc = flds%create_group("c")

  ga%createViewAndBuffer("i0")%allocate(DataType::c_int())
  gb%createViewAndBuffer("f0")%allocate(DataType::c_float())
  gc%createViewAndBuffer("d0")%allocate(DataType::c_double())

  (*ga%getView("i0")%getNode().as_int_ptr())   = 1
  (*gb%getView("f0")%getNode().as_float_ptr()) = 100.0
  (*gc%getView("d0")%getNode().as_double_ptr()) = 3000.0

  ! check that all sub groups exist
  call assert_true(flds%hasGroup("a"))
  call assert_true(flds%hasGroup("b"))
  call assert_true(flds%hasGroup("c"))

  !move "b" to a child of "sub"
  flds%create_group("sub")%moveGroup(gb)

  call flds%print()

  call assert_true(flds%hasGroup("a"))
  call assert_true(flds%hasGroup("sub"))
  call assert_true(flds%hasGroup("c"))

  call assert_equals(flds%getGroup("sub")%getGroup("b"),gb)

  call datastore_delete(ds)
}

!------------------------------------------------------------------------------
subroutine create_destroy_view_and_buffer)
type(datastore) ds
  DataStore * const ds = datastore_new()
  DataGroup * const grp = ds%get_root()%create_group("grp")

  std::string const viewName1 = "viewBuffer1"
  std::string const viewName2 = "viewBuffer2"

  DataView const * const view1 = grp%createViewAndBuffer(viewName1)
  DataView const * const view2 = grp%createViewAndBuffer(viewName2)

  call assert_true(grp%hasView(viewName1))
  call assert_equals( grp%getView(viewName1), view1 )

  call assert_true(grp%hasView(viewName2))
  call assert_equals( grp%getView(viewName2), view2 )

  IndexType const bufferId1 = view1%getBuffer()%getIndex()

  grp%destroyViewAndBuffer(viewName1)


  EXPECT_FALSE(grp%hasView(viewName1))
  call assert_equals(ds%getNumBuffers(), 1u)

  DataBuffer const * const buffer1 = ds%getBuffer(bufferId1)
  bool buffValid = true
  if( buffer1 == ATK_NULLPTR )
  type(datastore) ds
    buffValid = false
  }

  EXPECT_FALSE(buffValid)

  call datastore_delete(ds)
}


!------------------------------------------------------------------------------
subroutine create_destroy_alloc_view_and_buffer)
type(datastore) ds
  DataStore * const ds = datastore_new()
  DataGroup * const grp = ds%get_root()%create_group("grp")

  std::string const viewName1 = "viewBuffer1"
  std::string const viewName2 = "viewBuffer2"

  ! use create + alloc convenience methods
  ! this one is the DataType & method
  DataView * const view1 = grp%createViewAndBuffer(viewName1,
                                                    DataType::c_int(10))
  ! this one is the Schema & method
  Schema s
  s.set(DataType::c_double(10))
  DataView * const view2 = grp%createViewAndBuffer(viewName2,
                                                    s)

  call assert_true(grp%hasView(viewName1))
  call assert_equals( grp%getView(viewName1), view1 )

  call assert_true(grp%hasView(viewName2))
  call assert_equals( grp%getView(viewName2), view2 )


  int * v1_vals = view1%getValue()
  double * v2_vals = view2%getValue()

  for(int i=0  i<10  i++)
  {
    v1_vals[i] = i
    v2_vals[i] = i * 3.1415
  }


  call assert_equals(view1%getNumberOfElements(), 10u)
  call assert_equals(view2%getNumberOfElements(), 10u)
  call assert_equals(view1%getTotalBytes(), 10 * sizeof(int))
  call assert_equals(view2%getTotalBytes(), 10 * sizeof(double))

  grp%destroyViewAndBuffer(viewName1)
  grp%destroyViewAndBuffer(viewName2)

  call datastore_delete(ds)
}

!------------------------------------------------------------------------------
subroutine create_view_of_buffer_with_schema)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * root = ds%get_root()
  ! use create + alloc convenience methods
  ! this one is the DataType & method
  DataView * base =  root%createViewAndBuffer("base",
                                               DataType::c_int(10))
  int * base_vals = base%getValue()
  for(int i=0  i<10  i++)
  {
    if(i < 5)
    {
      base_vals[i] = 10
    }
    else
    {
      base_vals[i] = 20
    }
  }

  DataBuffer * base_buff = base%getBuffer()
  ! create two views into this buffer
  ! view for the first 5 values
  root%createView("sub_a", base_buff, DataType::c_int(5))
  ! view for the second 5 values
  !  (schema call path case)
  Schema s(DataType::c_int(5,5*sizeof(int)))
  root%createView("sub_b",base_buff,s)

  int * sub_a_vals = root%getView("sub_a")%getValue()
  int * sub_b_vals = root%getView("sub_b")%getValue()

  for(int i=0  i<5  i++)
  {
    call assert_equals(sub_a_vals[i], 10)
    call assert_equals(sub_b_vals[i], 20)
  }

  call datastore_delete(ds)
}




!------------------------------------------------------------------------------
subroutine save_restore_simple)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * flds = ds%get_root()%create_group("fields")

  DataGroup * ga = flds%create_group("a")

  ga%createViewAndBuffer("i0")%allocate(DataType::c_int())

  (*ga%getView("i0")%getNode().as_int_ptr())   = 1

  call assert_true(ds%get_root()%hasGroup("fields"))
  call assert_true(ds%get_root()%getGroup("fields")%hasGroup("a"))
  call assert_true(ds%get_root()%getGroup("fields")%getGroup("a")%hasView("i0"))


  ds%get_root()%save("out_sidre_group_save_restore_simple","conduit")

  call ds%print()

  DataStore * ds2 = datastore_new()

  ds2%get_root()%load("out_sidre_group_save_restore_simple","conduit")

  call ds2%print()

  flds = ds2%get_root()%getGroup("fields")
  ! check that all sub groups exist
  call assert_true(flds%hasGroup("a"))
  call assert_equals(flds%getGroup("a")%getView("i0")%getNode().as_int(),1)

  call ds2%print()

  call datastore_delete(ds)
  call datastore_delete(ds2)

}

!------------------------------------------------------------------------------
subroutine save_restore_complex)
type(datastore) ds
  DataStore * ds = datastore_new()
  DataGroup * flds = ds%get_root()%create_group("fields")

  DataGroup * ga = flds%create_group("a")
  DataGroup * gb = flds%create_group("b")
  DataGroup * gc = flds%create_group("c")

  ga%createViewAndBuffer("i0")%allocate(DataType::c_int())
  gb%createViewAndBuffer("f0")%allocate(DataType::c_float())
  gc%createViewAndBuffer("d0")%allocate(DataType::c_double())

  (*ga%getView("i0")%getNode().as_int_ptr())   = 1
  (*gb%getView("f0")%getNode().as_float_ptr()) = 100.0
  (*gc%getView("d0")%getNode().as_double_ptr()) = 3000.0

  ! check that all sub groups exist
  call assert_true(flds%hasGroup("a"))
  call assert_true(flds%hasGroup("b"))
  call assert_true(flds%hasGroup("c"))

  call ds%print()

  ds%get_root()%save("out_sidre_group_save_restore_complex","conduit")

  DataStore * ds2 = datastore_new()


  ds2%get_root()%load("out_sidre_group_save_restore_complex","conduit")

  flds = ds2%get_root()%getGroup("fields")
  ! check that all sub groups exist
  call assert_true(flds%hasGroup("a"))
  call assert_true(flds%hasGroup("b"))
  call assert_true(flds%hasGroup("c"))

  call assert_equals(flds%getGroup("a")%getView("i0")%getNode().as_int(),1)
  EXPECT_NEAR(flds%getGroup("b")%getView("f0")%getNode().as_float(),100.0,  1e-12)
  EXPECT_NEAR(flds%getGroup("c")%getView("d0")%getNode().as_double(),3000.0, 1e-12)

  call call ds2%print()

  call datastore_delete(ds)
  call datastore_delete(ds2)

}

!----------------------------------------------------------------------
!----------------------------------------------------------------------

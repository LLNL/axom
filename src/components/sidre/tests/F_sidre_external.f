!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!/

module sidre_external
  use iso_c_binding
  use fruit
  use sidre_mod
  implicit none

contains

!------------------------------------------------------------------------------
! Test DataBuffer::declareExternal()
!------------------------------------------------------------------------------
  subroutine declare_external_buffer
    type(datastore) ds
    type(databuffer) dbuff_0, dbuff_1, dbuff_2
    integer(C_INT), allocatable :: idata(:)
    real(C_DOUBLE), allocatable :: ddata(:)
    integer, parameter :: len = 11
    integer ii

    ds = datastore_new()

    allocate(idata(len))
    allocate(ddata(len))

    do ii = 1, len
       idata(ii) = ii
       ddata(ii) = idata(ii) * 2.0
    enddo

    dbuff_0 = ds%createBuffer()
    dbuff_1 = ds%createBuffer()
    dbuff_2 = ds%createBuffer()

    call dbuff_0%allocate(ATK_C_DOUBLE_T, len)
    call dbuff_1%declare(ATK_C_INT, len)
    call dbuff_1%setExternalData(idata)
    call dbuff_2%declare(ATK_C_DOUBLE_T, len)
    call dbuff_2%setExternalData(ddata)

    call assert_equals(dbuff_0%isExternal(), .false.)
    call assert_equals(dbuff_1%isExternal(), .true.)
    call assert_equals(dbuff_2%isExternal(), .true.)

    call assert_equals(dbuff_0%getTotalBytes(), sizeof(CONDUIT_NATIVE_DOUBLE)*len)
    call assert_equals(dbuff_1%getTotalBytes(), sizeof(CONDUIT_NATIVE_INT)*len)
    call assert_equals(dbuff_2%getTotalBytes(), sizeof(CONDUIT_NATIVE_DOUBLE)*len)

    call ds%print()

    call ds%delete()
    deallocate(idata)
    deallocate(ddata)
  end subroutine declare_external_buffer

!------------------------------------------------------------------------------
! Test DataGroup::createExternalView()
!------------------------------------------------------------------------------
  subroutine create_external_view
    type(datastore) ds
    type(dataroot) root
    type(dataview) iview, dview
    integer(C_INT), allocatable :: idata(:)
    real(C_DOUBLE), allocatable :: ddata(:)
    integer(C_INT), allocatable :: idata_chk(:)
    real(C_DOUBLE), allocatable :: ddata_chk(:)
    integer, parameter :: len = 11
    integer ii

    ds = datastore_new()
    root = ds%getRoot()

    allocate(idata(len))
    allocate(ddata(len))

    do ii = 1, len
       idata(ii) = ii
       ddata(ii) = idata(ii) * 2.0
    enddo

    iview = root%createExternalView("idata", idata, ATK_C_INT_T, len)
    dview = root%createExternalView("ddata", ddata, ATK_C_DOUBLE_T, len)
    call assert_true(root%getNumViews() .eq. 2)

    call iview%print()
    call dview%print()

    call iview%getValue(idata_chk)
    do ii = 1, len
       call assert_equals(idata_chk(ii), idata(ii))
    enddo

    call dview%getValue(ddata_chk)
    do ii = 1, len
       call assert_equals(ddata_chk(ii), ddata(ii))
    enddo

    delete ds
    deallocate(idata)
    deallocate(ddata)
  end subroutine create_external_view

!------------------------------------------------------------------------------
! Test DataGroup::save(), DataGroup::load() with external buffers
!------------------------------------------------------------------------------
  subroutine save_load_external_view
    type(datastore) ds, ds2
    type(datagroup) root, root2
    type(dataview)  iview, dview
    integer(C_INT), allocatable :: idata(:)
    real(C_DOUBLE), allocatable :: ddata(:)
    integer(C_INT), pointer :: idata_chk(:)
    real(C_DOUBLE), pointer :: ddata_chk(:)
    integer, parameter :: len = 11
    integer ii

    ds = datastore_new()
    root = ds%getRoot()

    allocate(idata(len))
    allocate(ddata(len))

    do ii = 1, len
       idata(ii) = ii
       ddata(ii) = idata(ii) * 2.0
    enddo

    iview = root%createExternalView("idata", idata, ATK_INT_C_T, len)
    dview = root%createExternalView("ddata", ddata, ATK_DOUBLE_C_T, len)
    call assert_equals(root%getNumViews(), 2)
    call assert_equals(root%getView("idata")%getBuffer()%isExternal(), .true.)
    call assert_equals(root%getView("ddata")%getBuffer()%isExternal(), .true.)

    call iview%print()
    call dview%print()

    call root%save("out_sidre_external_save_restore_external_view", "conduit")

    ds%print()


    ds2 = datastore_new()
    root2 = ds2%getRoot()

    call ds2%load("out_sidre_external_save_restore_external_view","conduit")

    call ds2%print()

    call assert_equals(root2%getNumViews(), 2)
    call assert_equals(root2%getView("idata")%getBuffer()%isExternal(), .false.)
    call assert_equals(root2%getView("ddata")%getBuffer()%isExternal(), .false.)

    call root2%getView("idata")%getValue(idata_chk)
    do ii = 1, len
       call assert_equals(idata_chk(ii), idata(ii))
    enddo

    call root2%getView("ddata")%getValue(ddata_chk)
    do ii = 1, len
       call assert_equals(ddata_chk(ii), ddata(ii))
    enddo

    call ds%delete()
    call ds2%delete()
    deallocate(idata)
    deallocate(ddata)
  end subroutine save_load_external_view

!----------------------------------------------------------------------
end module sidre_external
!----------------------------------------------------------------------

function fortran_test() bind(C,name="fortran_test")
  use fruit
  use sidre_group
  implicit none
  integer(C_INT) fortran_test
  logical ok

  call init_fruit

  call declare_external_buffer
  call create_external_view
  call save_load_external_view


  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (ok) then
     fortran_test = 0
  else
     fortran_test = 1
  endif
end function fortran_test

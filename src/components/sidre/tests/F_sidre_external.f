!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!

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
    integer(C_INT), allocatable, target :: idata(:)
    real(C_DOUBLE), allocatable, target :: ddata(:)
    integer, parameter :: len = 11
    integer ii

    call set_case_name("declare_external_buffer")

    ds = datastore_new()

    allocate(idata(len))
    allocate(ddata(len))

    do ii = 1, len
       idata(ii) = ii
       ddata(ii) = idata(ii) * 2.0
    enddo

    dbuff_0 = ds%create_buffer()
    dbuff_1 = ds%create_buffer()
    dbuff_2 = ds%create_buffer()

    call dbuff_0%allocate(SIDRE_DOUBLE_ID, len)
    call dbuff_1%declare(SIDRE_INT_ID, len)
    call dbuff_1%set_external_data(c_loc(idata))
    call dbuff_2%declare(SIDRE_DOUBLE_ID, len)
    call dbuff_2%set_external_data(c_loc(ddata))

    call assert_equals(dbuff_0%is_external(), .false.)
    call assert_equals(dbuff_1%is_external(), .true.)
    call assert_equals(dbuff_2%is_external(), .true.)

!--    call assert_equals(dbuff_0%get_total_bytes(), sizeof(CONDUIT_NATIVE_DOUBLE)*len)
!--    call assert_equals(dbuff_1%get_total_bytes(), sizeof(CONDUIT_NATIVE_INT)*len)
!--    call assert_equals(dbuff_2%get_total_bytes(), sizeof(CONDUIT_NATIVE_DOUBLE)*len)

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
    type(datagroup) root
    type(dataview) iview, dview
    integer(C_INT), allocatable, target :: idata(:)
    real(C_DOUBLE), allocatable, target :: ddata(:)
    integer(C_INT), pointer :: idata_chk(:)
    real(C_DOUBLE), pointer :: ddata_chk(:)
    integer(C_LONG), parameter :: len = 11
    integer ii

    call set_case_name("create_external_view")

    ds = datastore_new()
    root = ds%get_root()

    allocate(idata(len))
    allocate(ddata(len))

    do ii = 1, len
       idata(ii) = ii
       ddata(ii) = idata(ii) * 2.0
    enddo

    iview = root%create_view_external("idata", c_loc(idata))
    call iview%apply(SIDRE_INT_ID, len)
    dview = root%create_view_external("ddata", c_loc(ddata))
    call dview%apply(SIDRE_DOUBLE_ID, len)
    call assert_true(root%get_num_views() .eq. 2)

    call iview%print()
    call dview%print()

    call iview%get_data(idata_chk)
    do ii = 1, len
       call assert_equals(idata_chk(ii), idata(ii))
    enddo

    call dview%get_data(ddata_chk)
    do ii = 1, len
       call assert_equals(ddata_chk(ii), ddata(ii))
    enddo

    call ds%delete()
    deallocate(idata)
    deallocate(ddata)
  end subroutine create_external_view

!------------------------------------------------------------------------------
! Test DataGroup::save(), DataGroup::load() with external buffers
!------------------------------------------------------------------------------
  subroutine save_load_external_view
    type(datastore) ds, ds2
    type(datagroup) root, root2
    type(dataview)  iview, dview, iview2, dview2
    type(databuffer) tmpbuff
    integer(C_INT), allocatable, target :: idata(:)
    real(C_DOUBLE), allocatable, target :: ddata(:)
    integer(C_INT), pointer :: idata_chk(:)
    real(C_DOUBLE), pointer :: ddata_chk(:)
    integer(C_LONG), parameter :: len = 11
    integer ii

    call set_case_name("save_load_external_view")

    ds = datastore_new()
    root = ds%get_root()

    allocate(idata(len))
    allocate(ddata(len))

    do ii = 1, len
       idata(ii) = ii
       ddata(ii) = idata(ii) * 2.0
    enddo

    iview = root%create_view_external("idata", c_loc(idata))
    call iview%apply_type_nelems(SIDRE_INT_ID, len)
    dview = root%create_view_external("ddata", c_loc(ddata))
    call dview%apply_type_nelems(SIDRE_DOUBLE_ID, len)

    call assert_true(root%get_num_views() .eq. 2)
    tmpbuff = iview%get_buffer()
    call assert_equals(tmpbuff%is_external(), .true.)
    tmpbuff = dview%get_buffer()
    call assert_equals(tmpbuff%is_external(), .true.)

    call iview%print()
    call dview%print()

    call root%save("out_sidre_external_save_restore_external_view", "conduit")

    call ds%print()


    ds2 = datastore_new()
    root2 = ds2%get_root()

    call root2%load("out_sidre_external_save_restore_external_view","conduit")

    call ds2%print()

    call assert_true(root2%get_num_views() .eq. 2)

    iview = root2%get_view("idata")
    dview = root2%get_view("ddata")

    tmpbuff = iview%get_buffer()
    call assert_equals(tmpbuff%is_external(), .false.)

    tmpbuff = dview%get_buffer()
    call assert_equals(tmpbuff%is_external(), .false.)

    call iview%get_data(idata_chk)
    do ii = 1, len
       call assert_equals(idata_chk(ii), idata(ii))
    enddo

    call dview%get_data(ddata_chk)
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

program fortran_test
  use fruit
  use sidre_external
  implicit none
  logical ok

  call init_fruit

  call declare_external_buffer
  call create_external_view
  call save_load_external_view


  call fruit_summary
  call fruit_finalize

  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif
end program fortran_test

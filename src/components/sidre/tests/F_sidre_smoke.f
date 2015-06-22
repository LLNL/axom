!
!  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
!  Produced at the Lawrence Livermore National Laboratory.
! 
!  All rights reserved.
! 
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!

module sidre_smoke
  use fruit
  use sidre_mod
  implicit none

contains
!------------------------------------------------------------------------------

  subroutine create_datastore
    type(datastore) ds

    ds = datastore_new()
    call datastore_delete(ds)

    call assert_true(.true.)
  end subroutine create_datastore

!------------------------------------------------------------------------------

  subroutine valid_invalid
    type(datastore) ds
    type(datagroup) root
    integer idx
    character(10) name

    ds = datastore_new()

    idx = 3;
    call assert_true(idx /= invalid_index)

    name = "foo"
    call assert_true(is_name_valid(name))

    root = ds%get_root()

    call assert_true(root%get_group_name(idx) == " ")
    call assert_true(root%get_group_index(name) == invalid_index)

    call datastore_delete(ds)
  end subroutine valid_invalid

!------------------------------------------------------------------------------
end module sidre_smoke
!------------------------------------------------------------------------------

program tester
  use fruit
  use sidre_smoke
  call init_fruit

  call create_datastore

  call fruit_summary
  call fruit_finalize
end program tester

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
  implicit none

contains
!------------------------------------------------------------------------------

  subroutine create_datastore
    use sidre_mod
    type(datastore) ds

    ds = datastore_new()
    call datastore_delete(ds)

    call assert_true(.true.)
  end subroutine create_datastore

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

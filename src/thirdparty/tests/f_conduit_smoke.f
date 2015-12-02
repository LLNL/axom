!------------------------------------------------------------------------------
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
! f_conduit_smoke.f
!
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
module f_conduit_smoke
!------------------------------------------------------------------------------

  use iso_c_binding
  use fruit
  use conduit
  implicit none

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine t_conduit_smoke
        type(C_PTR) cnode
        integer res
        !--------------
        ! c++ ~equiv:
        !--------------
        ! Node n;
        ! n.print_detailed();
        cnode = conduit_node_create()
        call conduit_node_print_detailed(cnode)
        call conduit_node_destroy(cnode)
    
    end subroutine t_conduit_smoke

!------------------------------------------------------------------------------
end module f_conduit_smoke
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
program fortran_test
!------------------------------------------------------------------------------
  use fruit
  use f_conduit_smoke
  implicit none
  logical ok
  
  call init_fruit

  !----------------------------------------------------------------------------
  ! call our test
  !----------------------------------------------------------------------------
  call t_conduit_smoke
  
  call fruit_summary
  call fruit_finalize
  call is_all_successful(ok)
  if (.not. ok) then
     call exit(1)
  endif

!------------------------------------------------------------------------------
end program fortran_test
!------------------------------------------------------------------------------



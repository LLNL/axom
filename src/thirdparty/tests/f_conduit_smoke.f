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
function fortran_test() bind(C,name="fortran_test")
!------------------------------------------------------------------------------
  use fruit
  use f_conduit_smoke
  implicit none
  integer(C_INT) fortran_test
  logical res
  
  call init_fruit

  !----------------------------------------------------------------------------
  ! call our test
  !----------------------------------------------------------------------------
  call t_conduit_smoke
  
  call fruit_summary
  call fruit_finalize
  call is_all_successful(res)
  if (res) then
     fortran_test = 0
  else
     fortran_test = 1
  endif

!------------------------------------------------------------------------------
end function fortran_test
!------------------------------------------------------------------------------



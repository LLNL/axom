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
! fruit_smoke.f
!
!------------------------------------------------------------------------------
module fruit_smoke
  use iso_c_binding
  use fruit
  implicit none

contains
!------------------------------------------------------------------------------

  subroutine simple_test
        call assert_equals (42, 42)
  end subroutine simple_test


!----------------------------------------------------------------------
end module fruit_smoke
!----------------------------------------------------------------------

function fortran_test() bind(C,name="fortran_test")
  use fruit
  use fruit_smoke
  implicit none
  integer(C_INT) fortran_test
  logical res

  call init_fruit
!----------
! Our tests
  call simple_test
!----------

  call fruit_summary
  call fruit_finalize
  call is_all_successful(res)
  if (res) then
     fortran_test = 0
  else
     fortran_test = 1
  endif
  
end function fortran_test


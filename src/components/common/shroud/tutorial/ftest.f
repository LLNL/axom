!
! Test tutorial wrapper
!

program tester
  use iso_c_binding
  use tutorial_mod

  type(class1) cptr
  real(C_DOUBLE) rv2
  logical rv3
  character(30) rv4

  call function1

  cptr = class1_new()
  call cptr%method1

  ! Arguments
  ! Integer and Real
  rv2 = function2(1.5d0, 2)

  ! logical
  rv3 = function3(.false.)

  ! character
  rv4 = function4("bird", "dog")
  print *, rv4



end program tester

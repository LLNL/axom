!
! Test tutorial wrapper
!

program tester
  use iso_c_binding
  use tutorial_mod

  type(class1) cptr
  real(C_DOUBLE) rv2

  call function1

  cptr = class1_new()
  call cptr%method1

  ! Arguments
  ! Integer and Real
  rv2 = function2(1.5d0, 2)


end program tester

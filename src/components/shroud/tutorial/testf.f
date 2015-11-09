!
! Test tutorial wrapper
!

program tester
  use iso_c_binding
  use tutorial_mod
  implicit none

  type(class1) cptr
  real(C_DOUBLE) rv2
  logical rv3
  character(30) rv4, rv4b

  call function1

  cptr = class1_new()
  call cptr%method1

  ! Arguments
  ! Integer and Real
  rv2 = function2(1.5d0, 2)

  ! logical
  rv3 = function3(.false.)

  ! character
  rv4 = function4a("bird", "dog")
  print *, rv4
  call function4b("bird", "dog", rv4b)
  print *, rv4b

  print *, "function5", function5()
  print *, "function5", function5(0.0d0)
  print *, "function5", function5(arg2=0)
  print *, "functino5", function5(2.0d0, 2)

  call function6_from_name("name")
  call function6_from_index(1)
  call function6("name")
  call function6(1)

  call function7(1)
  call function7(1.0d0)

  call function9(1.0)
  call function9(1.0d0)

end program tester

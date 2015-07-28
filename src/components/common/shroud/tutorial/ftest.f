!
! Test tutorial wrapper
!

program tester
  use tutorial_mod

  type(class1) cptr

  call function1

  cptr = class1_new()
  call cptr%method1



end program tester

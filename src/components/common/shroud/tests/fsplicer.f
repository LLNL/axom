

! splicer push class.exclass1

! splicer begin module_top
top of module splicer  1
! splicer end module_top

! splicer begin component_part
  component part 1a
  component part 1b
! splicer end   component_part

! splicer begin type_bound_procedure_part
  type bound procedure part 1
! splicer end   type_bound_procedure_part

! splicer push method

! splicer begin splicer_special
blah blah blah
! splicer end splicer_special

! splicer pop method


! splicer begin extra_methods
  insert extra methods here
! splicer end   extra_methods

! splicer pop class.exclass1


! splicer push class.exclass2
! splicer begin module_top
top of module splicer  2
! splicer end module_top
! splicer pop class.exclass2




# test a full path
! splicer begin  class.exclass1.method.extra_method2
  ! extra method 2
! splicer end    class.exclass1.method.extra_method2



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

! splicer pop class.exclass1.method



! splicer push class.exclass2
! splicer begin module_top
top of module splicer  2
! splicer end module_top
! splicer pop class.exclass2

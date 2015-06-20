

! splicer begin class.ExClass1.module_top
top of module splicer  1
! splicer end class.ExClass1.module_top

! splicer begin class.ExClass1.component_part
  component part 1a
  component part 1b
! splicer end   class.ExClass1.component_part

! splicer begin class.ExClass1.type_bound_procedure_part
  type bound procedure part 1
! splicer end   class.ExClass1.type_bound_procedure_part

! splicer begin class.ExClass1.method.splicer_special
blah blah blah
! splicer end class.ExClass1.method.splicer_special

! splicer begin class.ExClass1.extra_methods
  insert extra methods here
! splicer end   class.ExClass1.extra_methods


! splicer begin class.ExClass2.module_top
top of module splicer  2
! splicer end class.ExClass2.module_top




# test a full path
! splicer begin  class.ExClass1.method.extra_method2
  ! extra method 2
! splicer end    class.ExClass1.method.extra_method2

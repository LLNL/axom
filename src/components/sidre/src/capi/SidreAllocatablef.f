!
! SidreAllocatablef.f - Routines used by Fortran interface
! Uses cog to insert some generated code into this file.
!

![[[cog
!import cog
!import genfsidresplicer as gen
!gen.print_lines(cog.outl, gen.print_atk_size_allocatable)
!]]]

function atk_size_allocatable_int_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(IN) :: array
    integer(C_SIZE_T) :: rv
    rv = 1
end function atk_size_allocatable_int_scalar_ptr

function atk_size_allocatable_int_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(IN) :: array(:)
    integer(C_SIZE_T) :: rv
    rv = size(array)
end function atk_size_allocatable_int_1d_ptr

function atk_size_allocatable_long_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(IN) :: array
    integer(C_SIZE_T) :: rv
    rv = 1
end function atk_size_allocatable_long_scalar_ptr

function atk_size_allocatable_long_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(IN) :: array(:)
    integer(C_SIZE_T) :: rv
    rv = size(array)
end function atk_size_allocatable_long_1d_ptr

function atk_size_allocatable_float_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(IN) :: array
    integer(C_SIZE_T) :: rv
    rv = 1
end function atk_size_allocatable_float_scalar_ptr

function atk_size_allocatable_float_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(IN) :: array(:)
    integer(C_SIZE_T) :: rv
    rv = size(array)
end function atk_size_allocatable_float_1d_ptr

function atk_size_allocatable_double_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(IN) :: array
    integer(C_SIZE_T) :: rv
    rv = 1
end function atk_size_allocatable_double_scalar_ptr

function atk_size_allocatable_double_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(IN) :: array(:)
    integer(C_SIZE_T) :: rv
    rv = size(array)
end function atk_size_allocatable_double_1d_ptr
![[[end]]]


![[[cog
!gen.print_lines(cog.outl, gen.print_atk_address_allocatable)
!]]]

function atk_address_allocatable_int_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(IN), target :: array
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_int_scalar_ptr

function atk_address_allocatable_int_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_INT), allocatable, intent(IN), target :: array(:)
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_int_1d_ptr

function atk_address_allocatable_long_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(IN), target :: array
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_long_scalar_ptr

function atk_address_allocatable_long_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    integer(C_LONG), allocatable, intent(IN), target :: array(:)
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_long_1d_ptr

function atk_address_allocatable_float_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(IN), target :: array
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_float_scalar_ptr

function atk_address_allocatable_float_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_FLOAT), allocatable, intent(IN), target :: array(:)
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_float_1d_ptr

function atk_address_allocatable_double_scalar_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(IN), target :: array
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_double_scalar_ptr

function atk_address_allocatable_double_1d_ptr(array) result(rv)
    use iso_c_binding
    implicit none
    real(C_DOUBLE), allocatable, intent(IN), target :: array(:)
    type(C_PTR) :: rv
    rv = c_loc(array)
end function atk_address_allocatable_double_1d_ptr
![[[end]]]


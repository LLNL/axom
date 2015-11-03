! Fortran code that will be inserted into code via shroud splicer blocks

! splicer begin module_top
!
! Type parameters
! Must be kept in sync with SidreTypes.h
!
integer, parameter :: ATK_INT8_T = 3
integer, parameter :: ATK_INT16_T = 4
integer, parameter :: ATK_INT32_T = 5
integer, parameter :: ATK_INT64_T = 6
integer, parameter :: ATK_UINT8_T = 7
integer, parameter :: ATK_UINT16_T = 8
integer, parameter :: ATK_UINT32_T = 9
integer, parameter :: ATK_UINT64_T = 10
integer, parameter :: ATK_FLOAT32_T = 11
integer, parameter :: ATK_FLOAT64_T  = 12
integer, parameter :: ATK_CHAR8_STR_T = 13

integer, parameter :: ATK_C_INT_T = 14
integer, parameter :: ATK_C_LONG_T = 15
integer, parameter :: ATK_C_FLOAT_T = 16
integer, parameter :: ATK_C_DOUBLE_T = 17

integer, parameter :: invalid_index = -1

! splicer end module_top


# ATK_create_fortran_allocatable_view is not in api.yaml since it is not in src/core and
# only required for the fortran interface.

! splicer begin additional_interfaces
function ATK_create_fortran_allocatable_view(group, name, lname, addr, itype, rank) &
   bind(C,name="ATK_create_fortran_allocatable_view") &
   result(rv)
      use iso_c_binding
      type(C_PTR), value, intent(IN)    :: group
      character(kind=C_CHAR), intent(IN) :: name(*)
      integer(C_INT), value, intent(IN) :: lname
      type(C_PTR), value                :: addr
      integer(C_INT), value, intent(IN) :: itype
      integer(C_INT), value, intent(IN) :: rank
      type(C_PTR) rv
end function ATK_create_fortran_allocatable_view
! splicer end additional_interfaces

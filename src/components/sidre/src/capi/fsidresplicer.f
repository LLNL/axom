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


------------------ DataView

! splicer begin class.DataView.type_bound_procedure_part
procedure :: dataview_get_value_int_1d_ptr
procedure :: dataview_get_value_float_1d_ptr
generic :: get_value => &
    dataview_get_value_int_1d_ptr, &
    dataview_get_value_float_1d_ptr
! splicer end class.DataView.type_bound_procedure_part


! splicer begin class.DataView.additional_functions

subroutine dataview_get_value_int_1d_ptr(view, value)
    use iso_c_binding
    implicit none
    class(dataview), intent(IN) :: view
    integer(C_INT), pointer, intent(OUT) :: value(:)

    type(C_PTR) cptr
    integer(C_SIZE_T) nelems

    cptr = view%get_data_pointer()
    nelems = view%get_number_of_elements()
    call c_f_pointer(cptr, value, [ nelems ])

end subroutine dataview_get_value_int_1d_ptr

subroutine dataview_get_value_float_1d_ptr(view, value)
    use iso_c_binding
    implicit none
    class(dataview), intent(IN) :: view
    real(C_FLOAT), pointer, intent(OUT) :: value(:)

    type(C_PTR) cptr
    integer(C_SIZE_T) nelems

    cptr = view%get_data_pointer()
    nelems = view%get_number_of_elements()
    call c_f_pointer(cptr, value, [ nelems ])

end subroutine dataview_get_value_float_1d_ptr

! splicer end class.DataView.additional_functions

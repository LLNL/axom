! wrapfstrings.f
! This is generated code, do not edit
!>
!! \file wrapfstrings.f
!! \brief Shroud generated wrapper for strings library
!<
module strings_mod
    use fstr_mod
    ! splicer begin module_use
    ! splicer end module_use
    implicit none
    
    
    interface
        
        pure function str_get_name() &
                result(rv) &
                bind(C, name="STR_get_name")
            use iso_c_binding
            implicit none
            type(C_PTR) rv
        end function str_get_name
        
        subroutine str_get_name_bufferify(output, Loutput) &
                bind(C, name="STR_get_name_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(OUT) :: output(*)
            integer(C_INT), value, intent(IN) :: Loutput
        end subroutine str_get_name_bufferify
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains
    
    ! void getName(string_result_as_arg & output+intent(out)+len(Loutput)) const
    ! string_to_buffer_and_len - string_to_buffer_and_len
    ! function_index=2
    !>
    !! \brief return a string as argument
    !!
    !<
    subroutine get_name(output)
        use iso_c_binding
        implicit none
        character(*), intent(OUT) :: output
        ! splicer begin get_name
        call str_get_name_bufferify(  &
            output,  &
            len(output))
        ! splicer end get_name
    end subroutine get_name
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module strings_mod

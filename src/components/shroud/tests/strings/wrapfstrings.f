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
        
        subroutine str_accept_string_const_reference(arg1) &
                bind(C, name="STR_accept_string_const_reference")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: arg1(*)
        end subroutine str_accept_string_const_reference
        
        subroutine str_accept_string_const_reference_bufferify(arg1, Larg1) &
                bind(C, name="STR_accept_string_const_reference_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: arg1(*)
            integer(C_INT), value, intent(IN) :: Larg1
        end subroutine str_accept_string_const_reference_bufferify
        
        subroutine str_accept_string_reference(arg1) &
                bind(C, name="STR_accept_string_reference")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: arg1(*)
        end subroutine str_accept_string_reference
        
        subroutine str_accept_string_reference_bufferify(arg1, Larg1) &
                bind(C, name="STR_accept_string_reference_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: arg1(*)
            integer(C_INT), value, intent(IN) :: Larg1
        end subroutine str_accept_string_reference_bufferify
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains
    
    ! void getName(string_result_as_arg & output+intent(out)+len(Loutput))+pure
    ! string_to_buffer_and_len - string_to_buffer_and_len
    ! function_index=4
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
    
    ! void acceptStringConstReference(const std::string & arg1+intent(in)+len_trim(Larg1))
    ! string_to_buffer_and_len - string_to_buffer_and_len
    ! function_index=6
    !>
    !! \brief Accept a const string reference
    !!
    !! Save contents of arg1.
    !! arg1 is assumed to be intent(IN) since it is const
    !! Will copy in.
    !<
    subroutine accept_string_const_reference(arg1)
        use iso_c_binding
        implicit none
        character(*), intent(IN) :: arg1
        ! splicer begin accept_string_const_reference
        call str_accept_string_const_reference_bufferify(  &
            arg1,  &
            len_trim(arg1))
        ! splicer end accept_string_const_reference
    end subroutine accept_string_const_reference
    
    ! void acceptStringReference(std::string & arg1+intent(in)+len_trim(Larg1))
    ! string_to_buffer_and_len - string_to_buffer_and_len
    ! function_index=8
    !>
    !! \brief Accept a string reference
    !!
    !! Append "dog" to the end of arg1.
    !! arg1 is assumed to be intent(INOUT)
    !! Must copy in and copy out.
    !<
    subroutine accept_string_reference(arg1)
        use iso_c_binding
        implicit none
        character(*), intent(IN) :: arg1
        ! splicer begin accept_string_reference
        call str_accept_string_reference_bufferify(  &
            arg1,  &
            len_trim(arg1))
        ! splicer end accept_string_reference
    end subroutine accept_string_reference
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module strings_mod

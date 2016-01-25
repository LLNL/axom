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
        
        subroutine str_pass_char_ptr(dest, Ndest, src) &
                bind(C, name="STR_pass_char_ptr")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(OUT) :: dest(*)
            integer(C_INT), value, intent(IN) :: Ndest
            character(kind=C_CHAR), intent(IN) :: src(*)
        end subroutine str_pass_char_ptr
        
        subroutine str_pass_char_ptr_bufferify(dest, Ndest, src, Lsrc) &
                bind(C, name="STR_pass_char_ptr_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(OUT) :: dest(*)
            integer(C_INT), value, intent(IN) :: Ndest
            character(kind=C_CHAR), intent(IN) :: src(*)
            integer(C_INT), value, intent(IN) :: Lsrc
        end subroutine str_pass_char_ptr_bufferify
        
        pure function str_get_char1() &
                result(rv) &
                bind(C, name="STR_get_char1")
            use iso_c_binding
            implicit none
            type(C_PTR) rv
        end function str_get_char1
        
        subroutine str_get_char1_bufferify(SH_F_rv, LSH_F_rv) &
                bind(C, name="STR_get_char1_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv(*)
            integer(C_INT), value, intent(IN) :: LSH_F_rv
        end subroutine str_get_char1_bufferify
        
        function str_get_char2() &
                result(rv) &
                bind(C, name="STR_get_char2")
            use iso_c_binding
            implicit none
            type(C_PTR) rv
        end function str_get_char2
        
        subroutine str_get_char2_bufferify(SH_F_rv, LSH_F_rv) &
                bind(C, name="STR_get_char2_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv(*)
            integer(C_INT), value, intent(IN) :: LSH_F_rv
        end subroutine str_get_char2_bufferify
        
        function str_get_char3() &
                result(rv) &
                bind(C, name="STR_get_char3")
            use iso_c_binding
            implicit none
            type(C_PTR) rv
        end function str_get_char3
        
        subroutine str_get_char3_bufferify(output, Loutput) &
                bind(C, name="STR_get_char3_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(OUT) :: output(*)
            integer(C_INT), value, intent(IN) :: Loutput
        end subroutine str_get_char3_bufferify
        
        pure function str_get_string1() &
                result(rv) &
                bind(C, name="STR_get_string1")
            use iso_c_binding
            implicit none
            type(C_PTR) rv
        end function str_get_string1
        
        subroutine str_get_string1_bufferify(SH_F_rv, LSH_F_rv) &
                bind(C, name="STR_get_string1_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv(*)
            integer(C_INT), value, intent(IN) :: LSH_F_rv
        end subroutine str_get_string1_bufferify
        
        function str_get_string2() &
                result(rv) &
                bind(C, name="STR_get_string2")
            use iso_c_binding
            implicit none
            type(C_PTR) rv
        end function str_get_string2
        
        subroutine str_get_string2_bufferify(SH_F_rv, LSH_F_rv) &
                bind(C, name="STR_get_string2_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv(*)
            integer(C_INT), value, intent(IN) :: LSH_F_rv
        end subroutine str_get_string2_bufferify
        
        function str_get_string3() &
                result(rv) &
                bind(C, name="STR_get_string3")
            use iso_c_binding
            implicit none
            type(C_PTR) rv
        end function str_get_string3
        
        subroutine str_get_string3_bufferify(output, Loutput) &
                bind(C, name="STR_get_string3_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(OUT) :: output(*)
            integer(C_INT), value, intent(IN) :: Loutput
        end subroutine str_get_string3_bufferify
        
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
        
        subroutine str_accept_string_reference(arg1, Narg1) &
                bind(C, name="STR_accept_string_reference")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(INOUT) :: arg1(*)
            integer(C_INT), value, intent(IN) :: Narg1
        end subroutine str_accept_string_reference
        
        subroutine str_accept_string_reference_bufferify(arg1, Larg1, Narg1) &
                bind(C, name="STR_accept_string_reference_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(INOUT) :: arg1(*)
            integer(C_INT), value, intent(IN) :: Larg1
            integer(C_INT), value, intent(IN) :: Narg1
        end subroutine str_accept_string_reference_bufferify
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains
    
    ! void passCharPtr(char * dest+intent(out)+len(Ndest), const char * src+intent(in))
    ! string_to_buffer_and_len
    ! function_index=0
    !>
    !! \brief strcpy like behavior
    !!
    !! dest is marked intent(OUT) to override the intent(INOUT) default
    !! This avoid a copy-in on dest.
    !<
    subroutine pass_char_ptr(dest, src)
        use iso_c_binding
        implicit none
        character(*), intent(OUT) :: dest
        character(*), intent(IN) :: src
        ! splicer begin pass_char_ptr
        call str_pass_char_ptr_bufferify(  &
            dest,  &
            len(dest, kind=C_INT),  &
            src,  &
            len_trim(src, kind=C_INT))
        ! splicer end pass_char_ptr
    end subroutine pass_char_ptr
    
    ! const string_result_fstr * getChar1()+pure
    ! function_index=11
    !>
    !! \brief return a 'const char *' as character(*)
    !!
    !<
    function get_char1() result(rv)
        use iso_c_binding
        implicit none
        character(kind=C_CHAR, len=strlen_ptr(str_get_char1())) :: rv
        ! splicer begin get_char1
        rv = fstr(str_get_char1())
        ! splicer end get_char1
    end function get_char1
    
    ! const char * getChar2()
    ! string_to_buffer_and_len
    ! function_index=2
    !>
    !! \brief return 'const char *' with fixed size (len=30)
    !!
    !<
    function get_char2() result(rv)
        use iso_c_binding
        implicit none
        character(kind=C_CHAR, len=30) :: rv
        ! splicer begin get_char2
        call str_get_char2_bufferify(  &
            rv,  &
            len(rv, kind=C_INT))
        ! splicer end get_char2
    end function get_char2
    
    ! void getChar3(char_result_as_arg * output+intent(out)+len(Loutput))
    ! string_to_buffer_and_len - string_to_buffer_and_len
    ! function_index=14
    !>
    !! \brief return a 'const char *' as argument
    !!
    !<
    subroutine get_char3(output)
        use iso_c_binding
        implicit none
        character(*), intent(OUT) :: output
        ! splicer begin get_char3
        call str_get_char3_bufferify(  &
            output,  &
            len(output, kind=C_INT))
        ! splicer end get_char3
    end subroutine get_char3
    
    ! const string_result_fstr & getString1()+pure
    ! function_index=16
    !>
    !! \brief return a 'const string&' as character(*)
    !!
    !<
    function get_string1() result(rv)
        use iso_c_binding
        implicit none
        character(kind=C_CHAR, len=strlen_ptr(str_get_string1())) :: rv
        ! splicer begin get_string1
        rv = fstr(str_get_string1())
        ! splicer end get_string1
    end function get_string1
    
    ! const string & getString2()
    ! string_to_buffer_and_len
    ! function_index=5
    !>
    !! \brief return 'const string&' with fixed size (len=30)
    !!
    !<
    function get_string2() result(rv)
        use iso_c_binding
        implicit none
        character(kind=C_CHAR, len=30) :: rv
        ! splicer begin get_string2
        call str_get_string2_bufferify(  &
            rv,  &
            len(rv, kind=C_INT))
        ! splicer end get_string2
    end function get_string2
    
    ! void getString3(string_result_as_arg & output+intent(out)+len(Loutput))
    ! string_to_buffer_and_len - string_to_buffer_and_len
    ! function_index=19
    !>
    !! \brief return a 'const string&' as argument
    !!
    !<
    subroutine get_string3(output)
        use iso_c_binding
        implicit none
        character(*), intent(OUT) :: output
        ! splicer begin get_string3
        call str_get_string3_bufferify(  &
            output,  &
            len(output, kind=C_INT))
        ! splicer end get_string3
    end subroutine get_string3
    
    ! void acceptStringConstReference(const std::string & arg1+intent(in))
    ! string_to_buffer_and_len
    ! function_index=7
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
            len_trim(arg1, kind=C_INT))
        ! splicer end accept_string_const_reference
    end subroutine accept_string_const_reference
    
    ! void acceptStringReference(std::string & arg1+intent(inout)+len(Narg1))
    ! string_to_buffer_and_len
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
        character(*), intent(INOUT) :: arg1
        ! splicer begin accept_string_reference
        call str_accept_string_reference_bufferify(  &
            arg1,  &
            len_trim(arg1, kind=C_INT),  &
            len(arg1, kind=C_INT))
        ! splicer end accept_string_reference
    end subroutine accept_string_reference
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module strings_mod

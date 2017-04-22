! wrapfstrings.f
! This is generated code, do not edit
!>
!! \file wrapfstrings.f
!! \brief Shroud generated wrapper for strings library
!<
! splicer begin file_top
! splicer end file_top
module strings_mod
    ! splicer begin module_use
    ! splicer end module_use
    implicit none


    interface

        subroutine pass_char(status) &
                bind(C, name="STR_pass_char")
            use iso_c_binding, only : C_CHAR
            implicit none
            character(kind=C_CHAR), value, intent(IN) :: status
        end subroutine pass_char

        function c_return_char() &
                result(SH_rv) &
                bind(C, name="STR_return_char")
            use iso_c_binding, only : C_CHAR
            implicit none
            character(kind=C_CHAR) :: SH_rv
        end function c_return_char

        subroutine c_return_char_bufferify(SH_F_rv) &
                bind(C, name="STR_return_char_bufferify")
            use iso_c_binding, only : C_CHAR
            implicit none
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv
        end subroutine c_return_char_bufferify

        subroutine c_pass_char_ptr(dest, src) &
                bind(C, name="STR_pass_char_ptr")
            use iso_c_binding, only : C_CHAR
            implicit none
            character(kind=C_CHAR), intent(OUT) :: dest(*)
            character(kind=C_CHAR), intent(IN) :: src(*)
        end subroutine c_pass_char_ptr

        subroutine c_pass_char_ptr_bufferify(dest, Ndest, src, Lsrc) &
                bind(C, name="STR_pass_char_ptr_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(OUT) :: dest(*)
            integer(C_INT), value, intent(IN) :: Ndest
            character(kind=C_CHAR), intent(IN) :: src(*)
            integer(C_INT), value, intent(IN) :: Lsrc
        end subroutine c_pass_char_ptr_bufferify

        pure function c_get_char1() &
                result(SH_rv) &
                bind(C, name="STR_get_char1")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR) SH_rv
        end function c_get_char1

        subroutine c_get_char1_bufferify(SH_F_rv, NSH_F_rv) &
                bind(C, name="STR_get_char1_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv(*)
            integer(C_INT), value, intent(IN) :: NSH_F_rv
        end subroutine c_get_char1_bufferify

        function c_get_char2() &
                result(SH_rv) &
                bind(C, name="STR_get_char2")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR) SH_rv
        end function c_get_char2

        subroutine c_get_char2_bufferify(SH_F_rv, NSH_F_rv) &
                bind(C, name="STR_get_char2_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv(*)
            integer(C_INT), value, intent(IN) :: NSH_F_rv
        end subroutine c_get_char2_bufferify

        function c_get_char3() &
                result(SH_rv) &
                bind(C, name="STR_get_char3")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR) SH_rv
        end function c_get_char3

        subroutine c_get_char3_bufferify(output, Noutput) &
                bind(C, name="STR_get_char3_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(OUT) :: output(*)
            integer(C_INT), value, intent(IN) :: Noutput
        end subroutine c_get_char3_bufferify

        pure function c_get_string1() &
                result(SH_rv) &
                bind(C, name="STR_get_string1")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR) SH_rv
        end function c_get_string1

        subroutine c_get_string1_bufferify(SH_F_rv, NSH_F_rv) &
                bind(C, name="STR_get_string1_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv(*)
            integer(C_INT), value, intent(IN) :: NSH_F_rv
        end subroutine c_get_string1_bufferify

        function c_get_string2() &
                result(SH_rv) &
                bind(C, name="STR_get_string2")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR) SH_rv
        end function c_get_string2

        subroutine c_get_string2_bufferify(SH_F_rv, NSH_F_rv) &
                bind(C, name="STR_get_string2_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv(*)
            integer(C_INT), value, intent(IN) :: NSH_F_rv
        end subroutine c_get_string2_bufferify

        function c_get_string3() &
                result(SH_rv) &
                bind(C, name="STR_get_string3")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR) SH_rv
        end function c_get_string3

        subroutine c_get_string3_bufferify(output, Noutput) &
                bind(C, name="STR_get_string3_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(OUT) :: output(*)
            integer(C_INT), value, intent(IN) :: Noutput
        end subroutine c_get_string3_bufferify

        subroutine c_accept_string_const_reference(arg1) &
                bind(C, name="STR_accept_string_const_reference")
            use iso_c_binding, only : C_CHAR
            implicit none
            character(kind=C_CHAR), intent(IN) :: arg1(*)
        end subroutine c_accept_string_const_reference

        subroutine c_accept_string_const_reference_bufferify(arg1, Larg1) &
                bind(C, name="STR_accept_string_const_reference_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: arg1(*)
            integer(C_INT), value, intent(IN) :: Larg1
        end subroutine c_accept_string_const_reference_bufferify

        subroutine c_accept_string_reference(arg1) &
                bind(C, name="STR_accept_string_reference")
            use iso_c_binding, only : C_CHAR
            implicit none
            character(kind=C_CHAR), intent(INOUT) :: arg1(*)
        end subroutine c_accept_string_reference

        subroutine c_accept_string_reference_bufferify(arg1, Larg1, Narg1) &
                bind(C, name="STR_accept_string_reference_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(INOUT) :: arg1(*)
            integer(C_INT), value, intent(IN) :: Larg1
            integer(C_INT), value, intent(IN) :: Narg1
        end subroutine c_accept_string_reference_bufferify

        subroutine c_explicit1(name) &
                bind(C, name="STR_explicit1")
            use iso_c_binding, only : C_CHAR
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
        end subroutine c_explicit1

        subroutine c_explicit1_buffer(name, AAlen) &
                bind(C, name="STR_explicit1_BUFFER")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: AAlen
        end subroutine c_explicit1_buffer

        subroutine c_explicit2(name) &
                bind(C, name="STR_explicit2")
            use iso_c_binding, only : C_CHAR
            implicit none
            character(kind=C_CHAR), intent(OUT) :: name(*)
        end subroutine c_explicit2

        subroutine c_explicit2_bufferify(name, AAtrim) &
                bind(C, name="STR_explicit2_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(OUT) :: name(*)
            integer(C_INT), value, intent(IN) :: AAtrim
        end subroutine c_explicit2_bufferify

        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

    private fstr, fstr_ptr, fstr_arr, strlen_arr, strlen_ptr

    interface fstr
      module procedure fstr_ptr, fstr_arr
    end interface

    interface
       pure function strlen_ptr(s) result(result) bind(c,name="strlen")
         use, intrinsic :: iso_c_binding
         integer(c_int) :: result
         type(c_ptr), value, intent(in) :: s
       end function strlen_ptr
    end interface

contains

    ! char_scalar returnChar()
    ! string_to_buffer_and_len
    ! function_index=1
    !>
    !! \brief return a char argument (non-pointer)
    !!
    !<
    function return_char() result(SH_rv)
        character :: SH_rv
        ! splicer begin return_char
        call c_return_char_bufferify(SH_rv)
        ! splicer end return_char
    end function return_char

    ! void passCharPtr(char * dest+intent(out), const char * src+intent(in))
    ! string_to_buffer_and_len
    ! function_index=2
    !>
    !! \brief strcpy like behavior
    !!
    !! dest is marked intent(OUT) to override the intent(INOUT) default
    !! This avoid a copy-in on dest.
    !<
    subroutine pass_char_ptr(dest, src)
        use iso_c_binding, only : C_INT
        character(*), intent(OUT) :: dest
        character(*), intent(IN) :: src
        ! splicer begin pass_char_ptr
        call c_pass_char_ptr_bufferify(  &
            dest,  &
            len(dest, kind=C_INT),  &
            src,  &
            len_trim(src, kind=C_INT))
        ! splicer end pass_char_ptr
    end subroutine pass_char_ptr

    ! const string_result_fstr * getChar1()+pure
    ! function_index=16
    !>
    !! \brief return a 'const char *' as character(*)
    !!
    !<
    function get_char1() result(SH_rv)
        use iso_c_binding, only : C_CHAR
        character(kind=C_CHAR, len=strlen_ptr(c_get_char1())) :: SH_rv
        ! splicer begin get_char1
        SH_rv = fstr(c_get_char1())
        ! splicer end get_char1
    end function get_char1

    ! const char * getChar2()
    ! string_to_buffer_and_len
    ! function_index=4
    !>
    !! \brief return 'const char *' with fixed size (len=30)
    !!
    !<
    function get_char2() result(SH_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        character(kind=C_CHAR, len=30) :: SH_rv
        ! splicer begin get_char2
        call c_get_char2_bufferify(  &
            SH_rv,  &
            len(SH_rv, kind=C_INT))
        ! splicer end get_char2
    end function get_char2

    ! void getChar3(char * output+intent(out)+len(Noutput))
    ! string_to_buffer_and_len - string_to_buffer_and_len
    ! function_index=19
    !>
    !! \brief return a 'const char *' as argument
    !!
    !<
    subroutine get_char3(output)
        use iso_c_binding, only : C_INT
        character(*), intent(OUT) :: output
        ! splicer begin get_char3
        call c_get_char3_bufferify(  &
            output,  &
            len(output, kind=C_INT))
        ! splicer end get_char3
    end subroutine get_char3

    ! const string_result_fstr & getString1()+pure
    ! function_index=21
    !>
    !! \brief return a 'const string&' as character(*)
    !!
    !<
    function get_string1() result(SH_rv)
        use iso_c_binding, only : C_CHAR
        character(kind=C_CHAR, len=strlen_ptr(c_get_string1())) :: SH_rv
        ! splicer begin get_string1
        SH_rv = fstr(c_get_string1())
        ! splicer end get_string1
    end function get_string1

    ! const string & getString2()
    ! string_to_buffer_and_len
    ! function_index=7
    !>
    !! \brief return 'const string&' with fixed size (len=30)
    !!
    !<
    function get_string2() result(SH_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        character(kind=C_CHAR, len=30) :: SH_rv
        ! splicer begin get_string2
        call c_get_string2_bufferify(  &
            SH_rv,  &
            len(SH_rv, kind=C_INT))
        ! splicer end get_string2
    end function get_string2

    ! void getString3(string & output+intent(out)+len(Noutput))
    ! string_to_buffer_and_len - string_to_buffer_and_len
    ! function_index=24
    !>
    !! \brief return a 'const string&' as argument
    !!
    !<
    subroutine get_string3(output)
        use iso_c_binding, only : C_INT
        character(*), intent(OUT) :: output
        ! splicer begin get_string3
        call c_get_string3_bufferify(  &
            output,  &
            len(output, kind=C_INT))
        ! splicer end get_string3
    end subroutine get_string3

    ! void acceptStringConstReference(const std::string & arg1+intent(in))
    ! string_to_buffer_and_len
    ! function_index=9
    !>
    !! \brief Accept a const string reference
    !!
    !! Save contents of arg1.
    !! arg1 is assumed to be intent(IN) since it is const
    !! Will copy in.
    !<
    subroutine accept_string_const_reference(arg1)
        use iso_c_binding, only : C_INT
        character(*), intent(IN) :: arg1
        ! splicer begin accept_string_const_reference
        call c_accept_string_const_reference_bufferify(  &
            arg1,  &
            len_trim(arg1, kind=C_INT))
        ! splicer end accept_string_const_reference
    end subroutine accept_string_const_reference

    ! void acceptStringReference(std::string & arg1+intent(inout))
    ! string_to_buffer_and_len
    ! function_index=10
    !>
    !! \brief Accept a string reference
    !!
    !! Append "dog" to the end of arg1.
    !! arg1 is assumed to be intent(INOUT)
    !! Must copy in and copy out.
    !<
    subroutine accept_string_reference(arg1)
        use iso_c_binding, only : C_INT
        character(*), intent(INOUT) :: arg1
        ! splicer begin accept_string_reference
        call c_accept_string_reference_bufferify(  &
            arg1,  &
            len_trim(arg1, kind=C_INT),  &
            len(arg1, kind=C_INT))
        ! splicer end accept_string_reference
    end subroutine accept_string_reference

    ! void explicit1(char * name+intent(in)+len_trim(AAlen))
    ! string_to_buffer_and_len
    ! function_index=11
    subroutine explicit1(name)
        use iso_c_binding, only : C_INT
        character(*), intent(IN) :: name
        ! splicer begin explicit1
        call c_explicit1_buffer(  &
            name,  &
            len_trim(name, kind=C_INT))
        ! splicer end explicit1
    end subroutine explicit1

    ! void explicit2(char * name+intent(out)+len(AAtrim))
    ! string_to_buffer_and_len
    ! function_index=12
    subroutine explicit2(name)
        use iso_c_binding, only : C_INT
        character(*), intent(OUT) :: name
        ! splicer begin explicit2
        call c_explicit2_bufferify(  &
            name,  &
            len(name, kind=C_INT))
        ! splicer end explicit2
    end subroutine explicit2

    ! splicer begin additional_functions
    ! splicer end additional_functions

    ! Convert a null-terminated C "char *" pointer to a Fortran string.
    function fstr_ptr(s) result(fs)
      use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_f_pointer
      type(c_ptr), intent(in) :: s
      character(kind=c_char, len=strlen_ptr(s)) :: fs
      character(kind=c_char), pointer :: cptr(:)
      integer :: i
      call c_f_pointer(s, cptr, [len(fs)])
      do i=1, len(fs)
         fs(i:i) = cptr(i)
      enddo
    end function fstr_ptr

    ! Convert a null-terminated array of characters to a Fortran string.
    function fstr_arr(s) result(fs)
      use, intrinsic :: iso_c_binding, only : c_char, c_null_char
      character(kind=c_char, len=1), intent(in) :: s(*)
      character(kind=c_char, len=strlen_arr(s)) :: fs
      integer :: i
      do i = 1, len(fs)
         fs(i:i) = s(i)
      enddo
    end function fstr_arr

    ! Count the characters in a null-terminated array.
    pure function strlen_arr(s)
      use, intrinsic :: iso_c_binding, only : c_char, c_null_char
      character(kind=c_char, len=1), intent(in) :: s(*)
      integer :: i, strlen_arr
      i=1
      do
         if (s(i) == c_null_char) exit
         i = i+1
      enddo
      strlen_arr = i-1
    end function strlen_arr

end module strings_mod

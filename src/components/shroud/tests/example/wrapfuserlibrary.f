! wrapfuserlibrary.f
! This is generated code, do not edit
! blah blah
! yada yada
!
module userlibrary_mod
    use fstr_mod
    ! splicer begin module_use
    ! splicer end module_use
    implicit none
    
    
    interface
        
        subroutine local_function1() &
                bind(C, name="AA_local_function1")
            use iso_c_binding
            implicit none
        end subroutine local_function1
        
        function aa_is_name_valid(name) result(rv) &
                bind(C, name="AA_is_name_valid")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            logical(C_BOOL) :: rv
        end function aa_is_name_valid
        
        function aa_is_name_valid_bufferify(name, Lname) result(rv) &
                bind(C, name="AA_is_name_valid_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            logical(C_BOOL) :: rv
        end function aa_is_name_valid_bufferify
        
        subroutine aa_test_names(name) &
                bind(C, name="AA_test_names")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
        end subroutine aa_test_names
        
        subroutine aa_test_names_bufferify(name, Lname) &
                bind(C, name="AA_test_names_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
        end subroutine aa_test_names_bufferify
        
        subroutine aa_test_names_flag(name, flag) &
                bind(C, name="AA_test_names_flag")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: flag
        end subroutine aa_test_names_flag
        
        subroutine aa_test_names_flag_bufferify(name, Lname, flag) &
                bind(C, name="AA_test_names_flag_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            integer(C_INT), value, intent(IN) :: flag
        end subroutine aa_test_names_flag_bufferify
        
        subroutine aa_testoptional_0() &
                bind(C, name="AA_testoptional_0")
            use iso_c_binding
            implicit none
        end subroutine aa_testoptional_0
        
        subroutine aa_testoptional_1(i) &
                bind(C, name="AA_testoptional_1")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: i
        end subroutine aa_testoptional_1
        
        subroutine aa_testoptional_2(i, j) &
                bind(C, name="AA_testoptional_2")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: i
            integer(C_LONG), value, intent(IN) :: j
        end subroutine aa_testoptional_2
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface
    
    interface test_names
        module procedure test_names
        module procedure test_names_flag
    end interface test_names
    
    interface testoptional
        module procedure testoptional_0
        module procedure testoptional_1
        module procedure testoptional_2
    end interface testoptional

contains
    
    ! bool isNameValid(const std::string & name+intent(in))
    ! string_to_buffer_and_len
    ! function_index=37
    function is_name_valid(name) result(rv)
        use iso_c_binding
        implicit none
        character(*), intent(IN) :: name
        logical :: rv
        ! splicer begin is_name_valid
        rv = name .ne. " "
        ! splicer end is_name_valid
    end function is_name_valid
    
    ! void test_names(const std::string & name+intent(in))
    ! string_to_buffer_and_len
    ! function_index=38
    subroutine test_names(name)
        use iso_c_binding
        implicit none
        character(*), intent(IN) :: name
        ! splicer begin test_names
        call aa_test_names_bufferify(  &
            name,  &
            len_trim(name))
        ! splicer end test_names
    end subroutine test_names
    
    ! void test_names(const std::string & name+intent(in), int flag+intent(in)+value)
    ! string_to_buffer_and_len
    ! function_index=39
    subroutine test_names_flag(name, flag)
        use iso_c_binding
        implicit none
        character(*), intent(IN) :: name
        integer(C_INT), value, intent(IN) :: flag
        ! splicer begin test_names_flag
        call aa_test_names_flag_bufferify(  &
            name,  &
            len_trim(name),  &
            flag)
        ! splicer end test_names_flag
    end subroutine test_names_flag
    
    ! void testoptional()
    ! has_default_arg
    ! function_index=41
    subroutine testoptional_0()
        use iso_c_binding
        implicit none
        ! splicer begin testoptional_0
        call aa_testoptional_0()
        ! splicer end testoptional_0
    end subroutine testoptional_0
    
    ! void testoptional(int i+default(1)+intent(in)+value)
    ! has_default_arg
    ! function_index=42
    subroutine testoptional_1(i)
        use iso_c_binding
        implicit none
        integer(C_INT), value, intent(IN) :: i
        ! splicer begin testoptional_1
        call aa_testoptional_1(i)
        ! splicer end testoptional_1
    end subroutine testoptional_1
    
    ! void testoptional(int i+default(1)+intent(in)+value, long j+default(2)+intent(in)+value)
    ! function_index=40
    subroutine testoptional_2(i, j)
        use iso_c_binding
        implicit none
        integer(C_INT), value, intent(IN) :: i
        integer(C_LONG), value, intent(IN) :: j
        ! splicer begin testoptional_2
        call aa_testoptional_2(  &
            i,  &
            j)
        ! splicer end testoptional_2
    end subroutine testoptional_2
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module userlibrary_mod

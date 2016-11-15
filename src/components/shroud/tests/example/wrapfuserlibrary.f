! wrapfuserlibrary.f
! This is generated code, do not edit
! blah blah
! yada yada
!
!>
!! \file wrapfuserlibrary.f
!! \brief Shroud generated wrapper for UserLibrary library
!<
module userlibrary_mod
    ! splicer begin module_use
    ! splicer end module_use
    implicit none
    
    
    interface
        
        subroutine local_function1() &
                bind(C, name="AA_local_function1")
            use iso_c_binding
            implicit none
        end subroutine local_function1
        
        function c_is_name_valid(name) &
                result(rv) &
                bind(C, name="AA_is_name_valid")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            logical(C_BOOL) :: rv
        end function c_is_name_valid
        
        function c_is_name_valid_bufferify(name, Lname) &
                result(rv) &
                bind(C, name="AA_is_name_valid_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            logical(C_BOOL) :: rv
        end function c_is_name_valid_bufferify
        
        function c_is_initialized() &
                result(rv) &
                bind(C, name="AA_is_initialized")
            use iso_c_binding
            implicit none
            logical(C_BOOL) :: rv
        end function c_is_initialized
        
        subroutine c_test_names(name) &
                bind(C, name="AA_test_names")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
        end subroutine c_test_names
        
        subroutine c_test_names_bufferify(name, Lname) &
                bind(C, name="AA_test_names_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
        end subroutine c_test_names_bufferify
        
        subroutine c_test_names_flag(name, flag) &
                bind(C, name="AA_test_names_flag")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: flag
        end subroutine c_test_names_flag
        
        subroutine c_test_names_flag_bufferify(name, Lname, flag) &
                bind(C, name="AA_test_names_flag_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            integer(C_INT), value, intent(IN) :: flag
        end subroutine c_test_names_flag_bufferify
        
        subroutine c_testoptional_0() &
                bind(C, name="AA_testoptional_0")
            use iso_c_binding
            implicit none
        end subroutine c_testoptional_0
        
        subroutine c_testoptional_1(i) &
                bind(C, name="AA_testoptional_1")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: i
        end subroutine c_testoptional_1
        
        subroutine c_testoptional_2(i, j) &
                bind(C, name="AA_testoptional_2")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: i
            integer(C_LONG), value, intent(IN) :: j
        end subroutine c_testoptional_2
        
        subroutine testmpi(comm) &
                bind(C, name="AA_testmpi")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: comm
        end subroutine testmpi
        
        subroutine c_testgroup1(grp) &
                bind(C, name="AA_testgroup1")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: grp
        end subroutine c_testgroup1
        
        subroutine c_testgroup2(grp) &
                bind(C, name="AA_testgroup2")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: grp
        end subroutine c_testgroup2
        
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
    ! function_index=48
    function is_name_valid(name) result(rv)
        use iso_c_binding, only : C_BOOL, C_INT
        implicit none
        character(*), intent(IN) :: name
        logical :: rv
        ! splicer begin is_name_valid
        rv = name .ne. " "
        ! splicer end is_name_valid
    end function is_name_valid
    
    ! bool isInitialized()
    ! function_index=49
    function is_initialized() result(rv)
        use iso_c_binding, only : C_BOOL
        implicit none
        logical :: rv
        ! splicer begin is_initialized
        rv = c_is_initialized()
        ! splicer end is_initialized
    end function is_initialized
    
    ! void test_names(const std::string & name+intent(in))
    ! string_to_buffer_and_len
    ! function_index=50
    subroutine test_names(name)
        use iso_c_binding, only : C_INT
        implicit none
        character(*), intent(IN) :: name
        ! splicer begin test_names
        call c_test_names_bufferify(  &
            name,  &
            len_trim(name, kind=C_INT))
        ! splicer end test_names
    end subroutine test_names
    
    ! void test_names(const std::string & name+intent(in), int flag+intent(in)+value)
    ! string_to_buffer_and_len
    ! function_index=51
    subroutine test_names_flag(name, flag)
        use iso_c_binding, only : C_INT
        implicit none
        character(*), intent(IN) :: name
        integer(C_INT), value, intent(IN) :: flag
        ! splicer begin test_names_flag
        call c_test_names_flag_bufferify(  &
            name,  &
            len_trim(name, kind=C_INT),  &
            flag)
        ! splicer end test_names_flag
    end subroutine test_names_flag
    
    ! void testoptional()
    ! has_default_arg
    ! function_index=56
    subroutine testoptional_0()
        implicit none
        ! splicer begin testoptional_0
        call c_testoptional_0()
        ! splicer end testoptional_0
    end subroutine testoptional_0
    
    ! void testoptional(int i+default(1)+intent(in)+value)
    ! has_default_arg
    ! function_index=57
    subroutine testoptional_1(i)
        use iso_c_binding, only : C_INT
        implicit none
        integer(C_INT), value, intent(IN) :: i
        ! splicer begin testoptional_1
        call c_testoptional_1(i)
        ! splicer end testoptional_1
    end subroutine testoptional_1
    
    ! void testoptional(int i+default(1)+intent(in)+value, long j+default(2)+intent(in)+value)
    ! function_index=52
    subroutine testoptional_2(i, j)
        use iso_c_binding, only : C_LONG, C_INT
        implicit none
        integer(C_INT), value, intent(IN) :: i
        integer(C_LONG), value, intent(IN) :: j
        ! splicer begin testoptional_2
        call c_testoptional_2(  &
            i,  &
            j)
        ! splicer end testoptional_2
    end subroutine testoptional_2
    
    ! void testgroup1(DataGroup * grp+intent(in)+value)
    ! function_index=54
    subroutine testgroup1(grp)
        use sidre_mod, only : datagroup
        implicit none
        type(datagroup), value, intent(IN) :: grp
        ! splicer begin testgroup1
        call c_testgroup1(grp%get_instance())
        ! splicer end testgroup1
    end subroutine testgroup1
    
    ! void testgroup2(const DataGroup * grp+intent(in)+value)
    ! function_index=55
    subroutine testgroup2(grp)
        use sidre_mod, only : datagroup
        implicit none
        type(datagroup), value, intent(IN) :: grp
        ! splicer begin testgroup2
        call c_testgroup2(grp%get_instance())
        ! splicer end testgroup2
    end subroutine testgroup2
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module userlibrary_mod

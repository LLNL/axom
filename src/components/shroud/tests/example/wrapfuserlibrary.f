! wrapfuserlibrary.f
! This is generated code, do not edit
! blah blah
! yada yada
!
module userlibrary_mod
    use fstr_mod
    implicit none
    
    
    interface
        
        subroutine aa_local_function1() &
                bind(C, name="AA_local_function1")
            use iso_c_binding
            implicit none
        end subroutine aa_local_function1
        
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
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface
    
    interface test_names
        module procedure test_names
        module procedure test_names_flag
    end interface test_names

contains
    
    subroutine local_function1()
        use iso_c_binding
        implicit none
        ! splicer begin local_function1
        call aa_local_function1()
        ! splicer end local_function1
    end subroutine local_function1
    
    function is_name_valid(name) result(rv)
        use iso_c_binding
        implicit none
        character(*), intent(IN) :: name
        logical :: rv
        ! splicer begin is_name_valid
        rv = name .ne. " "
        ! splicer end is_name_valid
    end function is_name_valid
    
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
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module userlibrary_mod

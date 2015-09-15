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
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

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
        character(*) :: name
        logical :: rv
        ! splicer begin is_name_valid
        rv = name .ne. " "
        ! splicer end is_name_valid
    end function is_name_valid
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module userlibrary_mod

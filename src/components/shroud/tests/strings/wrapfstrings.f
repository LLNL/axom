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
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains
    
    ! const string & getName() const
    ! function_index=0
    !>
    !! \brief return a string as argument
    !!
    !<
    subroutine get_name(output)
        use iso_c_binding
        implicit none
        character(*), intent(OUT) :: output
        type(C_PTR) :: rv
        ! splicer begin get_name
        rv = str_get_name()
        call FccCopyPtr(output, len(output), rv)
        ! splicer end get_name
    end subroutine get_name
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module strings_mod

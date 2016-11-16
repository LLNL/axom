! wrapfdefault_library.f
! This is generated code, do not edit
!>
!! \file wrapfdefault_library.f
!! \brief Shroud generated wrapper for default_library library
!<
module default_library_mod
    ! splicer begin module_use
    ! splicer end module_use
    implicit none
    
    
    interface
        
        subroutine function1() &
                bind(C, name="DEF_function1")
            use iso_c_binding
            implicit none
        end subroutine function1
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module default_library_mod

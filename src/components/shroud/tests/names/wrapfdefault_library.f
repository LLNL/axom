! wrapfdefault_library.f
! This is generated code, do not edit
module default_library_mod
    use fstr_mod
    implicit none
    
    
    interface
        
        subroutine yyy_def_function1() &
                bind(C, name="YYY_DEF_function1")
            use iso_c_binding
            implicit none
        end subroutine yyy_def_function1
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

contains
    
    subroutine function1()
        use iso_c_binding
        implicit none
        ! splicer begin function1
        call yyy_def_function1()
        ! splicer end function1
    end subroutine function1
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module default_library_mod

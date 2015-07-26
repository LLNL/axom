! wrapftutorial.f
! This is generated code, do not edit
module tutorial_mod
    use fstr_mod
    implicit none
    
    
    interface
        
        subroutine tut_function1() &
                bind(C, name="TUT_function1")
            use iso_c_binding
            implicit none
        end subroutine tut_function1
    end interface

contains
    
    subroutine function1()
        use iso_c_binding
        implicit none
        ! splicer begin function1
        call tut_function1()
        ! splicer end function1
    end subroutine function1

end module tutorial_mod

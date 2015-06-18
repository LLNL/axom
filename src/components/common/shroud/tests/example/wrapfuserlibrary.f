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
    end interface

contains
    
    subroutine local_function1()
        use iso_c_binding
        implicit none
        ! splicer begin local_function1
        call aa_local_function1()
        ! splicer end local_function1
    end subroutine local_function1

end module userlibrary_mod

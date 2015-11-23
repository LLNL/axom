! top.f
! This is generated code, do not edit
module top_module
    use fstr_mod
    implicit none
    
    
    interface
        
        subroutine yyy_tes_function1() &
                bind(C, name="YYY_TES_function1")
            use iso_c_binding
            implicit none
        end subroutine yyy_tes_function1
        
        subroutine f_c_name_special() &
                bind(C, name="c_name_special")
            use iso_c_binding
            implicit none
        end subroutine f_c_name_special
        
        subroutine yyy_tes_function3a_0(i) &
                bind(C, name="YYY_TES_function3a_0")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: i
        end subroutine yyy_tes_function3a_0
        
        subroutine yyy_tes_function3a_1(i) &
                bind(C, name="YYY_TES_function3a_1")
            use iso_c_binding
            implicit none
            integer(C_LONG), value, intent(IN) :: i
        end subroutine yyy_tes_function3a_1
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface
    
    interface generic3
        module procedure F_name_function3a_int
        module procedure F_name_function3a_long
    end interface generic3

contains
    
    subroutine testnames_function1()
        use iso_c_binding
        implicit none
        ! splicer begin function1
        call yyy_tes_function1()
        ! splicer end function1
    end subroutine testnames_function1
    
    subroutine f_name_special()
        use iso_c_binding
        implicit none
        ! splicer begin function2
        call f_c_name_special()
        ! splicer end function2
    end subroutine f_name_special
    
    subroutine F_name_function3a_int(i)
        use iso_c_binding
        implicit none
        integer(C_INT) :: i
        ! splicer begin function3a_0
        call yyy_tes_function3a_0(i)
        ! splicer end function3a_0
    end subroutine F_name_function3a_int
    
    subroutine F_name_function3a_long(i)
        use iso_c_binding
        implicit none
        integer(C_LONG) :: i
        ! splicer begin function3a_1
        call yyy_tes_function3a_1(i)
        ! splicer end function3a_1
    end subroutine F_name_function3a_long
    
    ! splicer begin additional_functions
    ! splicer end additional_functions

end module top_module

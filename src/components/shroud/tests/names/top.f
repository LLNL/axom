! top.f
! This is generated code, do not edit
!>
!! \file top.f
!! \brief Shroud generated wrapper for testnames library
!<
! splicer begin file_top
! splicer end file_top
module top_module
    ! splicer begin module_use
    ! splicer end module_use
    implicit none


    interface

        subroutine yyy_tes_function1() &
                bind(C, name="YYY_TES_function1")
            implicit none
        end subroutine yyy_tes_function1

        subroutine f_c_name_special() &
                bind(C, name="c_name_special")
            implicit none
        end subroutine f_c_name_special

        subroutine yyy_tes_function3a_0(i) &
                bind(C, name="YYY_TES_function3a_0")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value, intent(IN) :: i
        end subroutine yyy_tes_function3a_0

        subroutine yyy_tes_function3a_1(i) &
                bind(C, name="YYY_TES_function3a_1")
            use iso_c_binding, only : C_LONG
            implicit none
            integer(C_LONG), value, intent(IN) :: i
        end subroutine yyy_tes_function3a_1

        function yyy_tes_function4() &
                result(RV) &
                bind(C, name="YYY_TES_function4")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT) :: RV
        end function yyy_tes_function4

        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
    end interface

    interface generic3
        module procedure F_name_function3a_int
        module procedure F_name_function3a_long
    end interface generic3

contains

    ! void function1()
    ! function_index=2
    subroutine testnames_function1()
        ! splicer begin function.function1
        call yyy_tes_function1()
        ! splicer end function.function1
    end subroutine testnames_function1

    ! void function2()
    ! function_index=3
    subroutine f_name_special()
        ! splicer begin function.function2
        call f_c_name_special()
        ! splicer end function.function2
    end subroutine f_name_special

    ! void function3a(int i+intent(in)+value)
    ! function_index=4
    subroutine F_name_function3a_int(i)
        use iso_c_binding, only : C_INT
        integer(C_INT), value, intent(IN) :: i
        ! splicer begin function.function3a_0
        call yyy_tes_function3a_0(i)
        ! splicer end function.function3a_0
    end subroutine F_name_function3a_int

    ! void function3a(long i+intent(in)+value)
    ! function_index=5
    subroutine F_name_function3a_long(i)
        use iso_c_binding, only : C_LONG
        integer(C_LONG), value, intent(IN) :: i
        ! splicer begin function.function3a_1
        call yyy_tes_function3a_1(i)
        ! splicer end function.function3a_1
    end subroutine F_name_function3a_long

    ! int function4()
    ! function_index=6
    function testnames_function4() result(RV)
        use iso_c_binding, only : C_INT
        integer(C_INT) :: RV
        ! splicer begin function.function4
        RV = yyy_tes_function4()
        ! splicer end function.function4
    end function testnames_function4

    ! splicer begin additional_functions
    ! splicer end additional_functions

end module top_module

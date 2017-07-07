! wrapfuserlibrary.f
! This is generated code, do not edit
! blah blah
! yada yada
!
!>
!! \file wrapfuserlibrary.f
!! \brief Shroud generated wrapper for UserLibrary library
!<
! splicer begin file_top
! splicer end file_top
module userlibrary_mod
    ! splicer begin module_use
    ! splicer end module_use
    implicit none


    interface

        subroutine local_function1() &
                bind(C, name="AA_local_function1")
            implicit none
        end subroutine local_function1

        function c_is_name_valid(name) &
                result(SH_rv) &
                bind(C, name="AA_is_name_valid")
            use iso_c_binding, only : C_BOOL, C_CHAR
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            logical(C_BOOL) :: SH_rv
        end function c_is_name_valid

        function c_is_name_valid_bufferify(name, Lname) &
                result(SH_rv) &
                bind(C, name="AA_is_name_valid_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            logical(C_BOOL) :: SH_rv
        end function c_is_name_valid_bufferify

        function c_is_initialized() &
                result(SH_rv) &
                bind(C, name="AA_is_initialized")
            use iso_c_binding, only : C_BOOL
            implicit none
            logical(C_BOOL) :: SH_rv
        end function c_is_initialized

        subroutine c_check_bool(arg1, arg2, arg3) &
                bind(C, name="AA_check_bool")
            use iso_c_binding, only : C_BOOL
            implicit none
            logical(C_BOOL), value, intent(IN) :: arg1
            logical(C_BOOL), intent(OUT) :: arg2
            logical(C_BOOL), intent(INOUT) :: arg3
        end subroutine c_check_bool

        subroutine c_test_names(name) &
                bind(C, name="AA_test_names")
            use iso_c_binding, only : C_CHAR
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
        end subroutine c_test_names

        subroutine c_test_names_bufferify(name, Lname) &
                bind(C, name="AA_test_names_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
        end subroutine c_test_names_bufferify

        subroutine c_test_names_flag(name, flag) &
                bind(C, name="AA_test_names_flag")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: flag
        end subroutine c_test_names_flag

        subroutine c_test_names_flag_bufferify(name, Lname, flag) &
                bind(C, name="AA_test_names_flag_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            integer(C_INT), value, intent(IN) :: flag
        end subroutine c_test_names_flag_bufferify

        subroutine c_testoptional_0() &
                bind(C, name="AA_testoptional_0")
            implicit none
        end subroutine c_testoptional_0

        subroutine c_testoptional_1(i) &
                bind(C, name="AA_testoptional_1")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value, intent(IN) :: i
        end subroutine c_testoptional_1

        subroutine c_testoptional_2(i, j) &
                bind(C, name="AA_testoptional_2")
            use iso_c_binding, only : C_INT, C_LONG
            implicit none
            integer(C_INT), value, intent(IN) :: i
            integer(C_LONG), value, intent(IN) :: j
        end subroutine c_testoptional_2

        subroutine testmpi(comm) &
                bind(C, name="AA_testmpi")
            use iso_c_binding, only : C_INT
            implicit none
            integer(C_INT), value, intent(IN) :: comm
        end subroutine testmpi

        subroutine c_testgroup1(grp) &
                bind(C, name="AA_testgroup1")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: grp
        end subroutine c_testgroup1

        subroutine c_testgroup2(grp) &
                bind(C, name="AA_testgroup2")
            use iso_c_binding, only : C_PTR
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
    function is_name_valid(name) result(SH_rv)
        use iso_c_binding, only : C_BOOL, C_INT
        character(*), intent(IN) :: name
        logical :: SH_rv
        ! splicer begin function.is_name_valid
        rv = name .ne. " "
        ! splicer end function.is_name_valid
    end function is_name_valid

    ! bool isInitialized()
    ! function_index=49
    function is_initialized() result(SH_rv)
        use iso_c_binding, only : C_BOOL
        logical :: SH_rv
        ! splicer begin function.is_initialized
        SH_rv = c_is_initialized()
        ! splicer end function.is_initialized
    end function is_initialized

    ! void checkBool(bool arg1+intent(in)+value, bool * arg2+intent(out), bool * arg3+intent(inout))
    ! function_index=50
    subroutine check_bool(arg1, arg2, arg3)
        use iso_c_binding, only : C_BOOL
        logical, value, intent(IN) :: arg1
        logical(C_BOOL) SH_arg1
        logical, intent(OUT) :: arg2
        logical(C_BOOL) SH_arg2
        logical, intent(INOUT) :: arg3
        logical(C_BOOL) SH_arg3
        SH_arg1 = arg1  ! coerce to C_BOOL
        SH_arg3 = arg3  ! coerce to C_BOOL
        ! splicer begin function.check_bool
        call c_check_bool(  &
            SH_arg1,  &
            SH_arg2,  &
            SH_arg3)
        ! splicer end function.check_bool
        arg2 = SH_arg2  ! coerce to logical
        arg3 = SH_arg3  ! coerce to logical
    end subroutine check_bool

    ! void test_names(const std::string & name+intent(in))
    ! string_to_buffer_and_len
    ! function_index=51
    subroutine test_names(name)
        use iso_c_binding, only : C_INT
        character(*), intent(IN) :: name
        ! splicer begin function.test_names
        call c_test_names_bufferify(  &
            name,  &
            len_trim(name, kind=C_INT))
        ! splicer end function.test_names
    end subroutine test_names

    ! void test_names(const std::string & name+intent(in), int flag+intent(in)+value)
    ! string_to_buffer_and_len
    ! function_index=52
    subroutine test_names_flag(name, flag)
        use iso_c_binding, only : C_INT
        character(*), intent(IN) :: name
        integer(C_INT), value, intent(IN) :: flag
        ! splicer begin function.test_names_flag
        call c_test_names_flag_bufferify(  &
            name,  &
            len_trim(name, kind=C_INT),  &
            flag)
        ! splicer end function.test_names_flag
    end subroutine test_names_flag

    ! void testoptional()
    ! has_default_arg
    ! function_index=57
    subroutine testoptional_0()
        ! splicer begin function.testoptional_0
        call c_testoptional_0()
        ! splicer end function.testoptional_0
    end subroutine testoptional_0

    ! void testoptional(int i+default(1)+intent(in)+value)
    ! has_default_arg
    ! function_index=58
    subroutine testoptional_1(i)
        use iso_c_binding, only : C_INT
        integer(C_INT), value, intent(IN) :: i
        ! splicer begin function.testoptional_1
        call c_testoptional_1(i)
        ! splicer end function.testoptional_1
    end subroutine testoptional_1

    ! void testoptional(int i+default(1)+intent(in)+value, long j+default(2)+intent(in)+value)
    ! function_index=53
    subroutine testoptional_2(i, j)
        use iso_c_binding, only : C_INT, C_LONG
        integer(C_INT), value, intent(IN) :: i
        integer(C_LONG), value, intent(IN) :: j
        ! splicer begin function.testoptional_2
        call c_testoptional_2(  &
            i,  &
            j)
        ! splicer end function.testoptional_2
    end subroutine testoptional_2

    ! void testgroup1(DataGroup * grp+intent(in)+value)
    ! function_index=55
    subroutine testgroup1(grp)
        use sidre_mod, only : group
        type(datagroup), value, intent(IN) :: grp
        ! splicer begin function.testgroup1
        call c_testgroup1(grp%get_instance())
        ! splicer end function.testgroup1
    end subroutine testgroup1

    ! void testgroup2(const DataGroup * grp+intent(in)+value)
    ! function_index=56
    subroutine testgroup2(grp)
        use sidre_mod, only : group
        type(datagroup), value, intent(IN) :: grp
        ! splicer begin function.testgroup2
        call c_testgroup2(grp%get_instance())
        ! splicer end function.testgroup2
    end subroutine testgroup2

    ! splicer begin additional_functions
    ! splicer end additional_functions

end module userlibrary_mod

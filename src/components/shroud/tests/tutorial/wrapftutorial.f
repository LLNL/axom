! wrapftutorial.f
! This is generated code, do not edit
module tutorial_mod
    use fstr_mod
    use, intrinsic :: iso_c_binding, only : C_PTR
    ! splicer begin module_use
    ! splicer end module_use
    ! splicer begin class.Class1.module_use
    ! splicer end class.Class1.module_use
    implicit none
    
    ! splicer begin module_top
    ! splicer end module_top
    
    ! splicer begin class.Class1.module_top
    ! splicer end class.Class1.module_top
    
    type class1
        type(C_PTR) voidptr
        ! splicer begin class.Class1.component_part
        ! splicer end class.Class1.component_part
    contains
        procedure :: delete => class1_delete
        procedure :: method1 => class1_method1
        ! splicer begin class.Class1.type_bound_procedure_part
        ! splicer end class.Class1.type_bound_procedure_part
    end type class1
    
    
    interface operator (.eq.)
        module procedure class1_eq
    end interface
    
    interface operator (.ne.)
        module procedure class1_ne
    end interface
    
    interface
        
        function tut_class1_new() result(rv) &
                bind(C, name="TUT_class1_new")
            use iso_c_binding
            implicit none
            type(C_PTR) :: rv
        end function tut_class1_new
        
        subroutine tut_class1_delete(self) &
                bind(C, name="TUT_class1_delete")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine tut_class1_delete
        
        subroutine tut_class1_method1(self) &
                bind(C, name="TUT_class1_method1")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine tut_class1_method1
        
        ! splicer begin class.Class1.additional_interfaces
        ! splicer end class.Class1.additional_interfaces
        
        subroutine tut_function1() &
                bind(C, name="TUT_function1")
            use iso_c_binding
            implicit none
        end subroutine tut_function1
        
        function tut_function2(arg1, arg2) result(rv) &
                bind(C, name="TUT_function2")
            use iso_c_binding
            implicit none
            real(C_DOUBLE), value, intent(IN) :: arg1
            integer(C_INT), value, intent(IN) :: arg2
            real(C_DOUBLE) :: rv
        end function tut_function2
        
        function tut_function3(arg) result(rv) &
                bind(C, name="TUT_function3")
            use iso_c_binding
            implicit none
            logical(C_BOOL), value, intent(IN) :: arg
            logical(C_BOOL) :: rv
        end function tut_function3
        
        pure function tut_function4a(arg1, arg2) result(rv) &
                bind(C, name="TUT_function4a")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: arg1(*)
            character(kind=C_CHAR), intent(IN) :: arg2(*)
            type(C_PTR) rv
        end function tut_function4a
        
        pure function tut_function4a_bufferify(arg1, Larg1, arg2, Larg2) result(rv) &
                bind(C, name="TUT_function4a_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: arg1(*)
            integer(C_INT), value, intent(IN) :: Larg1
            character(kind=C_CHAR), intent(IN) :: arg2(*)
            integer(C_INT), value, intent(IN) :: Larg2
            type(C_PTR) rv
        end function tut_function4a_bufferify
        
        function tut_function4b(arg1, arg2) result(rv) &
                bind(C, name="TUT_function4b")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: arg1(*)
            character(kind=C_CHAR), intent(IN) :: arg2(*)
            type(C_PTR) rv
        end function tut_function4b
        
        function tut_function4b_bufferify(arg1, Larg1, arg2, Larg2) result(rv) &
                bind(C, name="TUT_function4b_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: arg1(*)
            integer(C_INT), value, intent(IN) :: Larg1
            character(kind=C_CHAR), intent(IN) :: arg2(*)
            integer(C_INT), value, intent(IN) :: Larg2
            type(C_PTR) rv
        end function tut_function4b_bufferify
        
        function tut_function5(arg1, arg2) result(rv) &
                bind(C, name="TUT_function5")
            use iso_c_binding
            implicit none
            real(C_DOUBLE), value, intent(IN) :: arg1
            integer(C_INT), value, intent(IN) :: arg2
            real(C_DOUBLE) :: rv
        end function tut_function5
        
        subroutine tut_function6_from_name(name) &
                bind(C, name="TUT_function6_from_name")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
        end subroutine tut_function6_from_name
        
        subroutine tut_function6_bufferify(name, Lname) &
                bind(C, name="TUT_function6_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
        end subroutine tut_function6_bufferify
        
        subroutine tut_function6_from_index(indx) &
                bind(C, name="TUT_function6_from_index")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: indx
        end subroutine tut_function6_from_index
        
        subroutine tut_function7_int(arg) &
                bind(C, name="TUT_function7_int")
            use iso_c_binding
            implicit none
            integer(C_INT), value, intent(IN) :: arg
        end subroutine tut_function7_int
        
        subroutine tut_function7_double(arg) &
                bind(C, name="TUT_function7_double")
            use iso_c_binding
            implicit none
            real(C_DOUBLE), value, intent(IN) :: arg
        end subroutine tut_function7_double
        
        function tut_function8_int() result(rv) &
                bind(C, name="TUT_function8_int")
            use iso_c_binding
            implicit none
            integer(C_INT) :: rv
        end function tut_function8_int
        
        function tut_function8_double() result(rv) &
                bind(C, name="TUT_function8_double")
            use iso_c_binding
            implicit none
            real(C_DOUBLE) :: rv
        end function tut_function8_double
        
        subroutine tut_function9(arg) &
                bind(C, name="TUT_function9")
            use iso_c_binding
            implicit none
            real(C_DOUBLE), value, intent(IN) :: arg
        end subroutine tut_function9
        
        pure function tut_last_function_called() result(rv) &
                bind(C, name="TUT_last_function_called")
            use iso_c_binding
            implicit none
            type(C_PTR) rv
        end function tut_last_function_called
        
        ! splicer begin additional_interfaces
        subroutine all_test1(array)
          implicit none
          integer, dimension(:), allocatable :: array
        end subroutine all_test1
        ! splicer end additional_interfaces
    end interface
    
    interface function6
        module procedure function6_from_name
        module procedure function6_from_index
    end interface function6
    
    interface function7
        module procedure function7_int
        module procedure function7_double
    end interface function7
    
    interface function9
        module procedure function9_float
        module procedure function9_double
    end interface function9

contains
    
    function class1_new() result(rv)
        use iso_c_binding
        implicit none
        type(class1) :: rv
        ! splicer begin class.Class1.method.new
        rv%voidptr = tut_class1_new()
        ! splicer end class.Class1.method.new
    end function class1_new
    
    subroutine class1_delete(obj)
        use iso_c_binding
        implicit none
        class(class1) :: obj
        ! splicer begin class.Class1.method.delete
        call tut_class1_delete(obj%voidptr)
        obj%voidptr = C_NULL_PTR
        ! splicer end class.Class1.method.delete
    end subroutine class1_delete
    
    subroutine class1_method1(obj)
        use iso_c_binding
        implicit none
        class(class1) :: obj
        ! splicer begin class.Class1.method.method1
        call tut_class1_method1(obj%voidptr)
        ! splicer end class.Class1.method.method1
    end subroutine class1_method1
    
    ! splicer begin class.Class1.additional_functions
    ! splicer end class.Class1.additional_functions
    
    subroutine function1()
        use iso_c_binding
        implicit none
        ! splicer begin function1
        call tut_function1()
        ! splicer end function1
    end subroutine function1
    
    function function2(arg1, arg2) result(rv)
        use iso_c_binding
        implicit none
        real(C_DOUBLE) :: arg1
        integer(C_INT) :: arg2
        real(C_DOUBLE) :: rv
        ! splicer begin function2
        rv = tut_function2(  &
            arg1,  &
            arg2)
        ! splicer end function2
    end function function2
    
    function function3(arg) result(rv)
        use iso_c_binding
        implicit none
        logical :: arg
        logical :: rv
        logical(C_BOOL) tmp_arg
        tmp_arg = arg  ! coerce to C_BOOL
        ! splicer begin function3
        rv = tut_function3(tmp_arg)
        ! splicer end function3
    end function function3
    
    function function4a(arg1, arg2) result(rv)
        use iso_c_binding
        implicit none
        character(*) :: arg1
        character(*) :: arg2
        character(kind=C_CHAR, len=strlen_ptr(tut_function4a_bufferify(arg1, len_trim(arg1), arg2, len_trim(arg2)))) :: rv
        ! splicer begin function4a
        rv = fstr(tut_function4a_bufferify(  &
            arg1,  &
            len_trim(arg1),  &
            arg2,  &
            len_trim(arg2)))
        ! splicer end function4a
    end function function4a
    
    subroutine function4b(arg1, arg2, output)
        use iso_c_binding
        implicit none
        character(*) :: arg1
        character(*) :: arg2
        character(*), intent(OUT) :: output
        type(C_PTR) :: rv
        ! splicer begin function4b
        rv = tut_function4b_bufferify(  &
            arg1,  &
            len_trim(arg1),  &
            arg2,  &
            len_trim(arg2))
        call FccCopyPtr(output, len(output), rv)
        ! splicer end function4b
    end subroutine function4b
    
    function function5(arg1, arg2) result(rv)
        use iso_c_binding
        implicit none
        real(C_DOUBLE), optional :: arg1
        real(C_DOUBLE) :: tmp_arg1
        integer(C_INT), optional :: arg2
        integer(C_INT) :: tmp_arg2
        real(C_DOUBLE) :: rv
        if (present(arg1)) then
            tmp_arg1 = arg1
        else
            tmp_arg1 = 3.13
        endif
        if (present(arg2)) then
            tmp_arg2 = arg2
        else
            tmp_arg2 = 5
        endif
        ! splicer begin function5
        rv = tut_function5(  &
            tmp_arg1,  &
            tmp_arg2)
        ! splicer end function5
    end function function5
    
    subroutine function6_from_name(name)
        use iso_c_binding
        implicit none
        character(*) :: name
        ! splicer begin function6_from_name
        call tut_function6_bufferify(  &
            name,  &
            len_trim(name))
        ! splicer end function6_from_name
    end subroutine function6_from_name
    
    subroutine function6_from_index(indx)
        use iso_c_binding
        implicit none
        integer(C_INT) :: indx
        ! splicer begin function6_from_index
        call tut_function6_from_index(indx)
        ! splicer end function6_from_index
    end subroutine function6_from_index
    
    subroutine function7_int(arg)
        use iso_c_binding
        implicit none
        integer(C_INT) :: arg
        ! splicer begin function7_int
        call tut_function7_int(arg)
        ! splicer end function7_int
    end subroutine function7_int
    
    subroutine function7_double(arg)
        use iso_c_binding
        implicit none
        real(C_DOUBLE) :: arg
        ! splicer begin function7_double
        call tut_function7_double(arg)
        ! splicer end function7_double
    end subroutine function7_double
    
    function function8_int() result(rv)
        use iso_c_binding
        implicit none
        integer(C_INT) :: rv
        ! splicer begin function8_int
        rv = tut_function8_int()
        ! splicer end function8_int
    end function function8_int
    
    function function8_double() result(rv)
        use iso_c_binding
        implicit none
        real(C_DOUBLE) :: rv
        ! splicer begin function8_double
        rv = tut_function8_double()
        ! splicer end function8_double
    end function function8_double
    
    subroutine function9_float(arg)
        use iso_c_binding
        implicit none
        real(C_FLOAT) :: arg
        ! splicer begin function9_float
        call tut_function9(real(arg, C_DOUBLE))
        ! splicer end function9_float
    end subroutine function9_float
    
    subroutine function9_double(arg)
        use iso_c_binding
        implicit none
        real(C_DOUBLE) :: arg
        ! splicer begin function9_double
        call tut_function9(arg)
        ! splicer end function9_double
    end subroutine function9_double
    
    function last_function_called() result(rv)
        use iso_c_binding
        implicit none
        character(kind=C_CHAR, len=strlen_ptr(tut_last_function_called())) :: rv
        ! splicer begin last_function_called
        rv = fstr(tut_last_function_called())
        ! splicer end last_function_called
    end function last_function_called
    
    ! splicer begin additional_functions
    ! splicer end additional_functions
    
    function class1_eq(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(class1), intent(IN) ::a,b
        logical :: rv
        if (c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function class1_eq
    
    function class1_ne(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(class1), intent(IN) ::a,b
        logical :: rv
        if (.not. c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function class1_ne

end module tutorial_mod

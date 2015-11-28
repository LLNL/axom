! wrapfExClass2.f
! This is generated code, do not edit
! blah blah
! yada yada
!
module exclass2_mod
    use fstr_mod
    use exclass1_mod, only : exclass1
    use iso_c_binding, only : C_DOUBLE, C_FLOAT, C_INT, C_LONG
    ! splicer begin class.ExClass2.module_use
    ! splicer end class.ExClass2.module_use
    implicit none
    
    
    ! splicer begin class.ExClass2.module_top
    top of module splicer  2
    ! splicer end class.ExClass2.module_top
    
    type exclass2
        type(C_PTR) voidptr
        ! splicer begin class.ExClass2.component_part
        ! splicer end class.ExClass2.component_part
    contains
        procedure :: delete => exclass2_delete
        procedure :: get_name => exclass2_get_name
        procedure :: get_name_length => exclass2_get_name_length
        procedure :: get_class1 => exclass2_get_class1
        procedure :: declare_int => exclass2_declare_int
        procedure :: declare_long => exclass2_declare_long
        procedure :: destroyall => exclass2_destroyall
        procedure :: get_type_id => exclass2_get_type_id
        procedure :: set_value_int => exclass2_set_value_int
        procedure :: set_value_long => exclass2_set_value_long
        procedure :: set_value_float => exclass2_set_value_float
        procedure :: set_value_double => exclass2_set_value_double
        procedure :: get_value_int => exclass2_get_value_int
        procedure :: get_value_double => exclass2_get_value_double
        generic :: declare => &
            ! splicer begin class.ExClass2.generic.declare
            ! splicer end class.ExClass2.generic.declare
            declare_int,  &
            declare_long
        generic :: set_value => &
            ! splicer begin class.ExClass2.generic.set_value
            ! splicer end class.ExClass2.generic.set_value
            set_value_int,  &
            set_value_long,  &
            set_value_float,  &
            set_value_double
        ! splicer begin class.ExClass2.type_bound_procedure_part
        ! splicer end class.ExClass2.type_bound_procedure_part
    end type exclass2
    
    
    interface operator (.eq.)
        module procedure exclass2_eq
    end interface
    
    interface operator (.ne.)
        module procedure exclass2_ne
    end interface
    
    interface
        
        function aa_exclass2_ex_class2(name) result(rv) &
                bind(C, name="AA_exclass2_ex_class2")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR) :: rv
        end function aa_exclass2_ex_class2
        
        function aa_exclass2_ex_class2_bufferify(name, Lname) result(rv) &
                bind(C, name="AA_exclass2_ex_class2_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            type(C_PTR) :: rv
        end function aa_exclass2_ex_class2_bufferify
        
        subroutine aa_exclass2_delete(self) &
                bind(C, name="AA_exclass2_delete")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine aa_exclass2_delete
        
        pure function aa_exclass2_get_name(self) result(rv) &
                bind(C, name="AA_exclass2_get_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) rv
        end function aa_exclass2_get_name
        
        function aa_exclass2_get_name_length(self) result(rv) &
                bind(C, name="AA_exclass2_get_name_length")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: rv
        end function aa_exclass2_get_name_length
        
        function aa_exclass2_get_class1(self, in) result(rv) &
                bind(C, name="AA_exclass2_get_class1")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), value, intent(IN) :: in
            type(C_PTR) :: rv
        end function aa_exclass2_get_class1
        
        subroutine aa_exclass2_declare(self, type, len) &
                bind(C, name="AA_exclass2_declare")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
        end subroutine aa_exclass2_declare
        
        subroutine aa_exclass2_destroyall(self) &
                bind(C, name="AA_exclass2_destroyall")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine aa_exclass2_destroyall
        
        pure function aa_exclass2_get_type_id(self) result(rv) &
                bind(C, name="AA_exclass2_get_type_id")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: rv
        end function aa_exclass2_get_type_id
        
        subroutine aa_exclass2_set_value_int(self, value) &
                bind(C, name="AA_exclass2_set_value_int")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: value
        end subroutine aa_exclass2_set_value_int
        
        subroutine aa_exclass2_set_value_long(self, value) &
                bind(C, name="AA_exclass2_set_value_long")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_LONG), value, intent(IN) :: value
        end subroutine aa_exclass2_set_value_long
        
        subroutine aa_exclass2_set_value_float(self, value) &
                bind(C, name="AA_exclass2_set_value_float")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_FLOAT), value, intent(IN) :: value
        end subroutine aa_exclass2_set_value_float
        
        subroutine aa_exclass2_set_value_double(self, value) &
                bind(C, name="AA_exclass2_set_value_double")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_DOUBLE), value, intent(IN) :: value
        end subroutine aa_exclass2_set_value_double
        
        function aa_exclass2_get_value_int(self) result(rv) &
                bind(C, name="AA_exclass2_get_value_int")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: rv
        end function aa_exclass2_get_value_int
        
        function aa_exclass2_get_value_double(self) result(rv) &
                bind(C, name="AA_exclass2_get_value_double")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_DOUBLE) :: rv
        end function aa_exclass2_get_value_double
        
        ! splicer begin class.ExClass2.additional_interfaces
        ! splicer end class.ExClass2.additional_interfaces
    end interface

contains
    
    ! ExClass2 *ExClass2 (const string *name) +constructor
    ! string_to_buffer_and_len
    ! function_index=14
    function exclass2_ex_class2(name) result(rv)
        use iso_c_binding
        implicit none
        character(*), intent(IN) :: name
        type(exclass2) :: rv
        ! splicer begin class.ExClass2.method.ex_class2
        rv%voidptr = aa_exclass2_ex_class2_bufferify(  &
            name,  &
            len_trim(name))
        ! splicer end class.ExClass2.method.ex_class2
    end function exclass2_ex_class2
    
    ! void delete() +destructor
    ! function_index=15
    subroutine exclass2_delete(obj)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        ! splicer begin class.ExClass2.method.delete
        call aa_exclass2_delete(obj%voidptr)
        obj%voidptr = C_NULL_PTR
        ! splicer end class.ExClass2.method.delete
    end subroutine exclass2_delete
    
    ! const string& getName const
    ! function_index=16
    function exclass2_get_name(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        character(kind=C_CHAR, len=aa_exclass2_get_name_length(obj%voidptr)) :: rv
        ! splicer begin class.ExClass2.method.get_name
        rv = fstr(aa_exclass2_get_name(obj%voidptr))
        ! splicer end class.ExClass2.method.get_name
    end function exclass2_get_name
    
    ! function_index=17
    function exclass2_get_name_length(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        integer(C_INT) :: rv
        ! splicer begin class.ExClass2.method.get_name_length
        rv = aa_exclass2_get_name_length(obj%voidptr)
        ! splicer end class.ExClass2.method.get_name_length
    end function exclass2_get_name_length
    
    ! ExClass1 *get_class1(const ExClass1 *in)
    ! function_index=18
    function exclass2_get_class1(obj, in) result(rv)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        type(exclass1), value, intent(IN) :: in
        type(exclass1) :: rv
        ! splicer begin class.ExClass2.method.get_class1
        rv%voidptr = aa_exclass2_get_class1(  &
            obj%voidptr,  &
            in%voidptr)
        ! splicer end class.ExClass2.method.get_class1
    end function exclass2_get_class1
    
    ! void* declare(TypeID type, SidreLength len = 1)
    ! fortran_generic
    ! function_index=31
    subroutine exclass2_declare_int(obj, type, len)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN), optional :: len
        integer(C_INT) :: tmp_len
        if (present(len)) then
            tmp_len = len
        else
            tmp_len = 1
        endif
        ! splicer begin class.ExClass2.method.declare_int
        call aa_exclass2_declare(  &
            obj%voidptr,  &
            type,  &
            int(tmp_len, C_LONG))
        ! splicer end class.ExClass2.method.declare_int
    end subroutine exclass2_declare_int
    
    ! void* declare(TypeID type, SidreLength len = 1)
    ! fortran_generic
    ! function_index=32
    subroutine exclass2_declare_long(obj, type, len)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN), optional :: len
        integer(C_LONG) :: tmp_len
        if (present(len)) then
            tmp_len = len
        else
            tmp_len = 1
        endif
        ! splicer begin class.ExClass2.method.declare_long
        call aa_exclass2_declare(  &
            obj%voidptr,  &
            type,  &
            int(tmp_len, C_LONG))
        ! splicer end class.ExClass2.method.declare_long
    end subroutine exclass2_declare_long
    
    ! void destroyall()
    ! function_index=20
    subroutine exclass2_destroyall(obj)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        ! splicer begin class.ExClass2.method.destroyall
        call aa_exclass2_destroyall(obj%voidptr)
        ! splicer end class.ExClass2.method.destroyall
    end subroutine exclass2_destroyall
    
    ! TypeID getTypeID() const
    ! function_index=21
    function exclass2_get_type_id(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        integer(C_INT) :: rv
        ! splicer begin class.ExClass2.method.get_type_id
        rv = aa_exclass2_get_type_id(obj%voidptr)
        ! splicer end class.ExClass2.method.get_type_id
    end function exclass2_get_type_id
    
    ! void setValue(ValueType value)
    ! cpp_template
    ! function_index=24
    subroutine exclass2_set_value_int(obj, value)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        integer(C_INT), value, intent(IN) :: value
        ! splicer begin class.ExClass2.method.set_value_int
        call aa_exclass2_set_value_int(  &
            obj%voidptr,  &
            value)
        ! splicer end class.ExClass2.method.set_value_int
    end subroutine exclass2_set_value_int
    
    ! void setValue(ValueType value)
    ! cpp_template
    ! function_index=25
    subroutine exclass2_set_value_long(obj, value)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        integer(C_LONG), value, intent(IN) :: value
        ! splicer begin class.ExClass2.method.set_value_long
        call aa_exclass2_set_value_long(  &
            obj%voidptr,  &
            value)
        ! splicer end class.ExClass2.method.set_value_long
    end subroutine exclass2_set_value_long
    
    ! void setValue(ValueType value)
    ! cpp_template
    ! function_index=26
    subroutine exclass2_set_value_float(obj, value)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        real(C_FLOAT), value, intent(IN) :: value
        ! splicer begin class.ExClass2.method.set_value_float
        call aa_exclass2_set_value_float(  &
            obj%voidptr,  &
            value)
        ! splicer end class.ExClass2.method.set_value_float
    end subroutine exclass2_set_value_float
    
    ! void setValue(ValueType value)
    ! cpp_template
    ! function_index=27
    subroutine exclass2_set_value_double(obj, value)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        real(C_DOUBLE), value, intent(IN) :: value
        ! splicer begin class.ExClass2.method.set_value_double
        call aa_exclass2_set_value_double(  &
            obj%voidptr,  &
            value)
        ! splicer end class.ExClass2.method.set_value_double
    end subroutine exclass2_set_value_double
    
    ! ValueType getValue()
    ! cpp_template
    ! function_index=28
    function exclass2_get_value_int(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        integer(C_INT) :: rv
        ! splicer begin class.ExClass2.method.get_value_int
        rv = aa_exclass2_get_value_int(obj%voidptr)
        ! splicer end class.ExClass2.method.get_value_int
    end function exclass2_get_value_int
    
    ! ValueType getValue()
    ! cpp_template
    ! function_index=29
    function exclass2_get_value_double(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        real(C_DOUBLE) :: rv
        ! splicer begin class.ExClass2.method.get_value_double
        rv = aa_exclass2_get_value_double(obj%voidptr)
        ! splicer end class.ExClass2.method.get_value_double
    end function exclass2_get_value_double
    
    ! splicer begin class.ExClass2.additional_functions
    ! splicer end class.ExClass2.additional_functions
    
    function exclass2_eq(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(exclass2), intent(IN) ::a,b
        logical :: rv
        if (c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function exclass2_eq
    
    function exclass2_ne(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(exclass2), intent(IN) ::a,b
        logical :: rv
        if (.not. c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function exclass2_ne

end module exclass2_mod

! wrapfExClass2.f
! This is generated code, do not edit
! blah blah
! yada yada
!
!>
!! \file wrapfExClass2.f
!! \brief Shroud generated wrapper for ExClass2 class
!<
! splicer begin file_top
! splicer end file_top
module exclass2_mod
    use exclass1_mod, only : exclass1
    use iso_c_binding, only : C_DOUBLE, C_FLOAT, C_INT, C_LONG, C_PTR
    ! splicer begin class.ExClass2.module_use
    ! splicer end class.ExClass2.module_use
    implicit none


    ! splicer begin class.ExClass2.module_top
    top of module splicer  2
    ! splicer end class.ExClass2.module_top

    type exclass2
        type(C_PTR), private :: voidptr
        ! splicer begin class.ExClass2.component_part
        ! splicer end class.ExClass2.component_part
    contains
        procedure :: delete => exclass2_delete
        procedure :: get_name => exclass2_get_name
        procedure :: get_name2 => exclass2_get_name2
        procedure :: get_name3 => exclass2_get_name3
        procedure :: get_name4 => exclass2_get_name4
        procedure :: get_name_length => exclass2_get_name_length
        procedure :: get_class1 => exclass2_get_class1
        procedure :: declare_0_int => exclass2_declare_0_int
        procedure :: declare_0_long => exclass2_declare_0_long
        procedure :: declare_1_int => exclass2_declare_1_int
        procedure :: declare_1_long => exclass2_declare_1_long
        procedure :: destroyall => exclass2_destroyall
        procedure :: get_type_id => exclass2_get_type_id
        procedure :: set_value_int => exclass2_set_value_int
        procedure :: set_value_long => exclass2_set_value_long
        procedure :: set_value_float => exclass2_set_value_float
        procedure :: set_value_double => exclass2_set_value_double
        procedure :: get_value_int => exclass2_get_value_int
        procedure :: get_value_double => exclass2_get_value_double
        procedure :: yadda => exclass2_yadda
        procedure :: associated => exclass2_associated
        generic :: declare => &
            ! splicer begin class.ExClass2.generic.declare
            ! splicer end class.ExClass2.generic.declare
            declare_0_int,  &
            declare_0_long,  &
            declare_1_int,  &
            declare_1_long
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

        function c_exclass2_ex_class2(name) &
                result(SH_rv) &
                bind(C, name="AA_exclass2_ex_class2")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR) :: SH_rv
        end function c_exclass2_ex_class2

        function c_exclass2_ex_class2_bufferify(name, Lname) &
                result(SH_rv) &
                bind(C, name="AA_exclass2_ex_class2_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            type(C_PTR) :: SH_rv
        end function c_exclass2_ex_class2_bufferify

        subroutine c_exclass2_delete(self) &
                bind(C, name="AA_exclass2_delete")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine c_exclass2_delete

        pure function c_exclass2_get_name(self) &
                result(SH_rv) &
                bind(C, name="AA_exclass2_get_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) SH_rv
        end function c_exclass2_get_name

        subroutine c_exclass2_get_name_bufferify(self, SH_F_rv, NSH_F_rv) &
                bind(C, name="AA_exclass2_get_name_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv(*)
            integer(C_INT), value, intent(IN) :: NSH_F_rv
        end subroutine c_exclass2_get_name_bufferify

        function c_exclass2_get_name2(self) &
                result(SH_rv) &
                bind(C, name="AA_exclass2_get_name2")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) SH_rv
        end function c_exclass2_get_name2

        subroutine c_exclass2_get_name2_bufferify(self, SH_F_rv, NSH_F_rv) &
                bind(C, name="AA_exclass2_get_name2_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv(*)
            integer(C_INT), value, intent(IN) :: NSH_F_rv
        end subroutine c_exclass2_get_name2_bufferify

        pure function c_exclass2_get_name3(self) &
                result(SH_rv) &
                bind(C, name="AA_exclass2_get_name3")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) SH_rv
        end function c_exclass2_get_name3

        subroutine c_exclass2_get_name3_bufferify(self, SH_F_rv, NSH_F_rv) &
                bind(C, name="AA_exclass2_get_name3_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv(*)
            integer(C_INT), value, intent(IN) :: NSH_F_rv
        end subroutine c_exclass2_get_name3_bufferify

        function c_exclass2_get_name4(self) &
                result(SH_rv) &
                bind(C, name="AA_exclass2_get_name4")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) SH_rv
        end function c_exclass2_get_name4

        subroutine c_exclass2_get_name4_bufferify(self, SH_F_rv, NSH_F_rv) &
                bind(C, name="AA_exclass2_get_name4_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(OUT) :: SH_F_rv(*)
            integer(C_INT), value, intent(IN) :: NSH_F_rv
        end subroutine c_exclass2_get_name4_bufferify

        function c_exclass2_get_name_length(self) &
                result(SH_rv) &
                bind(C, name="AA_exclass2_get_name_length")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: SH_rv
        end function c_exclass2_get_name_length

        function c_exclass2_get_class1(self, in) &
                result(SH_rv) &
                bind(C, name="AA_exclass2_get_class1")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), value, intent(IN) :: in
            type(C_PTR) :: SH_rv
        end function c_exclass2_get_class1

        subroutine c_exclass2_declare_0(self, type) &
                bind(C, name="AA_exclass2_declare_0")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
        end subroutine c_exclass2_declare_0

        subroutine c_exclass2_declare_1(self, type, len) &
                bind(C, name="AA_exclass2_declare_1")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
        end subroutine c_exclass2_declare_1

        subroutine c_exclass2_destroyall(self) &
                bind(C, name="AA_exclass2_destroyall")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine c_exclass2_destroyall

        pure function c_exclass2_get_type_id(self) &
                result(SH_rv) &
                bind(C, name="AA_exclass2_get_type_id")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: SH_rv
        end function c_exclass2_get_type_id

        subroutine c_exclass2_set_value_int(self, value) &
                bind(C, name="AA_exclass2_set_value_int")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: value
        end subroutine c_exclass2_set_value_int

        subroutine c_exclass2_set_value_long(self, value) &
                bind(C, name="AA_exclass2_set_value_long")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_LONG), value, intent(IN) :: value
        end subroutine c_exclass2_set_value_long

        subroutine c_exclass2_set_value_float(self, value) &
                bind(C, name="AA_exclass2_set_value_float")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_FLOAT), value, intent(IN) :: value
        end subroutine c_exclass2_set_value_float

        subroutine c_exclass2_set_value_double(self, value) &
                bind(C, name="AA_exclass2_set_value_double")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_DOUBLE), value, intent(IN) :: value
        end subroutine c_exclass2_set_value_double

        function c_exclass2_get_value_int(self) &
                result(SH_rv) &
                bind(C, name="AA_exclass2_get_value_int")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: SH_rv
        end function c_exclass2_get_value_int

        function c_exclass2_get_value_double(self) &
                result(SH_rv) &
                bind(C, name="AA_exclass2_get_value_double")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_DOUBLE) :: SH_rv
        end function c_exclass2_get_value_double

        ! splicer begin class.ExClass2.additional_interfaces
        ! splicer end class.ExClass2.additional_interfaces
    end interface

contains

    ! ExClass2 * ExClass2(const string * name+intent(in))+constructor
    ! string_to_buffer_and_len
    ! function_index=18
    function exclass2_ex_class2(name) result(SH_rv)
        use iso_c_binding, only : C_INT
        implicit none
        character(*), intent(IN) :: name
        type(exclass2) :: SH_rv
        ! splicer begin class.ExClass2.method.ex_class2
        SH_rv%voidptr = c_exclass2_ex_class2_bufferify(  &
            name,  &
            len_trim(name, kind=C_INT))
        ! splicer end class.ExClass2.method.ex_class2
    end function exclass2_ex_class2

    ! void delete()+destructor
    ! function_index=19
    subroutine exclass2_delete(obj)
        use iso_c_binding, only : C_NULL_PTR
        implicit none
        class(exclass2) :: obj
        ! splicer begin class.ExClass2.method.delete
        call c_exclass2_delete(obj%voidptr)
        obj%voidptr = C_NULL_PTR
        ! splicer end class.ExClass2.method.delete
    end subroutine exclass2_delete

    ! const string & getName() const
    ! string_to_buffer_and_len
    ! function_index=20
    function exclass2_get_name(obj) result(SH_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        implicit none
        class(exclass2) :: obj
        character(kind=C_CHAR, len=aa_exclass2_get_name_length(obj%voidptr)) :: SH_rv
        ! splicer begin class.ExClass2.method.get_name
        call c_exclass2_get_name_bufferify(  &
            obj%voidptr,  &
            SH_rv,  &
            len(SH_rv, kind=C_INT))
        ! splicer end class.ExClass2.method.get_name
    end function exclass2_get_name

    ! const string & getName2()
    ! string_to_buffer_and_len
    ! function_index=21
    function exclass2_get_name2(obj) result(SH_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        implicit none
        class(exclass2) :: obj
        character(kind=C_CHAR, len=strlen_ptr(c_exclass2_get_name2_bufferify(  &
            obj%voidptr,  &
            SH_rv,  &
            len(SH_rv, kind=C_INT)))) :: SH_rv
        ! splicer begin class.ExClass2.method.get_name2
        call c_exclass2_get_name2_bufferify(  &
            obj%voidptr,  &
            SH_rv,  &
            len(SH_rv, kind=C_INT))
        ! splicer end class.ExClass2.method.get_name2
    end function exclass2_get_name2

    ! string & getName3() const
    ! string_to_buffer_and_len
    ! function_index=22
    function exclass2_get_name3(obj) result(SH_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        implicit none
        class(exclass2) :: obj
        character(kind=C_CHAR, len=strlen_ptr(c_exclass2_get_name3_bufferify(  &
            obj%voidptr,  &
            SH_rv,  &
            len(SH_rv, kind=C_INT)))) :: SH_rv
        ! splicer begin class.ExClass2.method.get_name3
        call c_exclass2_get_name3_bufferify(  &
            obj%voidptr,  &
            SH_rv,  &
            len(SH_rv, kind=C_INT))
        ! splicer end class.ExClass2.method.get_name3
    end function exclass2_get_name3

    ! string & getName4()
    ! string_to_buffer_and_len
    ! function_index=23
    function exclass2_get_name4(obj) result(SH_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        implicit none
        class(exclass2) :: obj
        character(kind=C_CHAR, len=strlen_ptr(c_exclass2_get_name4_bufferify(  &
            obj%voidptr,  &
            SH_rv,  &
            len(SH_rv, kind=C_INT)))) :: SH_rv
        ! splicer begin class.ExClass2.method.get_name4
        call c_exclass2_get_name4_bufferify(  &
            obj%voidptr,  &
            SH_rv,  &
            len(SH_rv, kind=C_INT))
        ! splicer end class.ExClass2.method.get_name4
    end function exclass2_get_name4

    ! const int GetNameLength()
    ! function_index=24
    !>
    !! \brief helper function for Fortran
    !!
    !<
    function exclass2_get_name_length(obj) result(SH_rv)
        use iso_c_binding, only : C_INT
        implicit none
        class(exclass2) :: obj
        integer(C_INT) :: SH_rv
        ! splicer begin class.ExClass2.method.get_name_length
        SH_rv = c_exclass2_get_name_length(obj%voidptr)
        ! splicer end class.ExClass2.method.get_name_length
    end function exclass2_get_name_length

    ! ExClass1 * get_class1(const ExClass1 * in+intent(in)+value)
    ! function_index=25
    function exclass2_get_class1(obj, in) result(SH_rv)
        use exclass1_mod, only : exclass1
        implicit none
        class(exclass2) :: obj
        type(exclass1), value, intent(IN) :: in
        type(exclass1) :: SH_rv
        ! splicer begin class.ExClass2.method.get_class1
        SH_rv%voidptr = c_exclass2_get_class1(  &
            obj%voidptr,  &
            in%yadda())
        ! splicer end class.ExClass2.method.get_class1
    end function exclass2_get_class1

    ! void * declare(TypeID type+intent(in)+value)
    ! fortran_generic - has_default_arg
    ! function_index=43
    subroutine exclass2_declare_0_int(obj, type)
        use iso_c_binding, only : C_INT
        implicit none
        class(exclass2) :: obj
        integer(C_INT), value, intent(IN) :: type
        ! splicer begin class.ExClass2.method.declare_0_int
        call c_exclass2_declare_0(  &
            obj%voidptr,  &
            type)
        ! splicer end class.ExClass2.method.declare_0_int
    end subroutine exclass2_declare_0_int

    ! void * declare(TypeID type+intent(in)+value)
    ! fortran_generic - has_default_arg
    ! function_index=44
    subroutine exclass2_declare_0_long(obj, type)
        use iso_c_binding, only : C_INT
        implicit none
        class(exclass2) :: obj
        integer(C_INT), value, intent(IN) :: type
        ! splicer begin class.ExClass2.method.declare_0_long
        call c_exclass2_declare_0(  &
            obj%voidptr,  &
            type)
        ! splicer end class.ExClass2.method.declare_0_long
    end subroutine exclass2_declare_0_long

    ! void * declare(TypeID type+intent(in)+value, int len+default(1)+intent(in)+value)
    ! fortran_generic
    ! function_index=45
    subroutine exclass2_declare_1_int(obj, type, len)
        use iso_c_binding, only : C_LONG, C_INT
        implicit none
        class(exclass2) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: len
        ! splicer begin class.ExClass2.method.declare_1_int
        call c_exclass2_declare_1(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.ExClass2.method.declare_1_int
    end subroutine exclass2_declare_1_int

    ! void * declare(TypeID type+intent(in)+value, long len+default(1)+intent(in)+value)
    ! fortran_generic
    ! function_index=46
    subroutine exclass2_declare_1_long(obj, type, len)
        use iso_c_binding, only : C_LONG, C_INT
        implicit none
        class(exclass2) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: len
        ! splicer begin class.ExClass2.method.declare_1_long
        call c_exclass2_declare_1(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.ExClass2.method.declare_1_long
    end subroutine exclass2_declare_1_long

    ! void destroyall()
    ! function_index=27
    subroutine exclass2_destroyall(obj)
        implicit none
        class(exclass2) :: obj
        ! splicer begin class.ExClass2.method.destroyall
        call c_exclass2_destroyall(obj%voidptr)
        ! splicer end class.ExClass2.method.destroyall
    end subroutine exclass2_destroyall

    ! TypeID getTypeID() const
    ! function_index=28
    function exclass2_get_type_id(obj) result(SH_rv)
        use iso_c_binding, only : C_INT
        implicit none
        class(exclass2) :: obj
        integer(C_INT) :: SH_rv
        ! splicer begin class.ExClass2.method.get_type_id
        SH_rv = c_exclass2_get_type_id(obj%voidptr)
        ! splicer end class.ExClass2.method.get_type_id
    end function exclass2_get_type_id

    ! void setValue(int value+intent(in)+value)
    ! cpp_template
    ! function_index=32
    subroutine exclass2_set_value_int(obj, value)
        use iso_c_binding, only : C_INT
        implicit none
        class(exclass2) :: obj
        integer(C_INT), value, intent(IN) :: value
        ! splicer begin class.ExClass2.method.set_value_int
        call c_exclass2_set_value_int(  &
            obj%voidptr,  &
            value)
        ! splicer end class.ExClass2.method.set_value_int
    end subroutine exclass2_set_value_int

    ! void setValue(long value+intent(in)+value)
    ! cpp_template
    ! function_index=33
    subroutine exclass2_set_value_long(obj, value)
        use iso_c_binding, only : C_LONG
        implicit none
        class(exclass2) :: obj
        integer(C_LONG), value, intent(IN) :: value
        ! splicer begin class.ExClass2.method.set_value_long
        call c_exclass2_set_value_long(  &
            obj%voidptr,  &
            value)
        ! splicer end class.ExClass2.method.set_value_long
    end subroutine exclass2_set_value_long

    ! void setValue(float value+intent(in)+value)
    ! cpp_template
    ! function_index=34
    subroutine exclass2_set_value_float(obj, value)
        use iso_c_binding, only : C_FLOAT
        implicit none
        class(exclass2) :: obj
        real(C_FLOAT), value, intent(IN) :: value
        ! splicer begin class.ExClass2.method.set_value_float
        call c_exclass2_set_value_float(  &
            obj%voidptr,  &
            value)
        ! splicer end class.ExClass2.method.set_value_float
    end subroutine exclass2_set_value_float

    ! void setValue(double value+intent(in)+value)
    ! cpp_template
    ! function_index=35
    subroutine exclass2_set_value_double(obj, value)
        use iso_c_binding, only : C_DOUBLE
        implicit none
        class(exclass2) :: obj
        real(C_DOUBLE), value, intent(IN) :: value
        ! splicer begin class.ExClass2.method.set_value_double
        call c_exclass2_set_value_double(  &
            obj%voidptr,  &
            value)
        ! splicer end class.ExClass2.method.set_value_double
    end subroutine exclass2_set_value_double

    ! int getValue()
    ! cpp_template
    ! function_index=36
    function exclass2_get_value_int(obj) result(SH_rv)
        use iso_c_binding, only : C_INT
        implicit none
        class(exclass2) :: obj
        integer(C_INT) :: SH_rv
        ! splicer begin class.ExClass2.method.get_value_int
        SH_rv = c_exclass2_get_value_int(obj%voidptr)
        ! splicer end class.ExClass2.method.get_value_int
    end function exclass2_get_value_int

    ! double getValue()
    ! cpp_template
    ! function_index=37
    function exclass2_get_value_double(obj) result(SH_rv)
        use iso_c_binding, only : C_DOUBLE
        implicit none
        class(exclass2) :: obj
        real(C_DOUBLE) :: SH_rv
        ! splicer begin class.ExClass2.method.get_value_double
        SH_rv = c_exclass2_get_value_double(obj%voidptr)
        ! splicer end class.ExClass2.method.get_value_double
    end function exclass2_get_value_double

    function exclass2_yadda(obj) result (voidptr)
        use iso_c_binding, only: C_PTR
        implicit none
        class(exclass2), intent(IN) :: obj
        type(C_PTR) :: voidptr
        voidptr = obj%voidptr
    end function exclass2_yadda

    function exclass2_associated(obj) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        class(exclass2), intent(IN) :: obj
        logical rv
        rv = c_associated(obj%voidptr)
    end function exclass2_associated

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

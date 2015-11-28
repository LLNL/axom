! wrapfExClass1.f
! This is generated code, do not edit
! blah blah
! yada yada
!
module exclass1_mod
    use fstr_mod
    use exclass2_mod, only : exclass2
    use iso_c_binding, only : C_INT, C_LONG
    ! splicer begin class.ExClass1.module_use
    ! splicer end class.ExClass1.module_use
    implicit none
    
    
    ! splicer begin class.ExClass1.module_top
    top of module splicer  1
    ! splicer end class.ExClass1.module_top
    
    type exclass1
        type(C_PTR) voidptr
        ! splicer begin class.ExClass1.component_part
          component part 1a
          component part 1b
        ! splicer end class.ExClass1.component_part
    contains
        procedure :: delete => exclass1_delete
        procedure :: increment_count => exclass1_increment_count
        procedure :: get_name => exclass1_get_name
        procedure :: get_name_length => exclass1_get_name_length
        procedure :: get_name_error_check => exclass1_get_name_error_check
        procedure :: get_name_arg => exclass1_get_name_arg
        procedure :: get_root => exclass1_get_root
        procedure :: get_value_from_int => exclass1_get_value_from_int
        procedure :: get_value_1 => exclass1_get_value_1
        procedure :: get_addr => exclass1_get_addr
        procedure :: has_addr => exclass1_has_addr
        procedure :: splicer_special => exclass1_splicer_special
        generic :: get_value => &
            ! splicer begin class.ExClass1.generic.get_value
            ! splicer end class.ExClass1.generic.get_value
            get_value_from_int,  &
            get_value_1
        ! splicer begin class.ExClass1.type_bound_procedure_part
          type bound procedure part 1
        ! splicer end class.ExClass1.type_bound_procedure_part
    end type exclass1
    
    
    interface operator (.eq.)
        module procedure exclass1_eq
    end interface
    
    interface operator (.ne.)
        module procedure exclass1_ne
    end interface
    
    interface
        
        function aa_exclass1_new(name) result(rv) &
                bind(C, name="AA_exclass1_new")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR) :: rv
        end function aa_exclass1_new
        
        function aa_exclass1_new_bufferify(name, Lname) result(rv) &
                bind(C, name="AA_exclass1_new_bufferify")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            type(C_PTR) :: rv
        end function aa_exclass1_new_bufferify
        
        subroutine aa_exclass1_delete(self) &
                bind(C, name="AA_exclass1_delete")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine aa_exclass1_delete
        
        function aa_exclass1_increment_count(self, incr) result(rv) &
                bind(C, name="AA_exclass1_increment_count")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: incr
            integer(C_INT) :: rv
        end function aa_exclass1_increment_count
        
        pure function aa_exclass1_get_name(self) result(rv) &
                bind(C, name="AA_exclass1_get_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) rv
        end function aa_exclass1_get_name
        
        pure function aa_exclass1_get_name_length(self) result(rv) &
                bind(C, name="AA_exclass1_get_name_length")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: rv
        end function aa_exclass1_get_name_length
        
        pure function aa_exclass1_get_name_error_check(self) result(rv) &
                bind(C, name="AA_exclass1_get_name_error_check")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) rv
        end function aa_exclass1_get_name_error_check
        
        pure function aa_exclass1_get_name_arg(self) result(rv) &
                bind(C, name="AA_exclass1_get_name_arg")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) rv
        end function aa_exclass1_get_name_arg
        
        function aa_exclass1_get_root(self) result(rv) &
                bind(C, name="AA_exclass1_get_root")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: rv
        end function aa_exclass1_get_root
        
        function aa_exclass1_get_value_from_int(self, value) result(rv) &
                bind(C, name="AA_exclass1_get_value_from_int")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: value
            integer(C_INT) :: rv
        end function aa_exclass1_get_value_from_int
        
        function aa_exclass1_get_value_1(self, value) result(rv) &
                bind(C, name="AA_exclass1_get_value_1")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_LONG), value, intent(IN) :: value
            integer(C_LONG) :: rv
        end function aa_exclass1_get_value_1
        
        function aa_exclass1_get_addr(self) result(rv) &
                bind(C, name="AA_exclass1_get_addr")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: rv
        end function aa_exclass1_get_addr
        
        function aa_exclass1_has_addr(self, in) result(rv) &
                bind(C, name="AA_exclass1_has_addr")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL), value, intent(IN) :: in
            logical(C_BOOL) :: rv
        end function aa_exclass1_has_addr
        
        subroutine aa_exclass1_splicer_special(self) &
                bind(C, name="AA_exclass1_splicer_special")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine aa_exclass1_splicer_special
        
        ! splicer begin class.ExClass1.additional_interfaces
        ! splicer end class.ExClass1.additional_interfaces
    end interface

contains
    
    ! ExClass1 *new ( const string *name ) +constructor
    ! string_to_buffer_and_len
    ! function_index=0
    function exclass1_new(name) result(rv)
        use iso_c_binding
        implicit none
        character(*), intent(IN) :: name
        type(exclass1) :: rv
        ! splicer begin class.ExClass1.method.new
        rv%voidptr = aa_exclass1_new_bufferify(  &
            name,  &
            len_trim(name))
        ! splicer end class.ExClass1.method.new
    end function exclass1_new
    
    ! void delete() +destructor
    ! function_index=1
    subroutine exclass1_delete(obj)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        ! splicer begin class.ExClass1.method.delete
        call aa_exclass1_delete(obj%voidptr)
        obj%voidptr = C_NULL_PTR
        ! splicer end class.ExClass1.method.delete
    end subroutine exclass1_delete
    
    ! int incrementCount(int incr)
    ! function_index=2
    function exclass1_increment_count(obj, incr) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        integer(C_INT), value, intent(IN) :: incr
        integer(C_INT) :: rv
        ! splicer begin class.ExClass1.method.increment_count
        rv = aa_exclass1_increment_count(  &
            obj%voidptr,  &
            incr)
        ! splicer end class.ExClass1.method.increment_count
    end function exclass1_increment_count
    
    ! const string& getName() const
    ! function_index=3
    function exclass1_get_name(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        character(kind=C_CHAR, len=aa_exclass1_get_name_length(obj%voidptr)) :: rv
        ! splicer begin class.ExClass1.method.get_name
        rv = fstr(aa_exclass1_get_name(obj%voidptr))
        ! splicer end class.ExClass1.method.get_name
    end function exclass1_get_name
    
    ! int GetNameLength() const
    ! function_index=4
    function exclass1_get_name_length(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        integer(C_INT) :: rv
        ! splicer begin class.ExClass1.method.get_name_length
        rv = aa_exclass1_get_name_length(obj%voidptr)
        ! splicer end class.ExClass1.method.get_name_length
    end function exclass1_get_name_length
    
    ! const string& getNameErrorCheck() const
    ! function_index=5
    function exclass1_get_name_error_check(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        character(kind=C_CHAR, len=strlen_ptr(aa_exclass1_get_name_error_check(obj%voidptr))) :: rv
        ! splicer begin class.ExClass1.method.get_name_error_check
        rv = fstr(aa_exclass1_get_name_error_check(obj%voidptr))
        ! splicer end class.ExClass1.method.get_name_error_check
    end function exclass1_get_name_error_check
    
    ! const string& getNameArg() const
    ! function_index=6
    subroutine exclass1_get_name_arg(obj, name)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        character(*), intent(OUT) :: name
        type(C_PTR) :: rv
        ! splicer begin class.ExClass1.method.get_name_arg
        rv = aa_exclass1_get_name_arg(obj%voidptr)
        call FccCopyPtr(name, len(name), rv)
        ! splicer end class.ExClass1.method.get_name_arg
    end subroutine exclass1_get_name_arg
    
    ! ExClass2 *getRoot()
    ! function_index=7
    function exclass1_get_root(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        type(exclass2) :: rv
        ! splicer begin class.ExClass1.method.get_root
        rv%voidptr = aa_exclass1_get_root(obj%voidptr)
        ! splicer end class.ExClass1.method.get_root
    end function exclass1_get_root
    
    ! int  getValue(int value)
    ! function_index=8
    function exclass1_get_value_from_int(obj, value) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        integer(C_INT), value, intent(IN) :: value
        integer(C_INT) :: rv
        ! splicer begin class.ExClass1.method.get_value_from_int
        rv = aa_exclass1_get_value_from_int(  &
            obj%voidptr,  &
            value)
        ! splicer end class.ExClass1.method.get_value_from_int
    end function exclass1_get_value_from_int
    
    ! long getValue(long value)
    ! function_index=9
    function exclass1_get_value_1(obj, value) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        integer(C_LONG), value, intent(IN) :: value
        integer(C_LONG) :: rv
        ! splicer begin class.ExClass1.method.get_value_1
        rv = aa_exclass1_get_value_1(  &
            obj%voidptr,  &
            value)
        ! splicer end class.ExClass1.method.get_value_1
    end function exclass1_get_value_1
    
    ! void *getAddr()
    ! function_index=10
    function exclass1_get_addr(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        type(C_PTR) :: rv
        ! splicer begin class.ExClass1.method.get_addr
        rv = aa_exclass1_get_addr(obj%voidptr)
        ! splicer end class.ExClass1.method.get_addr
    end function exclass1_get_addr
    
    ! bool hasAddr(bool in)
    ! function_index=11
    function exclass1_has_addr(obj, in) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        logical, value, intent(IN) :: in
        logical(C_BOOL) tmp_in
        logical :: rv
        tmp_in = in  ! coerce to C_BOOL
        ! splicer begin class.ExClass1.method.has_addr
        rv = aa_exclass1_has_addr(  &
            obj%voidptr,  &
            tmp_in)
        ! splicer end class.ExClass1.method.has_addr
    end function exclass1_has_addr
    
    ! void SplicerSpecial
    ! function_index=12
    subroutine exclass1_splicer_special(obj)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        ! splicer begin class.ExClass1.method.splicer_special
        blah blah blah
        ! splicer end class.ExClass1.method.splicer_special
    end subroutine exclass1_splicer_special
    
    ! splicer begin class.ExClass1.additional_functions
    ! splicer end class.ExClass1.additional_functions
    
    function exclass1_eq(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(exclass1), intent(IN) ::a,b
        logical :: rv
        if (c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function exclass1_eq
    
    function exclass1_ne(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(exclass1), intent(IN) ::a,b
        logical :: rv
        if (.not. c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function exclass1_ne

end module exclass1_mod

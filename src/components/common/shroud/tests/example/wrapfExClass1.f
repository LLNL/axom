! blah blah
! yada yada
!
module exclass1_mod
    use fstr_mod
    use exclass2_mod, only : exclass2
    use iso_c_binding
    
    type exclass1
        type(C_PTR) obj
    contains
        procedure :: increment_count => exclass1_increment_count
        procedure :: get_name => exclass1_get_name
        procedure :: get_name_length => exclass1_get_name_length
        procedure :: get_root => exclass1_get_root
        procedure :: get_addr => exclass1_get_addr
        procedure :: has_addr => exclass1_has_addr
        generic :: get_value => get_value_from_int, get_value_1
    end type exclass1
    
    interface
        
        function aa_exclass1_new(name) result(rv) &
                bind(C, name="AA_exclass1_new")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function aa_exclass1_new
        
        subroutine aa_exclass1_delete(self) &
                bind(C, name="AA_exclass1_delete")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
        end subroutine aa_exclass1_delete
        
        function aa_exclass1_increment_count(self, incr) result(rv) &
                bind(C, name="AA_exclass1_increment_count")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_INT), value :: incr
            integer(C_INT) :: rv
        end function aa_exclass1_increment_count
        
        pure function aa_exclass1_get_name(self) result(rv) &
                bind(C, name="AA_exclass1_get_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) rv
        end function aa_exclass1_get_name
        
        pure function aa_exclass1_get_name_length(self) result(rv) &
                bind(C, name="AA_exclass1_get_name_length")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_INT) :: rv
        end function aa_exclass1_get_name_length
        
        function aa_exclass1_get_root(self) result(rv) &
                bind(C, name="AA_exclass1_get_root")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function aa_exclass1_get_root
        
        function aa_exclass1_get_value_from_int(self, value) result(rv) &
                bind(C, name="AA_exclass1_get_value_from_int")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_INT), value :: value
            integer(C_INT) :: rv
        end function aa_exclass1_get_value_from_int
        
        function aa_exclass1_get_value_1(self, value) result(rv) &
                bind(C, name="AA_exclass1_get_value_1")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_LONG), value :: value
            integer(C_LONG) :: rv
        end function aa_exclass1_get_value_1
        
        function aa_exclass1_get_addr(self) result(rv) &
                bind(C, name="AA_exclass1_get_addr")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function aa_exclass1_get_addr
        
        function aa_exclass1_has_addr(self, in) result(rv) &
                bind(C, name="AA_exclass1_has_addr")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            logical(C_BOOL), value :: in
            logical(C_BOOL) :: rv
        end function aa_exclass1_has_addr
    end interface

contains
    
    function exclass1_new(name) result(rv)
        use iso_c_binding
        implicit none
        character(*) :: name
        type(exclass1) :: rv
        ! splicer begin
        rv%obj = aa_exclass1_new(trim(name) // C_NULL_CHAR)
        ! splicer end
    end function exclass1_new
    
    subroutine exclass1_delete(obj)
        use iso_c_binding
        implicit none
        type(exclass1) :: obj
        ! splicer begin
        call aa_exclass1_delete(obj%obj)
        obj%obj = C_NULL_PTR
        ! splicer end
    end subroutine exclass1_delete
    
    function exclass1_increment_count(obj, incr) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        integer(C_INT) :: incr
        integer(C_INT) :: rv
        ! splicer begin
        rv = aa_exclass1_increment_count(obj%obj, incr)
        ! splicer end
    end function exclass1_increment_count
    
    function exclass1_get_name(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        character(kind=C_CHAR, len=aa_exclass1_get_name_length(obj%obj)) :: rv
        type(C_PTR) :: rv_ptr
        ! splicer begin
        rv = fstr(aa_exclass1_get_name(obj%obj))
        ! splicer end
    end function exclass1_get_name
    
    function exclass1_get_name_length(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        integer(C_INT) :: rv
        ! splicer begin
        rv = aa_exclass1_get_name_length(obj%obj)
        ! splicer end
    end function exclass1_get_name_length
    
    function exclass1_get_root(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        type(exclass2) :: rv
        ! splicer begin
        rv%obj = aa_exclass1_get_root(obj%obj)
        ! splicer end
    end function exclass1_get_root
    
    function exclass1_get_value_from_int(obj, value) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        integer(C_INT) :: value
        integer(C_INT) :: rv
        ! splicer begin
        rv = aa_exclass1_get_value_from_int(obj%obj, value)
        ! splicer end
    end function exclass1_get_value_from_int
    
    function exclass1_get_value_1(obj, value) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        integer(C_LONG) :: value
        integer(C_LONG) :: rv
        ! splicer begin
        rv = aa_exclass1_get_value_1(obj%obj, value)
        ! splicer end
    end function exclass1_get_value_1
    
    function exclass1_get_addr(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        type(C_PTR) :: rv
        ! splicer begin
        rv = aa_exclass1_get_addr(obj%obj)
        ! splicer end
    end function exclass1_get_addr
    
    function exclass1_has_addr(obj, in) result(rv)
        use iso_c_binding
        implicit none
        class(exclass1) :: obj
        logical :: in
        logical :: rv
        ! splicer begin
        rv = bool2logical(aa_exclass1_has_addr(obj%obj, logical2bool(in)))
        ! splicer end
    end function exclass1_has_addr

end module exclass1_mod

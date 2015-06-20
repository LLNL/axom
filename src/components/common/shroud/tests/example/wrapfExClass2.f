! wrapfExClass2.f
! This is generated code, do not edit
! blah blah
! yada yada
!
module exclass2_mod
    use fstr_mod
    use exclass1_mod, only : exclass1
    use iso_c_binding, only : C_INT, C_LONG
    implicit none
    
    
    ! splicer begin class.ExClass2.module_top
    top of module splicer  2
    ! splicer end class.ExClass2.module_top
    
    type exclass2
        type(C_PTR) obj
        ! splicer begin class.ExClass2.component_part
        ! splicer end class.ExClass2.component_part
    contains
        procedure :: get_name => exclass2_get_name
        procedure :: get_name_length => exclass2_get_name_length
        procedure :: get_class1 => exclass2_get_class1
        procedure :: declare => exclass2_declare
        procedure :: destroyall => exclass2_destroyall
        procedure :: get_type_id => exclass2_get_type_id
        ! splicer begin class.ExClass2.type_bound_procedure_part
        ! splicer end class.ExClass2.type_bound_procedure_part
    end type exclass2
    
    interface
        
        function aa_exclass2_ex_class2(name) result(rv) &
                bind(C, name="AA_exclass2_ex_class2")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function aa_exclass2_ex_class2
        
        subroutine aa_exclass2_ex_class1(self) &
                bind(C, name="AA_exclass2_ex_class1")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine aa_exclass2_ex_class1
        
        pure function aa_exclass2_get_name(self) result(rv) &
                bind(C, name="AA_exclass2_get_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) rv
        end function aa_exclass2_get_name
        
        pure function aa_exclass2_get_name_length(self) result(rv) &
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
        
        function aa_exclass2_get_type_id(self) result(rv) &
                bind(C, name="AA_exclass2_get_type_id")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: rv
        end function aa_exclass2_get_type_id
    end interface

contains
    
    function exclass2_ex_class2(name) result(rv)
        use iso_c_binding
        implicit none
        character(*) :: name
        type(exclass2) :: rv
        ! splicer begin class.ExClass2.method.ex_class2
        rv%obj = aa_exclass2_ex_class2(trim(name) // C_NULL_CHAR)
        ! splicer end class.ExClass2.method.ex_class2
    end function exclass2_ex_class2
    
    subroutine exclass2_ex_class1(obj)
        use iso_c_binding
        implicit none
        type(exclass2) :: obj
        ! splicer begin class.ExClass2.method.ex_class1
        call aa_exclass2_ex_class1(obj%obj)
        obj%obj = C_NULL_PTR
        ! splicer end class.ExClass2.method.ex_class1
    end subroutine exclass2_ex_class1
    
    function exclass2_get_name(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        character(kind=C_CHAR, len=aa_exclass2_get_name_length(obj%obj)) :: rv
        ! splicer begin class.ExClass2.method.get_name
        rv = fstr(aa_exclass2_get_name(obj%obj))
        ! splicer end class.ExClass2.method.get_name
    end function exclass2_get_name
    
    function exclass2_get_name_length(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        integer(C_INT) :: rv
        ! splicer begin class.ExClass2.method.get_name_length
        rv = aa_exclass2_get_name_length(obj%obj)
        ! splicer end class.ExClass2.method.get_name_length
    end function exclass2_get_name_length
    
    function exclass2_get_class1(obj, in) result(rv)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        type(exclass1) :: in
        type(exclass1) :: rv
        ! splicer begin class.ExClass2.method.get_class1
        rv%obj = aa_exclass2_get_class1(obj%obj, in%obj)
        ! splicer end class.ExClass2.method.get_class1
    end function exclass2_get_class1
    
    subroutine exclass2_declare(obj, type, len)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.ExClass2.method.declare
        call aa_exclass2_declare(obj%obj, type, len)
        ! splicer end class.ExClass2.method.declare
    end subroutine exclass2_declare
    
    subroutine exclass2_destroyall(obj)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        ! splicer begin class.ExClass2.method.destroyall
        call aa_exclass2_destroyall(obj%obj)
        ! splicer end class.ExClass2.method.destroyall
    end subroutine exclass2_destroyall
    
    function exclass2_get_type_id(obj) result(rv)
        use iso_c_binding
        implicit none
        class(exclass2) :: obj
        integer(C_INT) :: rv
        ! splicer begin class.ExClass2.method.get_type_id
        rv = aa_exclass2_get_type_id(obj%obj)
        ! splicer end class.ExClass2.method.get_type_id
    end function exclass2_get_type_id
    ! splicer begin class.ExClass2.extra_methods
    ! splicer end class.ExClass2.extra_methods

end module exclass2_mod

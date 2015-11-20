! wrapfNames.f
! This is generated code, do not edit
module names_mod
    use fstr_mod
    ! splicer begin class.Names.module_use
    ! splicer end class.Names.module_use
    implicit none
    
    
    ! splicer begin class.Names.module_top
    ! splicer end class.Names.module_top
    
    type names
        type(C_PTR) voidptr
        ! splicer begin class.Names.component_part
        ! splicer end class.Names.component_part
    contains
        procedure :: method1 => names_method1
        procedure :: method2 => names_method2
        ! splicer begin class.Names.type_bound_procedure_part
        ! splicer end class.Names.type_bound_procedure_part
    end type names
    
    
    interface operator (.eq.)
        module procedure names_eq
    end interface
    
    interface operator (.ne.)
        module procedure names_ne
    end interface
    
    interface
        
        subroutine xxx_def_names_method1(self) &
                bind(C, name="XXX_DEF_names_method1")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine xxx_def_names_method1
        
        subroutine xxx_def_names_method2(self) &
                bind(C, name="XXX_DEF_names_method2")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine xxx_def_names_method2
        
        ! splicer begin class.Names.additional_interfaces
        ! splicer end class.Names.additional_interfaces
    end interface

contains
    
    subroutine names_method1(obj)
        use iso_c_binding
        implicit none
        class(names) :: obj
        ! splicer begin class.Names.method.method1
        call xxx_def_names_method1(obj%voidptr)
        ! splicer end class.Names.method.method1
    end subroutine names_method1
    
    subroutine names_method2(obj)
        use iso_c_binding
        implicit none
        class(names) :: obj
        ! splicer begin class.Names.method.method2
        call xxx_def_names_method2(obj%voidptr)
        ! splicer end class.Names.method.method2
    end subroutine names_method2
    
    ! splicer begin class.Names.additional_functions
    ! splicer end class.Names.additional_functions
    
    function names_eq(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(names), intent(IN) ::a,b
        logical :: rv
        if (c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function names_eq
    
    function names_ne(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(names), intent(IN) ::a,b
        logical :: rv
        if (.not. c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function names_ne

end module names_mod

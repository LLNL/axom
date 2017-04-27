! wrapfClass1.f
! This is generated code, do not edit
!>
!! \file wrapfClass1.f
!! \brief Shroud generated wrapper for Class1 class
!<
module class1_mod
    use iso_c_binding, only : C_INT, C_PTR
    implicit none



    type class1
        type(C_PTR), private :: voidptr
    contains
        procedure :: method1 => class1_method1
        procedure :: get_instance => class1_get_instance
        procedure :: set_instance => class1_set_instance
        procedure :: associated => class1_associated
    end type class1


    interface operator (.eq.)
        module procedure class1_eq
    end interface

    interface operator (.ne.)
        module procedure class1_ne
    end interface

    interface

        subroutine c_class1_method1(self, arg1) &
                bind(C, name="DEF_class1_method1")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: arg1
        end subroutine c_class1_method1

    end interface

contains

    subroutine class1_method1(obj, arg1)
        use iso_c_binding, only : C_INT
        class(class1) :: obj
        integer(C_INT), value, intent(IN) :: arg1
        call c_class1_method1(  &
            obj%voidptr,  &
            arg1)
    end subroutine class1_method1

    function class1_get_instance(obj) result (voidptr)
        use iso_c_binding, only: C_PTR
        implicit none
        class(class1), intent(IN) :: obj
        type(C_PTR) :: voidptr
        voidptr = obj%voidptr
    end function class1_get_instance

    subroutine class1_set_instance(obj, voidptr)
        use iso_c_binding, only: C_PTR
        implicit none
        class(class1), intent(INOUT) :: obj
        type(C_PTR), intent(IN) :: voidptr
        obj%voidptr = voidptr
    end subroutine class1_set_instance

    function class1_associated(obj) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        class(class1), intent(IN) :: obj
        logical rv
        rv = c_associated(obj%voidptr)
    end function class1_associated


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

end module class1_mod

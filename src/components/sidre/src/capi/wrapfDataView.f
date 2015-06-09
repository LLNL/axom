!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!
module dataview_mod
    use fstr_mod
    use databuffer_mod, only : databuffer
    use datagroup_mod, only : datagroup
    use iso_c_binding
    
    type dataview
        type(C_PTR) obj
    contains
        procedure :: declare => dataview_declare
        procedure :: allocate => dataview_allocate
        procedure :: has_buffer => dataview_has_buffer
        procedure :: get_name => dataview_get_name
        procedure :: get_buffer => dataview_get_buffer
        procedure :: get_data => dataview_get_data
        procedure :: get_owning_group => dataview_get_owning_group
    end type dataview
    
    interface
        
        function atk_dataview_declare(self, type, len) result(rv) bind(C, name="ATK_dataview_declare")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_INT), value :: type
            integer(C_LONG), value :: len
            type(C_PTR) :: rv
        end function atk_dataview_declare
        
        function atk_dataview_allocate(self, type, len) result(rv) bind(C, name="ATK_dataview_allocate")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_INT), value :: type
            integer(C_LONG), value :: len
            type(C_PTR) :: rv
        end function atk_dataview_allocate
        
        function atk_dataview_has_buffer(self) result(rv) bind(C, name="ATK_dataview_has_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            logical(C_BOOL) :: rv
        end function atk_dataview_has_buffer
        
        pure function atk_dataview_get_name(self) result(rv) bind(C, name="ATK_dataview_get_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) rv
        end function atk_dataview_get_name
        
        function atk_dataview_get_buffer(self) result(rv) bind(C, name="ATK_dataview_get_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function atk_dataview_get_buffer
        
        function atk_dataview_get_data(self) result(rv) bind(C, name="ATK_dataview_get_data")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function atk_dataview_get_data
        
        function atk_dataview_get_owning_group(self) result(rv) bind(C, name="ATK_dataview_get_owning_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function atk_dataview_get_owning_group
    end interface

contains
    
    function dataview_declare(obj, type, len) result(rv)
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        type(dataview) :: rv
        ! splicer begin
        rv%obj = atk_dataview_declare(obj%obj, type, len)
        ! splicer end
    end function dataview_declare
    
    function dataview_allocate(obj, type, len) result(rv)
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        type(dataview) :: rv
        ! splicer begin
        rv%obj = atk_dataview_allocate(obj%obj, type, len)
        ! splicer end
    end function dataview_allocate
    
    function dataview_has_buffer(obj) result(rv)
        implicit none
        class(dataview) :: obj
        logical :: rv
        ! splicer begin
        rv = bool2logical(atk_dataview_has_buffer(obj%obj))
        ! splicer end
    end function dataview_has_buffer
    
    function dataview_get_name(obj) result(rv)
        implicit none
        class(dataview) :: obj
        character(kind=C_CHAR, len=1) :: rv
        type(C_PTR) :: rv_ptr
        ! splicer begin
        rv = fstr(atk_dataview_get_name(obj%obj))
        ! splicer end
    end function dataview_get_name
    
    function dataview_get_buffer(obj) result(rv)
        implicit none
        class(dataview) :: obj
        type(databuffer) :: rv
        ! splicer begin
        rv%obj = atk_dataview_get_buffer(obj%obj)
        ! splicer end
    end function dataview_get_buffer
    
    function dataview_get_data(obj) result(rv)
        implicit none
        class(dataview) :: obj
        type(C_PTR) :: rv
        ! splicer begin
        rv = atk_dataview_get_data(obj%obj)
        ! splicer end
    end function dataview_get_data
    
    function dataview_get_owning_group(obj) result(rv)
        implicit none
        class(dataview) :: obj
        type(datagroup) :: rv
        ! splicer begin
        rv%obj = atk_dataview_get_owning_group(obj%obj)
        ! splicer end
    end function dataview_get_owning_group

end module dataview_mod

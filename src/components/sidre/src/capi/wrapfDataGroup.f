!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!
module datagroup_mod
    use fstr_mod
    use dataview_mod, only : dataview
    
    type datagroup
        type(C_PTR) obj
    contains
        procedure :: get_name => datagroup_get_name
        procedure :: create_view_and_buffer => datagroup_create_view_and_buffer
        procedure :: create_group => datagroup_create_group
    end type datagroup
    
    interface
        
        pure function atk_datagroup_get_name(self) result(rv) bind(C, name="ATK_datagroup_get_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) rv
        end function atk_datagroup_get_name
        
        function atk_datagroup_create_view_and_buffer(self, name) result(rv) bind(C, name="ATK_datagroup_create_view_and_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_create_view_and_buffer
        
        function atk_datagroup_create_group(self, name) result(rv) bind(C, name="ATK_datagroup_create_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_create_group
    end interface

contains
    
    function datagroup_get_name(obj) result(rv)
        implicit none
        class(datagroup) :: obj
        character(kind=C_CHAR, len=1) :: rv
        type(C_PTR) :: rv_ptr
        ! splicer begin
        rv = fstr(atk_datagroup_get_name(obj%obj))
        ! splicer end
    end function datagroup_get_name
    
    function datagroup_create_view_and_buffer(obj, name) result(rv)
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(dataview) :: rv
        ! splicer begin
        rv%obj = atk_datagroup_create_view_and_buffer(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end
    end function datagroup_create_view_and_buffer
    
    function datagroup_create_group(obj, name) result(rv)
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(datagroup) :: rv
        ! splicer begin
        rv%obj = atk_datagroup_create_group(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end
    end function datagroup_create_group

end module datagroup_mod

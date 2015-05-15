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
        procedure :: create_view_and_buffer => datagroup_create_view_and_buffer
        procedure :: create_group => datagroup_create_group
    end type datagroup
    
    interface
        
        function ds_datagroup_create_view_and_buffer(self, name) result(rv) bind(C, name="DS_datagroup_create_view_and_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function ds_datagroup_create_view_and_buffer
        
        function ds_datagroup_create_group(self, name) result(rv) bind(C, name="DS_datagroup_create_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function ds_datagroup_create_group
    end interface

contains
    
    function datagroup_create_view_and_buffer(obj, name) result(rv)
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(dataview) :: rv
        ! splicer begin
        rv%obj = ds_datagroup_create_view_and_buffer(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end
    end function datagroup_create_view_and_buffer
    
    function datagroup_create_group(obj, name) result(rv)
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(datagroup) :: rv
        ! splicer begin
        rv%obj = ds_datagroup_create_group(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end
    end function datagroup_create_group

end module datagroup_mod

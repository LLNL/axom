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
    use datastore_mod, only : datastore
    use dataview_mod, only : dataview
    use iso_c_binding
    
    type datagroup
        type(C_PTR) obj
    contains
        procedure :: get_name => datagroup_get_name
        procedure :: get_parent => datagroup_get_parent
        procedure :: get_data_store => datagroup_get_data_store
        procedure :: has_view => datagroup_has_view
        procedure :: create_view_and_buffer => datagroup_create_view_and_buffer
        procedure :: move_view => datagroup_move_view
        procedure :: copy_view => datagroup_copy_view
        procedure :: destroy_view_and_buffer => datagroup_destroy_view_and_buffer
        procedure :: get_view => datagroup_get_view
        procedure :: get_view_index => datagroup_get_view_index
        procedure :: get_view_name => datagroup_get_view_name
        procedure :: get_num_views => datagroup_get_num_views
        procedure :: has_group => datagroup_has_group
        procedure :: create_group => datagroup_create_group
        procedure :: destroy_group => datagroup_destroy_group
        procedure :: get_group_index => datagroup_get_group_index
        procedure :: get_group_name => datagroup_get_group_name
        procedure :: get_num_groups => datagroup_get_num_groups
        procedure :: print => datagroup_print
    end type datagroup
    
    interface
        
        pure function atk_datagroup_get_name(self) result(rv) bind(C, name="ATK_datagroup_get_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) rv
        end function atk_datagroup_get_name
        
        pure function atk_datagroup_get_parent(self) result(rv) bind(C, name="ATK_datagroup_get_parent")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function atk_datagroup_get_parent
        
        pure function atk_datagroup_get_data_store(self) result(rv) bind(C, name="ATK_datagroup_get_data_store")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: rv
        end function atk_datagroup_get_data_store
        
        function atk_datagroup_has_view(self, name) result(rv) bind(C, name="ATK_datagroup_has_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            character(kind=C_CHAR) :: name(*)
            logical(C_BOOL) :: rv
        end function atk_datagroup_has_view
        
        function atk_datagroup_create_view_and_buffer(self, name) result(rv) bind(C, name="ATK_datagroup_create_view_and_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_create_view_and_buffer
        
        function atk_datagroup_move_view(self, view) result(rv) bind(C, name="ATK_datagroup_move_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: view
            type(C_PTR) :: rv
        end function atk_datagroup_move_view
        
        function atk_datagroup_copy_view(self, view) result(rv) bind(C, name="ATK_datagroup_copy_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            type(C_PTR) :: view
            type(C_PTR) :: rv
        end function atk_datagroup_copy_view
        
        subroutine atk_datagroup_destroy_view_and_buffer(name) bind(C, name="ATK_datagroup_destroy_view_and_buffer")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR) :: name(*)
        end subroutine atk_datagroup_destroy_view_and_buffer
        
        function atk_datagroup_get_view(self, name) result(rv) bind(C, name="ATK_datagroup_get_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_get_view
        
        function atk_datagroup_get_view_index(self, name) result(rv) bind(C, name="ATK_datagroup_get_view_index")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            character(kind=C_CHAR) :: name(*)
            integer(C_INT) :: rv
        end function atk_datagroup_get_view_index
        
        pure function atk_datagroup_get_view_name(self, idx) result(rv) bind(C, name="ATK_datagroup_get_view_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_INT), value :: idx
            type(C_PTR) rv
        end function atk_datagroup_get_view_name
        
        function atk_datagroup_get_num_views(self) result(rv) bind(C, name="ATK_datagroup_get_num_views")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_SIZE_T) :: rv
        end function atk_datagroup_get_num_views
        
        function atk_datagroup_has_group(self, name) result(rv) bind(C, name="ATK_datagroup_has_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            character(kind=C_CHAR) :: name(*)
            logical(C_BOOL) :: rv
        end function atk_datagroup_has_group
        
        function atk_datagroup_create_group(self, name) result(rv) bind(C, name="ATK_datagroup_create_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_create_group
        
        subroutine atk_datagroup_destroy_group(name) bind(C, name="ATK_datagroup_destroy_group")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR) :: name(*)
        end subroutine atk_datagroup_destroy_group
        
        function atk_datagroup_get_group_index(self, name) result(rv) bind(C, name="ATK_datagroup_get_group_index")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            character(kind=C_CHAR) :: name(*)
            integer(C_INT) :: rv
        end function atk_datagroup_get_group_index
        
        pure function atk_datagroup_get_group_name(self, idx) result(rv) bind(C, name="ATK_datagroup_get_group_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_INT), value :: idx
            type(C_PTR) rv
        end function atk_datagroup_get_group_name
        
        function atk_datagroup_get_num_groups(self) result(rv) bind(C, name="ATK_datagroup_get_num_groups")
            use iso_c_binding
            implicit none
            type(C_PTR), value :: self
            integer(C_SIZE_T) :: rv
        end function atk_datagroup_get_num_groups
        
        subroutine atk_datagroup_print() bind(C, name="ATK_datagroup_print")
            use iso_c_binding
            implicit none
        end subroutine atk_datagroup_print
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
    
    function datagroup_get_parent(obj) result(rv)
        implicit none
        class(datagroup) :: obj
        type(datagroup) :: rv
        ! splicer begin
        rv%obj = atk_datagroup_get_parent(obj%obj)
        ! splicer end
    end function datagroup_get_parent
    
    function datagroup_get_data_store(obj) result(rv)
        implicit none
        class(datagroup) :: obj
        type(datastore) :: rv
        ! splicer begin
        rv%obj = atk_datagroup_get_data_store(obj%obj)
        ! splicer end
    end function datagroup_get_data_store
    
    function datagroup_has_view(obj, name) result(rv)
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        logical :: rv
        ! splicer begin
        rv = bool2logical(atk_datagroup_has_view(obj%obj, trim(name) // C_NULL_CHAR))
        ! splicer end
    end function datagroup_has_view
    
    function datagroup_create_view_and_buffer(obj, name) result(rv)
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(dataview) :: rv
        ! splicer begin
        rv%obj = atk_datagroup_create_view_and_buffer(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end
    end function datagroup_create_view_and_buffer
    
    function datagroup_move_view(obj, view) result(rv)
        implicit none
        class(datagroup) :: obj
        type(dataview) :: view
        type(dataview) :: rv
        ! splicer begin
        rv%obj = atk_datagroup_move_view(obj%obj, view)
        ! splicer end
    end function datagroup_move_view
    
    function datagroup_copy_view(obj, view) result(rv)
        implicit none
        class(datagroup) :: obj
        type(dataview) :: view
        type(dataview) :: rv
        ! splicer begin
        rv%obj = atk_datagroup_copy_view(obj%obj, view)
        ! splicer end
    end function datagroup_copy_view
    
    subroutine datagroup_destroy_view_and_buffer(name)
        implicit none
        character(*) :: name
        ! splicer begin
        call atk_datagroup_destroy_view_and_buffer(trim(name) // C_NULL_CHAR)
        ! splicer end
    end subroutine datagroup_destroy_view_and_buffer
    
    function datagroup_get_view(obj, name) result(rv)
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(dataview) :: rv
        ! splicer begin
        rv%obj = atk_datagroup_get_view(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end
    end function datagroup_get_view
    
    function datagroup_get_view_index(obj, name) result(rv)
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        integer(C_INT) :: rv
        ! splicer begin
        rv = atk_datagroup_get_view_index(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end
    end function datagroup_get_view_index
    
    function datagroup_get_view_name(obj, idx) result(rv)
        implicit none
        class(datagroup) :: obj
        integer(C_INT) :: idx
        character(kind=C_CHAR, len=1) :: rv
        type(C_PTR) :: rv_ptr
        ! splicer begin
        rv = fstr(atk_datagroup_get_view_name(obj%obj, idx))
        ! splicer end
    end function datagroup_get_view_name
    
    function datagroup_get_num_views(obj) result(rv)
        implicit none
        class(datagroup) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin
        rv = atk_datagroup_get_num_views(obj%obj)
        ! splicer end
    end function datagroup_get_num_views
    
    function datagroup_has_group(obj, name) result(rv)
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        logical :: rv
        ! splicer begin
        rv = bool2logical(atk_datagroup_has_group(obj%obj, trim(name) // C_NULL_CHAR))
        ! splicer end
    end function datagroup_has_group
    
    function datagroup_create_group(obj, name) result(rv)
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(datagroup) :: rv
        ! splicer begin
        rv%obj = atk_datagroup_create_group(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end
    end function datagroup_create_group
    
    subroutine datagroup_destroy_group(name)
        implicit none
        character(*) :: name
        ! splicer begin
        call atk_datagroup_destroy_group(trim(name) // C_NULL_CHAR)
        ! splicer end
    end subroutine datagroup_destroy_group
    
    function datagroup_get_group_index(obj, name) result(rv)
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        integer(C_INT) :: rv
        ! splicer begin
        rv = atk_datagroup_get_group_index(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end
    end function datagroup_get_group_index
    
    function datagroup_get_group_name(obj, idx) result(rv)
        implicit none
        class(datagroup) :: obj
        integer(C_INT) :: idx
        character(kind=C_CHAR, len=1) :: rv
        type(C_PTR) :: rv_ptr
        ! splicer begin
        rv = fstr(atk_datagroup_get_group_name(obj%obj, idx))
        ! splicer end
    end function datagroup_get_group_name
    
    function datagroup_get_num_groups(obj) result(rv)
        implicit none
        class(datagroup) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin
        rv = atk_datagroup_get_num_groups(obj%obj)
        ! splicer end
    end function datagroup_get_num_groups
    
    subroutine datagroup_print()
        implicit none
        ! splicer begin
        call atk_datagroup_print()
        ! splicer end
    end subroutine datagroup_print

end module datagroup_mod

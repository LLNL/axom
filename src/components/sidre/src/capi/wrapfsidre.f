!
! Copyright (c) 2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory.
!
! All rights reserved.
!
! This source code cannot be distributed without permission and
! further review from Lawrence Livermore National Laboratory.
!
module sidre_mod
    use fstr_mod
    use, intrinsic :: iso_c_binding, only : C_PTR
    implicit none
    
    ! splicer begin module_top
    !
    ! Type parameters
    ! Must be kept in sync with SidreTypes.h
    !
    integer, parameter :: ATK_INT8_T = 3
    integer, parameter :: ATK_INT16_T = 4
    integer, parameter :: ATK_INT32_T = 5
    integer, parameter :: ATK_INT64_T = 6
    integer, parameter :: ATK_UINT8_T = 7
    integer, parameter :: ATK_UINT16_T = 8
    integer, parameter :: ATK_UINT32_T = 9
    integer, parameter :: ATK_UINT64_T = 10
    integer, parameter :: ATK_FLOAT32_T = 11
    integer, parameter :: ATK_FLOAT64_T  = 12
    integer, parameter :: ATK_CHAR8_STR_T = 13
    
    integer, parameter :: ATK_C_INT_T = 14
    integer, parameter :: ATK_C_LONG_T = 15
    integer, parameter :: ATK_C_FLOAT_T = 16
    integer, parameter :: ATK_C_DOUBLE_T = 17
    ! splicer end module_top
    
    ! splicer begin class.datastore.module_top
    ! splicer end class.datastore.module_top
    
    type datastore
        type(C_PTR) obj
        ! splicer begin class.datastore.component_part
        ! splicer end class.datastore.component_part
    contains
        procedure :: get_root => datastore_get_root
        procedure :: get_buffer => datastore_get_buffer
        procedure :: create_buffer => datastore_create_buffer
        procedure :: destroy_buffer => datastore_destroy_buffer
        procedure :: get_num_buffers => datastore_get_num_buffers
        procedure :: print => datastore_print
        ! splicer begin class.datastore.type_bound_procedure_part
        ! splicer end class.datastore.type_bound_procedure_part
    end type datastore
    
    ! splicer begin class.datagroup.module_top
    ! splicer end class.datagroup.module_top
    
    type datagroup
        type(C_PTR) obj
        ! splicer begin class.datagroup.component_part
        ! splicer end class.datagroup.component_part
    contains
        procedure :: get_name => datagroup_get_name
        procedure :: get_parent => datagroup_get_parent
        procedure :: get_data_store => datagroup_get_data_store
        procedure :: get_num_views => datagroup_get_num_views
        procedure :: get_num_groups => datagroup_get_num_groups
        procedure :: has_view => datagroup_has_view
        procedure :: create_view_and_buffer_simple => datagroup_create_view_and_buffer_simple
        procedure :: create_view_and_buffer_from_type => datagroup_create_view_and_buffer_from_type
        procedure :: create_opaque_view => datagroup_create_opaque_view
        procedure :: create_view => datagroup_create_view
        procedure :: create_external_view => datagroup_create_external_view
        procedure :: move_view => datagroup_move_view
        procedure :: copy_view => datagroup_copy_view
        procedure :: destroy_view_and_buffer => datagroup_destroy_view_and_buffer
        procedure :: get_view => datagroup_get_view
        procedure :: get_view_index => datagroup_get_view_index
        procedure :: get_view_name => datagroup_get_view_name
        procedure :: has_group => datagroup_has_group
        procedure :: create_group => datagroup_create_group
        procedure :: move_group => datagroup_move_group
        procedure :: destroy_group => datagroup_destroy_group
        procedure :: get_group => datagroup_get_group
        procedure :: get_group_index => datagroup_get_group_index
        procedure :: get_group_name => datagroup_get_group_name
        procedure :: print => datagroup_print
        procedure :: save => datagroup_save
        procedure :: load => datagroup_load
        generic :: create_view_and_buffer => create_view_and_buffer_simple, create_view_and_buffer_from_type
        ! splicer begin class.datagroup.type_bound_procedure_part
        ! splicer end class.datagroup.type_bound_procedure_part
    end type datagroup
    
    ! splicer begin class.databuffer.module_top
    ! splicer end class.databuffer.module_top
    
    type databuffer
        type(C_PTR) obj
        ! splicer begin class.databuffer.component_part
        ! splicer end class.databuffer.component_part
    contains
        procedure :: get_index => databuffer_get_index
        procedure :: get_num_views => databuffer_get_num_views
        procedure :: declare => databuffer_declare
        procedure :: declare_external => databuffer_declare_external
        procedure :: allocate_existing => databuffer_allocate_existing
        procedure :: allocate_from_type => databuffer_allocate_from_type
        procedure :: reallocate => databuffer_reallocate
        procedure :: is_external => databuffer_is_external
        procedure :: get_data => databuffer_get_data
        procedure :: get_total_bytes => databuffer_get_total_bytes
        generic :: allocate => allocate_existing, allocate_from_type
        ! splicer begin class.databuffer.type_bound_procedure_part
        ! splicer end class.databuffer.type_bound_procedure_part
    end type databuffer
    
    ! splicer begin class.dataview.module_top
    ! splicer end class.dataview.module_top
    
    type dataview
        type(C_PTR) obj
        ! splicer begin class.dataview.component_part
        ! splicer end class.dataview.component_part
    contains
        procedure :: declare => dataview_declare
        procedure :: allocate => dataview_allocate
        procedure :: reallocate => dataview_reallocate
        procedure :: has_buffer => dataview_has_buffer
        procedure :: is_opaque => dataview_is_opaque
        procedure :: get_name => dataview_get_name
        procedure :: get_opaque => dataview_get_opaque
        procedure :: get_buffer => dataview_get_buffer
        procedure :: get_data_pointer => dataview_get_data_pointer
        procedure :: get_owning_group => dataview_get_owning_group
        procedure :: get_type_id => dataview_get_type_id
        procedure :: get_total_bytes => dataview_get_total_bytes
        procedure :: get_number_of_elements => dataview_get_number_of_elements
        ! splicer begin class.dataview.type_bound_procedure_part
        ! splicer end class.dataview.type_bound_procedure_part
    end type dataview
    
    interface
        
        function atk_datastore_new() result(rv) &
                bind(C, name="ATK_datastore_new")
            use iso_c_binding
            implicit none
            type(C_PTR) :: rv
        end function atk_datastore_new
        
        subroutine atk_datastore_delete(self) &
                bind(C, name="ATK_datastore_delete")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine atk_datastore_delete
        
        function atk_datastore_get_root(self) result(rv) &
                bind(C, name="ATK_datastore_get_root")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: rv
        end function atk_datastore_get_root
        
        function atk_datastore_get_buffer(self, idx) result(rv) &
                bind(C, name="ATK_datastore_get_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            type(C_PTR) :: rv
        end function atk_datastore_get_buffer
        
        function atk_datastore_create_buffer(self) result(rv) &
                bind(C, name="ATK_datastore_create_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: rv
        end function atk_datastore_create_buffer
        
        subroutine atk_datastore_destroy_buffer(self, id) &
                bind(C, name="ATK_datastore_destroy_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: id
        end subroutine atk_datastore_destroy_buffer
        
        function atk_datastore_get_num_buffers(self) result(rv) &
                bind(C, name="ATK_datastore_get_num_buffers")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: rv
        end function atk_datastore_get_num_buffers
        
        subroutine atk_datastore_print(self) &
                bind(C, name="ATK_datastore_print")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine atk_datastore_print
        
        pure function atk_datagroup_get_name(self) result(rv) &
                bind(C, name="ATK_datagroup_get_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) rv
        end function atk_datagroup_get_name
        
        pure function atk_datagroup_get_parent(self) result(rv) &
                bind(C, name="ATK_datagroup_get_parent")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: rv
        end function atk_datagroup_get_parent
        
        pure function atk_datagroup_get_data_store(self) result(rv) &
                bind(C, name="ATK_datagroup_get_data_store")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: rv
        end function atk_datagroup_get_data_store
        
        function atk_datagroup_get_num_views(self) result(rv) &
                bind(C, name="ATK_datagroup_get_num_views")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: rv
        end function atk_datagroup_get_num_views
        
        function atk_datagroup_get_num_groups(self) result(rv) &
                bind(C, name="ATK_datagroup_get_num_groups")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: rv
        end function atk_datagroup_get_num_groups
        
        function atk_datagroup_has_view(self, name) result(rv) &
                bind(C, name="ATK_datagroup_has_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
            logical(C_BOOL) :: rv
        end function atk_datagroup_has_view
        
        function atk_datagroup_create_view_and_buffer_simple(self, name) result(rv) &
                bind(C, name="ATK_datagroup_create_view_and_buffer_simple")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_create_view_and_buffer_simple
        
        function atk_datagroup_create_view_and_buffer_from_type(self, name, type, len) result(rv) &
                bind(C, name="ATK_datagroup_create_view_and_buffer_from_type")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
            type(C_PTR) :: rv
        end function atk_datagroup_create_view_and_buffer_from_type
        
        function atk_datagroup_create_opaque_view(self, name, opaque_ptr) result(rv) &
                bind(C, name="ATK_datagroup_create_opaque_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR), intent(IN) :: opaque_ptr
            type(C_PTR) :: rv
        end function atk_datagroup_create_opaque_view
        
        function atk_datagroup_create_view(self, name, buff) result(rv) &
                bind(C, name="ATK_datagroup_create_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR), intent(IN) :: buff
            type(C_PTR) :: rv
        end function atk_datagroup_create_view
        
        function atk_datagroup_create_external_view(self, name, external_data, type, len) result(rv) &
                bind(C, name="ATK_datagroup_create_external_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR), intent(IN) :: external_data
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
            type(C_PTR) :: rv
        end function atk_datagroup_create_external_view
        
        function atk_datagroup_move_view(self, view) result(rv) &
                bind(C, name="ATK_datagroup_move_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), intent(IN) :: view
            type(C_PTR) :: rv
        end function atk_datagroup_move_view
        
        function atk_datagroup_copy_view(self, view) result(rv) &
                bind(C, name="ATK_datagroup_copy_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), intent(IN) :: view
            type(C_PTR) :: rv
        end function atk_datagroup_copy_view
        
        subroutine atk_datagroup_destroy_view_and_buffer(self, name) &
                bind(C, name="ATK_datagroup_destroy_view_and_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
        end subroutine atk_datagroup_destroy_view_and_buffer
        
        function atk_datagroup_get_view(self, name) result(rv) &
                bind(C, name="ATK_datagroup_get_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_get_view
        
        function atk_datagroup_get_view_index(self, name) result(rv) &
                bind(C, name="ATK_datagroup_get_view_index")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
            integer(C_INT) :: rv
        end function atk_datagroup_get_view_index
        
        pure function atk_datagroup_get_view_name(self, idx) result(rv) &
                bind(C, name="ATK_datagroup_get_view_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            type(C_PTR) rv
        end function atk_datagroup_get_view_name
        
        function atk_datagroup_has_group(self, name) result(rv) &
                bind(C, name="ATK_datagroup_has_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
            logical(C_BOOL) :: rv
        end function atk_datagroup_has_group
        
        function atk_datagroup_create_group(self, name) result(rv) &
                bind(C, name="ATK_datagroup_create_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_create_group
        
        function atk_datagroup_move_group(self, grp) result(rv) &
                bind(C, name="ATK_datagroup_move_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), intent(IN) :: grp
            type(C_PTR) :: rv
        end function atk_datagroup_move_group
        
        subroutine atk_datagroup_destroy_group(self, name) &
                bind(C, name="ATK_datagroup_destroy_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
        end subroutine atk_datagroup_destroy_group
        
        function atk_datagroup_get_group(self, name) result(rv) &
                bind(C, name="ATK_datagroup_get_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_get_group
        
        function atk_datagroup_get_group_index(self, name) result(rv) &
                bind(C, name="ATK_datagroup_get_group_index")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: name(*)
            integer(C_INT) :: rv
        end function atk_datagroup_get_group_index
        
        pure function atk_datagroup_get_group_name(self, idx) result(rv) &
                bind(C, name="ATK_datagroup_get_group_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            type(C_PTR) rv
        end function atk_datagroup_get_group_name
        
        subroutine atk_datagroup_print(self) &
                bind(C, name="ATK_datagroup_print")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine atk_datagroup_print
        
        subroutine atk_datagroup_save(self, obase, protocol) &
                bind(C, name="ATK_datagroup_save")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: obase(*)
            character(kind=C_CHAR) :: protocol(*)
        end subroutine atk_datagroup_save
        
        subroutine atk_datagroup_load(self, obase, protocol) &
                bind(C, name="ATK_datagroup_load")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR) :: obase(*)
            character(kind=C_CHAR) :: protocol(*)
        end subroutine atk_datagroup_load
        
        function atk_databuffer_get_index(self) result(rv) &
                bind(C, name="ATK_databuffer_get_index")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: rv
        end function atk_databuffer_get_index
        
        function atk_databuffer_get_num_views(self) result(rv) &
                bind(C, name="ATK_databuffer_get_num_views")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: rv
        end function atk_databuffer_get_num_views
        
        subroutine atk_databuffer_declare(self, type, len) &
                bind(C, name="ATK_databuffer_declare")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
        end subroutine atk_databuffer_declare
        
        subroutine atk_databuffer_declare_external(self, external_data, type, len) &
                bind(C, name="ATK_databuffer_declare_external")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), intent(IN) :: external_data
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
        end subroutine atk_databuffer_declare_external
        
        subroutine atk_databuffer_allocate_existing(self) &
                bind(C, name="ATK_databuffer_allocate_existing")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine atk_databuffer_allocate_existing
        
        subroutine atk_databuffer_allocate_from_type(self, type, len) &
                bind(C, name="ATK_databuffer_allocate_from_type")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
        end subroutine atk_databuffer_allocate_from_type
        
        subroutine atk_databuffer_reallocate(self, type, len) &
                bind(C, name="ATK_databuffer_reallocate")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
        end subroutine atk_databuffer_reallocate
        
        function atk_databuffer_is_external(self) result(rv) &
                bind(C, name="ATK_databuffer_is_external")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: rv
        end function atk_databuffer_is_external
        
        function atk_databuffer_get_data(self) result(rv) &
                bind(C, name="ATK_databuffer_get_data")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: rv
        end function atk_databuffer_get_data
        
        function atk_databuffer_get_total_bytes(self) result(rv) &
                bind(C, name="ATK_databuffer_get_total_bytes")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: rv
        end function atk_databuffer_get_total_bytes
        
        subroutine atk_dataview_declare(self, type, len) &
                bind(C, name="ATK_dataview_declare")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
        end subroutine atk_dataview_declare
        
        subroutine atk_dataview_allocate(self, type, len) &
                bind(C, name="ATK_dataview_allocate")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
        end subroutine atk_dataview_allocate
        
        subroutine atk_dataview_reallocate(self, type, len) &
                bind(C, name="ATK_dataview_reallocate")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
        end subroutine atk_dataview_reallocate
        
        function atk_dataview_has_buffer(self) result(rv) &
                bind(C, name="ATK_dataview_has_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: rv
        end function atk_dataview_has_buffer
        
        function atk_dataview_is_opaque(self) result(rv) &
                bind(C, name="ATK_dataview_is_opaque")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: rv
        end function atk_dataview_is_opaque
        
        pure function atk_dataview_get_name(self) result(rv) &
                bind(C, name="ATK_dataview_get_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) rv
        end function atk_dataview_get_name
        
        function atk_dataview_get_opaque(self) result(rv) &
                bind(C, name="ATK_dataview_get_opaque")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: rv
        end function atk_dataview_get_opaque
        
        function atk_dataview_get_buffer(self) result(rv) &
                bind(C, name="ATK_dataview_get_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: rv
        end function atk_dataview_get_buffer
        
        function atk_dataview_get_data_pointer(self) result(rv) &
                bind(C, name="ATK_dataview_get_data_pointer")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: rv
        end function atk_dataview_get_data_pointer
        
        function atk_dataview_get_owning_group(self) result(rv) &
                bind(C, name="ATK_dataview_get_owning_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: rv
        end function atk_dataview_get_owning_group
        
        function atk_dataview_get_type_id(self) result(rv) &
                bind(C, name="ATK_dataview_get_type_id")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: rv
        end function atk_dataview_get_type_id
        
        function atk_dataview_get_total_bytes(self) result(rv) &
                bind(C, name="ATK_dataview_get_total_bytes")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: rv
        end function atk_dataview_get_total_bytes
        
        function atk_dataview_get_number_of_elements(self) result(rv) &
                bind(C, name="ATK_dataview_get_number_of_elements")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: rv
        end function atk_dataview_get_number_of_elements
    end interface

contains
    
    function datastore_new() result(rv)
        use iso_c_binding
        implicit none
        type(datastore) :: rv
        ! splicer begin class.datastore.method.new
        rv%obj = atk_datastore_new()
        ! splicer end class.datastore.method.new
    end function datastore_new
    
    subroutine datastore_delete(obj)
        use iso_c_binding
        implicit none
        type(datastore) :: obj
        ! splicer begin class.datastore.method.delete
        call atk_datastore_delete(obj%obj)
        obj%obj = C_NULL_PTR
        ! splicer end class.datastore.method.delete
    end subroutine datastore_delete
    
    function datastore_get_root(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datastore) :: obj
        type(datagroup) :: rv
        ! splicer begin class.datastore.method.get_root
        rv%obj = atk_datastore_get_root(obj%obj)
        ! splicer end class.datastore.method.get_root
    end function datastore_get_root
    
    function datastore_get_buffer(obj, idx) result(rv)
        use iso_c_binding
        implicit none
        class(datastore) :: obj
        integer(C_INT) :: idx
        type(databuffer) :: rv
        ! splicer begin class.datastore.method.get_buffer
        rv%obj = atk_datastore_get_buffer(obj%obj, idx)
        ! splicer end class.datastore.method.get_buffer
    end function datastore_get_buffer
    
    function datastore_create_buffer(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datastore) :: obj
        type(databuffer) :: rv
        ! splicer begin class.datastore.method.create_buffer
        rv%obj = atk_datastore_create_buffer(obj%obj)
        ! splicer end class.datastore.method.create_buffer
    end function datastore_create_buffer
    
    subroutine datastore_destroy_buffer(obj, id)
        use iso_c_binding
        implicit none
        class(datastore) :: obj
        integer(C_INT) :: id
        ! splicer begin class.datastore.method.destroy_buffer
        call atk_datastore_destroy_buffer(obj%obj, id)
        ! splicer end class.datastore.method.destroy_buffer
    end subroutine datastore_destroy_buffer
    
    function datastore_get_num_buffers(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datastore) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.datastore.method.get_num_buffers
        rv = atk_datastore_get_num_buffers(obj%obj)
        ! splicer end class.datastore.method.get_num_buffers
    end function datastore_get_num_buffers
    
    subroutine datastore_print(obj)
        use iso_c_binding
        implicit none
        class(datastore) :: obj
        ! splicer begin class.datastore.method.print
        call atk_datastore_print(obj%obj)
        ! splicer end class.datastore.method.print
    end subroutine datastore_print
    ! splicer begin class.datastore.extra_methods
    !  extra methods for datastore
    ! splicer end class.datastore.extra_methods
    
    function datagroup_get_name(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(kind=C_CHAR, len=1) :: rv
        ! splicer begin class.datagroup.method.get_name
        rv = fstr(atk_datagroup_get_name(obj%obj))
        ! splicer end class.datagroup.method.get_name
    end function datagroup_get_name
    
    function datagroup_get_parent(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(datagroup) :: rv
        ! splicer begin class.datagroup.method.get_parent
        rv%obj = atk_datagroup_get_parent(obj%obj)
        ! splicer end class.datagroup.method.get_parent
    end function datagroup_get_parent
    
    function datagroup_get_data_store(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(datastore) :: rv
        ! splicer begin class.datagroup.method.get_data_store
        rv%obj = atk_datagroup_get_data_store(obj%obj)
        ! splicer end class.datagroup.method.get_data_store
    end function datagroup_get_data_store
    
    function datagroup_get_num_views(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.datagroup.method.get_num_views
        rv = atk_datagroup_get_num_views(obj%obj)
        ! splicer end class.datagroup.method.get_num_views
    end function datagroup_get_num_views
    
    function datagroup_get_num_groups(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.datagroup.method.get_num_groups
        rv = atk_datagroup_get_num_groups(obj%obj)
        ! splicer end class.datagroup.method.get_num_groups
    end function datagroup_get_num_groups
    
    function datagroup_has_view(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        logical :: rv
        ! splicer begin class.datagroup.method.has_view
        rv = booltological(atk_datagroup_has_view(obj%obj, trim(name) // C_NULL_CHAR))
        ! splicer end class.datagroup.method.has_view
    end function datagroup_has_view
    
    function datagroup_create_view_and_buffer_simple(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(dataview) :: rv
        ! splicer begin class.datagroup.method.create_view_and_buffer_simple
        rv%obj = atk_datagroup_create_view_and_buffer_simple(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end class.datagroup.method.create_view_and_buffer_simple
    end function datagroup_create_view_and_buffer_simple
    
    function datagroup_create_view_and_buffer_from_type(obj, name, type, len) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        integer(C_INT) :: type
        integer(C_LONG) :: len
        type(dataview) :: rv
        ! splicer begin class.datagroup.method.create_view_and_buffer_from_type
        rv%obj = atk_datagroup_create_view_and_buffer_from_type(obj%obj, trim(name) // C_NULL_CHAR, type, len)
        ! splicer end class.datagroup.method.create_view_and_buffer_from_type
    end function datagroup_create_view_and_buffer_from_type
    
    function datagroup_create_opaque_view(obj, name, opaque_ptr) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(C_PTR) :: opaque_ptr
        type(dataview) :: rv
        ! splicer begin class.datagroup.method.create_opaque_view
        rv%obj = atk_datagroup_create_opaque_view(obj%obj, trim(name) // C_NULL_CHAR, opaque_ptr)
        ! splicer end class.datagroup.method.create_opaque_view
    end function datagroup_create_opaque_view
    
    function datagroup_create_view(obj, name, buff) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(databuffer) :: buff
        type(dataview) :: rv
        ! splicer begin class.datagroup.method.create_view
        rv%obj = atk_datagroup_create_view(obj%obj, trim(name) // C_NULL_CHAR, buff%obj)
        ! splicer end class.datagroup.method.create_view
    end function datagroup_create_view
    
    function datagroup_create_external_view(obj, name, external_data, type, len) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(C_PTR) :: external_data
        integer(C_INT) :: type
        integer(C_LONG) :: len
        type(dataview) :: rv
        ! splicer begin class.datagroup.method.create_external_view
        rv%obj = atk_datagroup_create_external_view(obj%obj, trim(name) // C_NULL_CHAR, external_data, type, len)
        ! splicer end class.datagroup.method.create_external_view
    end function datagroup_create_external_view
    
    function datagroup_move_view(obj, view) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(dataview) :: view
        type(dataview) :: rv
        ! splicer begin class.datagroup.method.move_view
        rv%obj = atk_datagroup_move_view(obj%obj, view%obj)
        ! splicer end class.datagroup.method.move_view
    end function datagroup_move_view
    
    function datagroup_copy_view(obj, view) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(dataview) :: view
        type(dataview) :: rv
        ! splicer begin class.datagroup.method.copy_view
        rv%obj = atk_datagroup_copy_view(obj%obj, view%obj)
        ! splicer end class.datagroup.method.copy_view
    end function datagroup_copy_view
    
    subroutine datagroup_destroy_view_and_buffer(obj, name)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        ! splicer begin class.datagroup.method.destroy_view_and_buffer
        call atk_datagroup_destroy_view_and_buffer(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end class.datagroup.method.destroy_view_and_buffer
    end subroutine datagroup_destroy_view_and_buffer
    
    function datagroup_get_view(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(dataview) :: rv
        ! splicer begin class.datagroup.method.get_view
        rv%obj = atk_datagroup_get_view(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end class.datagroup.method.get_view
    end function datagroup_get_view
    
    function datagroup_get_view_index(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        integer(C_INT) :: rv
        ! splicer begin class.datagroup.method.get_view_index
        rv = atk_datagroup_get_view_index(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end class.datagroup.method.get_view_index
    end function datagroup_get_view_index
    
    function datagroup_get_view_name(obj, idx) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        integer(C_INT) :: idx
        character(kind=C_CHAR, len=1) :: rv
        ! splicer begin class.datagroup.method.get_view_name
        rv = fstr(atk_datagroup_get_view_name(obj%obj, idx))
        ! splicer end class.datagroup.method.get_view_name
    end function datagroup_get_view_name
    
    function datagroup_has_group(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        logical :: rv
        ! splicer begin class.datagroup.method.has_group
        rv = booltological(atk_datagroup_has_group(obj%obj, trim(name) // C_NULL_CHAR))
        ! splicer end class.datagroup.method.has_group
    end function datagroup_has_group
    
    function datagroup_create_group(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(datagroup) :: rv
        ! splicer begin class.datagroup.method.create_group
        rv%obj = atk_datagroup_create_group(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end class.datagroup.method.create_group
    end function datagroup_create_group
    
    function datagroup_move_group(obj, grp) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(datagroup) :: grp
        type(datagroup) :: rv
        ! splicer begin class.datagroup.method.move_group
        rv%obj = atk_datagroup_move_group(obj%obj, grp%obj)
        ! splicer end class.datagroup.method.move_group
    end function datagroup_move_group
    
    subroutine datagroup_destroy_group(obj, name)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        ! splicer begin class.datagroup.method.destroy_group
        call atk_datagroup_destroy_group(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end class.datagroup.method.destroy_group
    end subroutine datagroup_destroy_group
    
    function datagroup_get_group(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(datagroup) :: rv
        ! splicer begin class.datagroup.method.get_group
        rv%obj = atk_datagroup_get_group(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end class.datagroup.method.get_group
    end function datagroup_get_group
    
    function datagroup_get_group_index(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        integer(C_INT) :: rv
        ! splicer begin class.datagroup.method.get_group_index
        rv = atk_datagroup_get_group_index(obj%obj, trim(name) // C_NULL_CHAR)
        ! splicer end class.datagroup.method.get_group_index
    end function datagroup_get_group_index
    
    function datagroup_get_group_name(obj, idx) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        integer(C_INT) :: idx
        character(kind=C_CHAR, len=1) :: rv
        ! splicer begin class.datagroup.method.get_group_name
        rv = fstr(atk_datagroup_get_group_name(obj%obj, idx))
        ! splicer end class.datagroup.method.get_group_name
    end function datagroup_get_group_name
    
    subroutine datagroup_print(obj)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        ! splicer begin class.datagroup.method.print
        call atk_datagroup_print(obj%obj)
        ! splicer end class.datagroup.method.print
    end subroutine datagroup_print
    
    subroutine datagroup_save(obj, obase, protocol)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: obase
        character(*) :: protocol
        ! splicer begin class.datagroup.method.save
        call atk_datagroup_save(obj%obj, trim(obase) // C_NULL_CHAR, trim(protocol) // C_NULL_CHAR)
        ! splicer end class.datagroup.method.save
    end subroutine datagroup_save
    
    subroutine datagroup_load(obj, obase, protocol)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: obase
        character(*) :: protocol
        ! splicer begin class.datagroup.method.load
        call atk_datagroup_load(obj%obj, trim(obase) // C_NULL_CHAR, trim(protocol) // C_NULL_CHAR)
        ! splicer end class.datagroup.method.load
    end subroutine datagroup_load
    ! splicer begin class.datagroup.extra_methods
    ! splicer end class.datagroup.extra_methods
    
    function databuffer_get_index(obj) result(rv)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: rv
        ! splicer begin class.databuffer.method.get_index
        rv = atk_databuffer_get_index(obj%obj)
        ! splicer end class.databuffer.method.get_index
    end function databuffer_get_index
    
    function databuffer_get_num_views(obj) result(rv)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.databuffer.method.get_num_views
        rv = atk_databuffer_get_num_views(obj%obj)
        ! splicer end class.databuffer.method.get_num_views
    end function databuffer_get_num_views
    
    subroutine databuffer_declare(obj, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.databuffer.method.declare
        call atk_databuffer_declare(obj%obj, type, len)
        ! splicer end class.databuffer.method.declare
    end subroutine databuffer_declare
    
    subroutine databuffer_declare_external(obj, external_data, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        type(C_PTR) :: external_data
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.databuffer.method.declare_external
        call atk_databuffer_declare_external(obj%obj, external_data, type, len)
        ! splicer end class.databuffer.method.declare_external
    end subroutine databuffer_declare_external
    
    subroutine databuffer_allocate_existing(obj)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        ! splicer begin class.databuffer.method.allocate_existing
        call atk_databuffer_allocate_existing(obj%obj)
        ! splicer end class.databuffer.method.allocate_existing
    end subroutine databuffer_allocate_existing
    
    subroutine databuffer_allocate_from_type(obj, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.databuffer.method.allocate_from_type
        call atk_databuffer_allocate_from_type(obj%obj, type, len)
        ! splicer end class.databuffer.method.allocate_from_type
    end subroutine databuffer_allocate_from_type
    
    subroutine databuffer_reallocate(obj, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.databuffer.method.reallocate
        call atk_databuffer_reallocate(obj%obj, type, len)
        ! splicer end class.databuffer.method.reallocate
    end subroutine databuffer_reallocate
    
    function databuffer_is_external(obj) result(rv)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        logical :: rv
        ! splicer begin class.databuffer.method.is_external
        rv = booltological(atk_databuffer_is_external(obj%obj))
        ! splicer end class.databuffer.method.is_external
    end function databuffer_is_external
    
    function databuffer_get_data(obj) result(rv)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        type(C_PTR) :: rv
        ! splicer begin class.databuffer.method.get_data
        rv = atk_databuffer_get_data(obj%obj)
        ! splicer end class.databuffer.method.get_data
    end function databuffer_get_data
    
    function databuffer_get_total_bytes(obj) result(rv)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.databuffer.method.get_total_bytes
        rv = atk_databuffer_get_total_bytes(obj%obj)
        ! splicer end class.databuffer.method.get_total_bytes
    end function databuffer_get_total_bytes
    ! splicer begin class.databuffer.extra_methods
    ! splicer end class.databuffer.extra_methods
    
    subroutine dataview_declare(obj, type, len)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.dataview.method.declare
        call atk_dataview_declare(obj%obj, type, len)
        ! splicer end class.dataview.method.declare
    end subroutine dataview_declare
    
    subroutine dataview_allocate(obj, type, len)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.dataview.method.allocate
        call atk_dataview_allocate(obj%obj, type, len)
        ! splicer end class.dataview.method.allocate
    end subroutine dataview_allocate
    
    subroutine dataview_reallocate(obj, type, len)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.dataview.method.reallocate
        call atk_dataview_reallocate(obj%obj, type, len)
        ! splicer end class.dataview.method.reallocate
    end subroutine dataview_reallocate
    
    function dataview_has_buffer(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        logical :: rv
        ! splicer begin class.dataview.method.has_buffer
        rv = booltological(atk_dataview_has_buffer(obj%obj))
        ! splicer end class.dataview.method.has_buffer
    end function dataview_has_buffer
    
    function dataview_is_opaque(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        logical :: rv
        ! splicer begin class.dataview.method.is_opaque
        rv = booltological(atk_dataview_is_opaque(obj%obj))
        ! splicer end class.dataview.method.is_opaque
    end function dataview_is_opaque
    
    function dataview_get_name(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        character(kind=C_CHAR, len=1) :: rv
        ! splicer begin class.dataview.method.get_name
        rv = fstr(atk_dataview_get_name(obj%obj))
        ! splicer end class.dataview.method.get_name
    end function dataview_get_name
    
    function dataview_get_opaque(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        type(C_PTR) :: rv
        ! splicer begin class.dataview.method.get_opaque
        rv = atk_dataview_get_opaque(obj%obj)
        ! splicer end class.dataview.method.get_opaque
    end function dataview_get_opaque
    
    function dataview_get_buffer(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        type(databuffer) :: rv
        ! splicer begin class.dataview.method.get_buffer
        rv%obj = atk_dataview_get_buffer(obj%obj)
        ! splicer end class.dataview.method.get_buffer
    end function dataview_get_buffer
    
    function dataview_get_data_pointer(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        type(C_PTR) :: rv
        ! splicer begin class.dataview.method.get_data_pointer
        rv = atk_dataview_get_data_pointer(obj%obj)
        ! splicer end class.dataview.method.get_data_pointer
    end function dataview_get_data_pointer
    
    function dataview_get_owning_group(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        type(datagroup) :: rv
        ! splicer begin class.dataview.method.get_owning_group
        rv%obj = atk_dataview_get_owning_group(obj%obj)
        ! splicer end class.dataview.method.get_owning_group
    end function dataview_get_owning_group
    
    function dataview_get_type_id(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: rv
        ! splicer begin class.dataview.method.get_type_id
        rv = atk_dataview_get_type_id(obj%obj)
        ! splicer end class.dataview.method.get_type_id
    end function dataview_get_type_id
    
    function dataview_get_total_bytes(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.dataview.method.get_total_bytes
        rv = atk_dataview_get_total_bytes(obj%obj)
        ! splicer end class.dataview.method.get_total_bytes
    end function dataview_get_total_bytes
    
    function dataview_get_number_of_elements(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.dataview.method.get_number_of_elements
        rv = atk_dataview_get_number_of_elements(obj%obj)
        ! splicer end class.dataview.method.get_number_of_elements
    end function dataview_get_number_of_elements
    ! splicer begin class.dataview.extra_methods
    ! splicer end class.dataview.extra_methods

end module sidre_mod

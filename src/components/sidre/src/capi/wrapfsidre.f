! wrapfsidre.f
! This is generated code, do not edit
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
    
    integer, parameter :: invalid_index = -1
    
    ! splicer end module_top
    
    ! splicer begin class.DataStore.module_top
    ! splicer end class.DataStore.module_top
    
    type datastore
        type(C_PTR) voidptr
        ! splicer begin class.DataStore.component_part
        ! splicer end class.DataStore.component_part
    contains
        procedure :: get_root => datastore_get_root
        procedure :: get_buffer => datastore_get_buffer
        procedure :: create_buffer => datastore_create_buffer
        procedure :: destroy_buffer => datastore_destroy_buffer
        procedure :: get_num_buffers => datastore_get_num_buffers
        procedure :: print => datastore_print
        ! splicer begin class.DataStore.type_bound_procedure_part
        ! splicer end class.DataStore.type_bound_procedure_part
    end type datastore
    
    ! splicer begin class.DataGroup.module_top
    ! splicer end class.DataGroup.module_top
    
    type datagroup
        type(C_PTR) voidptr
        ! splicer begin class.DataGroup.component_part
        ! splicer end class.DataGroup.component_part
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
        procedure :: get_group_name_length => datagroup_get_group_name_length
        procedure :: print => datagroup_print
        procedure :: save => datagroup_save
        procedure :: load => datagroup_load
        generic :: create_view_and_buffer => create_view_and_buffer_simple, create_view_and_buffer_from_type
        ! splicer begin class.DataGroup.type_bound_procedure_part
        ! splicer end class.DataGroup.type_bound_procedure_part
    end type datagroup
    
    ! splicer begin class.DataBuffer.module_top
    ! splicer end class.DataBuffer.module_top
    
    type databuffer
        type(C_PTR) voidptr
        ! splicer begin class.DataBuffer.component_part
        ! splicer end class.DataBuffer.component_part
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
        ! splicer begin class.DataBuffer.type_bound_procedure_part
        ! splicer end class.DataBuffer.type_bound_procedure_part
    end type databuffer
    
    ! splicer begin class.DataView.module_top
    ! splicer end class.DataView.module_top
    
    type dataview
        type(C_PTR) voidptr
        ! splicer begin class.DataView.component_part
        ! splicer end class.DataView.component_part
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
        ! splicer begin class.DataView.type_bound_procedure_part
        ! splicer end class.DataView.type_bound_procedure_part
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
        
        pure function atk_datastore_get_num_buffers(self) result(rv) &
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
        
        ! splicer begin class.DataStore.additional_interfaces
        ! splicer end class.DataStore.additional_interfaces
        
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
        
        pure function atk_datagroup_get_num_views(self) result(rv) &
                bind(C, name="ATK_datagroup_get_num_views")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: rv
        end function atk_datagroup_get_num_views
        
        pure function atk_datagroup_get_num_groups(self) result(rv) &
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
            character(kind=C_CHAR), intent(IN) :: name(*)
            logical(C_BOOL) :: rv
        end function atk_datagroup_has_view
        
        function atk_datagroup_create_view_and_buffer_simple(self, name) result(rv) &
                bind(C, name="ATK_datagroup_create_view_and_buffer_simple")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_create_view_and_buffer_simple
        
        function atk_datagroup_create_view_and_buffer_from_type(self, name, type, len) result(rv) &
                bind(C, name="ATK_datagroup_create_view_and_buffer_from_type")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
            type(C_PTR) :: rv
        end function atk_datagroup_create_view_and_buffer_from_type
        
        function atk_datagroup_create_opaque_view(self, name, opaque_ptr) result(rv) &
                bind(C, name="ATK_datagroup_create_opaque_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR), value, intent(IN) :: opaque_ptr
            type(C_PTR) :: rv
        end function atk_datagroup_create_opaque_view
        
        function atk_datagroup_create_view(self, name, buff) result(rv) &
                bind(C, name="ATK_datagroup_create_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR), value, intent(IN) :: buff
            type(C_PTR) :: rv
        end function atk_datagroup_create_view
        
        function atk_datagroup_create_external_view(self, name, external_data, type, len) result(rv) &
                bind(C, name="ATK_datagroup_create_external_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR), value, intent(IN) :: external_data
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
            type(C_PTR) :: rv
        end function atk_datagroup_create_external_view
        
        function atk_datagroup_move_view(self, view) result(rv) &
                bind(C, name="ATK_datagroup_move_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), value, intent(IN) :: view
            type(C_PTR) :: rv
        end function atk_datagroup_move_view
        
        function atk_datagroup_copy_view(self, view) result(rv) &
                bind(C, name="ATK_datagroup_copy_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), value, intent(IN) :: view
            type(C_PTR) :: rv
        end function atk_datagroup_copy_view
        
        subroutine atk_datagroup_destroy_view_and_buffer(self, name) &
                bind(C, name="ATK_datagroup_destroy_view_and_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
        end subroutine atk_datagroup_destroy_view_and_buffer
        
        function atk_datagroup_get_view(self, name) result(rv) &
                bind(C, name="ATK_datagroup_get_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_get_view
        
        pure function atk_datagroup_get_view_index(self, name) result(rv) &
                bind(C, name="ATK_datagroup_get_view_index")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
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
            character(kind=C_CHAR), intent(IN) :: name(*)
            logical(C_BOOL) :: rv
        end function atk_datagroup_has_group
        
        function atk_datagroup_create_group(self, name) result(rv) &
                bind(C, name="ATK_datagroup_create_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_create_group
        
        function atk_datagroup_move_group(self, grp) result(rv) &
                bind(C, name="ATK_datagroup_move_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), value, intent(IN) :: grp
            type(C_PTR) :: rv
        end function atk_datagroup_move_group
        
        subroutine atk_datagroup_destroy_group(self, name) &
                bind(C, name="ATK_datagroup_destroy_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
        end subroutine atk_datagroup_destroy_group
        
        function atk_datagroup_get_group(self, name) result(rv) &
                bind(C, name="ATK_datagroup_get_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_get_group
        
        pure function atk_datagroup_get_group_index(self, name) result(rv) &
                bind(C, name="ATK_datagroup_get_group_index")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
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
        
        pure function atk_datagroup_get_group_name_length(self, idx) result(rv) &
                bind(C, name="ATK_datagroup_get_group_name_length")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            integer(C_INT) :: rv
        end function atk_datagroup_get_group_name_length
        
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
            character(kind=C_CHAR), intent(IN) :: obase(*)
            character(kind=C_CHAR), intent(IN) :: protocol(*)
        end subroutine atk_datagroup_save
        
        subroutine atk_datagroup_load(self, obase, protocol) &
                bind(C, name="ATK_datagroup_load")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: obase(*)
            character(kind=C_CHAR), intent(IN) :: protocol(*)
        end subroutine atk_datagroup_load
        
        ! splicer begin class.DataGroup.additional_interfaces
        ! splicer end class.DataGroup.additional_interfaces
        
        pure function atk_databuffer_get_index(self) result(rv) &
                bind(C, name="ATK_databuffer_get_index")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: rv
        end function atk_databuffer_get_index
        
        pure function atk_databuffer_get_num_views(self) result(rv) &
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
            type(C_PTR), value, intent(IN) :: external_data
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
        
        pure function atk_databuffer_is_external(self) result(rv) &
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
        
        pure function atk_databuffer_get_total_bytes(self) result(rv) &
                bind(C, name="ATK_databuffer_get_total_bytes")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: rv
        end function atk_databuffer_get_total_bytes
        
        ! splicer begin class.DataBuffer.additional_interfaces
        ! splicer end class.DataBuffer.additional_interfaces
        
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
        
        pure function atk_dataview_has_buffer(self) result(rv) &
                bind(C, name="ATK_dataview_has_buffer")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: rv
        end function atk_dataview_has_buffer
        
        pure function atk_dataview_is_opaque(self) result(rv) &
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
        
        pure function atk_dataview_get_opaque(self) result(rv) &
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
        
        pure function atk_dataview_get_data_pointer(self) result(rv) &
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
        
        pure function atk_dataview_get_type_id(self) result(rv) &
                bind(C, name="ATK_dataview_get_type_id")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: rv
        end function atk_dataview_get_type_id
        
        pure function atk_dataview_get_total_bytes(self) result(rv) &
                bind(C, name="ATK_dataview_get_total_bytes")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: rv
        end function atk_dataview_get_total_bytes
        
        pure function atk_dataview_get_number_of_elements(self) result(rv) &
                bind(C, name="ATK_dataview_get_number_of_elements")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: rv
        end function atk_dataview_get_number_of_elements
        
        ! splicer begin class.DataView.additional_interfaces
        ! splicer end class.DataView.additional_interfaces
        
        function atk_is_name_valid(name) result(rv) &
                bind(C, name="ATK_is_name_valid")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            logical(C_BOOL) :: rv
        end function atk_is_name_valid
    end interface

contains
    
    function datastore_new() result(rv)
        use iso_c_binding
        implicit none
        type(datastore) :: rv
        ! splicer begin class.DataStore.method.new
        rv%voidptr = atk_datastore_new()
        ! splicer end class.DataStore.method.new
    end function datastore_new
    
    subroutine datastore_delete(obj)
        use iso_c_binding
        implicit none
        type(datastore) :: obj
        ! splicer begin class.DataStore.method.delete
        call atk_datastore_delete(obj%voidptr)
        obj%voidptr = C_NULL_PTR
        ! splicer end class.DataStore.method.delete
    end subroutine datastore_delete
    
    function datastore_get_root(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datastore) :: obj
        type(datagroup) :: rv
        ! splicer begin class.DataStore.method.get_root
        rv%voidptr = atk_datastore_get_root(obj%voidptr)
        ! splicer end class.DataStore.method.get_root
    end function datastore_get_root
    
    function datastore_get_buffer(obj, idx) result(rv)
        use iso_c_binding
        implicit none
        class(datastore) :: obj
        integer(C_INT) :: idx
        type(databuffer) :: rv
        ! splicer begin class.DataStore.method.get_buffer
        rv%voidptr = atk_datastore_get_buffer(obj%voidptr, idx)
        ! splicer end class.DataStore.method.get_buffer
    end function datastore_get_buffer
    
    function datastore_create_buffer(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datastore) :: obj
        type(databuffer) :: rv
        ! splicer begin class.DataStore.method.create_buffer
        rv%voidptr = atk_datastore_create_buffer(obj%voidptr)
        ! splicer end class.DataStore.method.create_buffer
    end function datastore_create_buffer
    
    subroutine datastore_destroy_buffer(obj, id)
        use iso_c_binding
        implicit none
        class(datastore) :: obj
        integer(C_INT) :: id
        ! splicer begin class.DataStore.method.destroy_buffer
        call atk_datastore_destroy_buffer(obj%voidptr, id)
        ! splicer end class.DataStore.method.destroy_buffer
    end subroutine datastore_destroy_buffer
    
    function datastore_get_num_buffers(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datastore) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.DataStore.method.get_num_buffers
        rv = atk_datastore_get_num_buffers(obj%voidptr)
        ! splicer end class.DataStore.method.get_num_buffers
    end function datastore_get_num_buffers
    
    subroutine datastore_print(obj)
        use iso_c_binding
        implicit none
        class(datastore) :: obj
        ! splicer begin class.DataStore.method.print
        call atk_datastore_print(obj%voidptr)
        ! splicer end class.DataStore.method.print
    end subroutine datastore_print
    
    ! splicer begin class.DataStore.additional_functions
    ! splicer end class.DataStore.additional_functions
    
    function datagroup_get_name(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(kind=C_CHAR, len=1) :: rv
        ! splicer begin class.DataGroup.method.get_name
        rv = fstr(atk_datagroup_get_name(obj%voidptr))
        ! splicer end class.DataGroup.method.get_name
    end function datagroup_get_name
    
    function datagroup_get_parent(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(datagroup) :: rv
        ! splicer begin class.DataGroup.method.get_parent
        rv%voidptr = atk_datagroup_get_parent(obj%voidptr)
        ! splicer end class.DataGroup.method.get_parent
    end function datagroup_get_parent
    
    function datagroup_get_data_store(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(datastore) :: rv
        ! splicer begin class.DataGroup.method.get_data_store
        rv%voidptr = atk_datagroup_get_data_store(obj%voidptr)
        ! splicer end class.DataGroup.method.get_data_store
    end function datagroup_get_data_store
    
    function datagroup_get_num_views(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.DataGroup.method.get_num_views
        rv = atk_datagroup_get_num_views(obj%voidptr)
        ! splicer end class.DataGroup.method.get_num_views
    end function datagroup_get_num_views
    
    function datagroup_get_num_groups(obj) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.DataGroup.method.get_num_groups
        rv = atk_datagroup_get_num_groups(obj%voidptr)
        ! splicer end class.DataGroup.method.get_num_groups
    end function datagroup_get_num_groups
    
    function datagroup_has_view(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        logical :: rv
        ! splicer begin class.DataGroup.method.has_view
        rv = booltological(atk_datagroup_has_view(obj%voidptr, trim(name) // C_NULL_CHAR))
        ! splicer end class.DataGroup.method.has_view
    end function datagroup_has_view
    
    function datagroup_create_view_and_buffer_simple(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.create_view_and_buffer_simple
        rv%voidptr = atk_datagroup_create_view_and_buffer_simple(obj%voidptr, trim(name) // C_NULL_CHAR)
        ! splicer end class.DataGroup.method.create_view_and_buffer_simple
    end function datagroup_create_view_and_buffer_simple
    
    function datagroup_create_view_and_buffer_from_type(obj, name, type, len) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        integer(C_INT) :: type
        integer(C_LONG) :: len
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.create_view_and_buffer_from_type
        rv%voidptr = atk_datagroup_create_view_and_buffer_from_type(obj%voidptr, trim(name) // C_NULL_CHAR, type, len)
        ! splicer end class.DataGroup.method.create_view_and_buffer_from_type
    end function datagroup_create_view_and_buffer_from_type
    
    function datagroup_create_opaque_view(obj, name, opaque_ptr) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(C_PTR) :: opaque_ptr
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.create_opaque_view
        rv%voidptr = atk_datagroup_create_opaque_view(obj%voidptr, trim(name) // C_NULL_CHAR, opaque_ptr)
        ! splicer end class.DataGroup.method.create_opaque_view
    end function datagroup_create_opaque_view
    
    function datagroup_create_view(obj, name, buff) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(databuffer) :: buff
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.create_view
        rv%voidptr = atk_datagroup_create_view(obj%voidptr, trim(name) // C_NULL_CHAR, buff%voidptr)
        ! splicer end class.DataGroup.method.create_view
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
        ! splicer begin class.DataGroup.method.create_external_view
        rv%voidptr = atk_datagroup_create_external_view(obj%voidptr, trim(name) // C_NULL_CHAR, external_data, type, len)
        ! splicer end class.DataGroup.method.create_external_view
    end function datagroup_create_external_view
    
    function datagroup_move_view(obj, view) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(dataview) :: view
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.move_view
        rv%voidptr = atk_datagroup_move_view(obj%voidptr, view%voidptr)
        ! splicer end class.DataGroup.method.move_view
    end function datagroup_move_view
    
    function datagroup_copy_view(obj, view) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(dataview) :: view
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.copy_view
        rv%voidptr = atk_datagroup_copy_view(obj%voidptr, view%voidptr)
        ! splicer end class.DataGroup.method.copy_view
    end function datagroup_copy_view
    
    subroutine datagroup_destroy_view_and_buffer(obj, name)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        ! splicer begin class.DataGroup.method.destroy_view_and_buffer
        call atk_datagroup_destroy_view_and_buffer(obj%voidptr, trim(name) // C_NULL_CHAR)
        ! splicer end class.DataGroup.method.destroy_view_and_buffer
    end subroutine datagroup_destroy_view_and_buffer
    
    function datagroup_get_view(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.get_view
        rv%voidptr = atk_datagroup_get_view(obj%voidptr, trim(name) // C_NULL_CHAR)
        ! splicer end class.DataGroup.method.get_view
    end function datagroup_get_view
    
    function datagroup_get_view_index(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        integer(C_INT) :: rv
        ! splicer begin class.DataGroup.method.get_view_index
        rv = atk_datagroup_get_view_index(obj%voidptr, trim(name) // C_NULL_CHAR)
        ! splicer end class.DataGroup.method.get_view_index
    end function datagroup_get_view_index
    
    function datagroup_get_view_name(obj, idx) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        integer(C_INT) :: idx
        character(kind=C_CHAR, len=1) :: rv
        ! splicer begin class.DataGroup.method.get_view_name
        rv = fstr(atk_datagroup_get_view_name(obj%voidptr, idx))
        ! splicer end class.DataGroup.method.get_view_name
    end function datagroup_get_view_name
    
    function datagroup_has_group(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        logical :: rv
        ! splicer begin class.DataGroup.method.has_group
        rv = booltological(atk_datagroup_has_group(obj%voidptr, trim(name) // C_NULL_CHAR))
        ! splicer end class.DataGroup.method.has_group
    end function datagroup_has_group
    
    function datagroup_create_group(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(datagroup) :: rv
        ! splicer begin class.DataGroup.method.create_group
        rv%voidptr = atk_datagroup_create_group(obj%voidptr, trim(name) // C_NULL_CHAR)
        ! splicer end class.DataGroup.method.create_group
    end function datagroup_create_group
    
    function datagroup_move_group(obj, grp) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(datagroup) :: grp
        type(datagroup) :: rv
        ! splicer begin class.DataGroup.method.move_group
        rv%voidptr = atk_datagroup_move_group(obj%voidptr, grp%voidptr)
        ! splicer end class.DataGroup.method.move_group
    end function datagroup_move_group
    
    subroutine datagroup_destroy_group(obj, name)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        ! splicer begin class.DataGroup.method.destroy_group
        call atk_datagroup_destroy_group(obj%voidptr, trim(name) // C_NULL_CHAR)
        ! splicer end class.DataGroup.method.destroy_group
    end subroutine datagroup_destroy_group
    
    function datagroup_get_group(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(datagroup) :: rv
        ! splicer begin class.DataGroup.method.get_group
        rv%voidptr = atk_datagroup_get_group(obj%voidptr, trim(name) // C_NULL_CHAR)
        ! splicer end class.DataGroup.method.get_group
    end function datagroup_get_group
    
    function datagroup_get_group_index(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        integer(C_INT) :: rv
        ! splicer begin class.DataGroup.method.get_group_index
        rv = atk_datagroup_get_group_index(obj%voidptr, trim(name) // C_NULL_CHAR)
        ! splicer end class.DataGroup.method.get_group_index
    end function datagroup_get_group_index
    
    function datagroup_get_group_name(obj, idx) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        integer(C_INT) :: idx
        character(kind=C_CHAR, len=atk_datagroup_get_group_name_length(obj%voidptr, idx)) :: rv
        ! splicer begin class.DataGroup.method.get_group_name
        rv = fstr(atk_datagroup_get_group_name(obj%voidptr, idx))
        ! splicer end class.DataGroup.method.get_group_name
    end function datagroup_get_group_name
    
    function datagroup_get_group_name_length(obj, idx) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        integer(C_INT) :: idx
        integer(C_INT) :: rv
        ! splicer begin class.DataGroup.method.get_group_name_length
        rv = atk_datagroup_get_group_name_length(obj%voidptr, idx)
        ! splicer end class.DataGroup.method.get_group_name_length
    end function datagroup_get_group_name_length
    
    subroutine datagroup_print(obj)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        ! splicer begin class.DataGroup.method.print
        call atk_datagroup_print(obj%voidptr)
        ! splicer end class.DataGroup.method.print
    end subroutine datagroup_print
    
    subroutine datagroup_save(obj, obase, protocol)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: obase
        character(*) :: protocol
        ! splicer begin class.DataGroup.method.save
        call atk_datagroup_save(obj%voidptr, trim(obase) // C_NULL_CHAR, trim(protocol) // C_NULL_CHAR)
        ! splicer end class.DataGroup.method.save
    end subroutine datagroup_save
    
    subroutine datagroup_load(obj, obase, protocol)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: obase
        character(*) :: protocol
        ! splicer begin class.DataGroup.method.load
        call atk_datagroup_load(obj%voidptr, trim(obase) // C_NULL_CHAR, trim(protocol) // C_NULL_CHAR)
        ! splicer end class.DataGroup.method.load
    end subroutine datagroup_load
    
    ! splicer begin class.DataGroup.additional_functions
    ! splicer end class.DataGroup.additional_functions
    
    function databuffer_get_index(obj) result(rv)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: rv
        ! splicer begin class.DataBuffer.method.get_index
        rv = atk_databuffer_get_index(obj%voidptr)
        ! splicer end class.DataBuffer.method.get_index
    end function databuffer_get_index
    
    function databuffer_get_num_views(obj) result(rv)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.DataBuffer.method.get_num_views
        rv = atk_databuffer_get_num_views(obj%voidptr)
        ! splicer end class.DataBuffer.method.get_num_views
    end function databuffer_get_num_views
    
    subroutine databuffer_declare(obj, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataBuffer.method.declare
        call atk_databuffer_declare(obj%voidptr, type, len)
        ! splicer end class.DataBuffer.method.declare
    end subroutine databuffer_declare
    
    subroutine databuffer_declare_external(obj, external_data, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        type(C_PTR) :: external_data
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataBuffer.method.declare_external
        call atk_databuffer_declare_external(obj%voidptr, external_data, type, len)
        ! splicer end class.DataBuffer.method.declare_external
    end subroutine databuffer_declare_external
    
    subroutine databuffer_allocate_existing(obj)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        ! splicer begin class.DataBuffer.method.allocate_existing
        call atk_databuffer_allocate_existing(obj%voidptr)
        ! splicer end class.DataBuffer.method.allocate_existing
    end subroutine databuffer_allocate_existing
    
    subroutine databuffer_allocate_from_type(obj, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataBuffer.method.allocate_from_type
        call atk_databuffer_allocate_from_type(obj%voidptr, type, len)
        ! splicer end class.DataBuffer.method.allocate_from_type
    end subroutine databuffer_allocate_from_type
    
    subroutine databuffer_reallocate(obj, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataBuffer.method.reallocate
        call atk_databuffer_reallocate(obj%voidptr, type, len)
        ! splicer end class.DataBuffer.method.reallocate
    end subroutine databuffer_reallocate
    
    function databuffer_is_external(obj) result(rv)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        logical :: rv
        ! splicer begin class.DataBuffer.method.is_external
        rv = booltological(atk_databuffer_is_external(obj%voidptr))
        ! splicer end class.DataBuffer.method.is_external
    end function databuffer_is_external
    
    function databuffer_get_data(obj) result(rv)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        type(C_PTR) :: rv
        ! splicer begin class.DataBuffer.method.get_data
        rv = atk_databuffer_get_data(obj%voidptr)
        ! splicer end class.DataBuffer.method.get_data
    end function databuffer_get_data
    
    function databuffer_get_total_bytes(obj) result(rv)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.DataBuffer.method.get_total_bytes
        rv = atk_databuffer_get_total_bytes(obj%voidptr)
        ! splicer end class.DataBuffer.method.get_total_bytes
    end function databuffer_get_total_bytes
    
    ! splicer begin class.DataBuffer.additional_functions
    ! splicer end class.DataBuffer.additional_functions
    
    subroutine dataview_declare(obj, type, len)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataView.method.declare
        call atk_dataview_declare(obj%voidptr, type, len)
        ! splicer end class.DataView.method.declare
    end subroutine dataview_declare
    
    subroutine dataview_allocate(obj, type, len)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataView.method.allocate
        call atk_dataview_allocate(obj%voidptr, type, len)
        ! splicer end class.DataView.method.allocate
    end subroutine dataview_allocate
    
    subroutine dataview_reallocate(obj, type, len)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataView.method.reallocate
        call atk_dataview_reallocate(obj%voidptr, type, len)
        ! splicer end class.DataView.method.reallocate
    end subroutine dataview_reallocate
    
    function dataview_has_buffer(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        logical :: rv
        ! splicer begin class.DataView.method.has_buffer
        rv = booltological(atk_dataview_has_buffer(obj%voidptr))
        ! splicer end class.DataView.method.has_buffer
    end function dataview_has_buffer
    
    function dataview_is_opaque(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        logical :: rv
        ! splicer begin class.DataView.method.is_opaque
        rv = booltological(atk_dataview_is_opaque(obj%voidptr))
        ! splicer end class.DataView.method.is_opaque
    end function dataview_is_opaque
    
    function dataview_get_name(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        character(kind=C_CHAR, len=1) :: rv
        ! splicer begin class.DataView.method.get_name
        rv = fstr(atk_dataview_get_name(obj%voidptr))
        ! splicer end class.DataView.method.get_name
    end function dataview_get_name
    
    function dataview_get_opaque(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        type(C_PTR) :: rv
        ! splicer begin class.DataView.method.get_opaque
        rv = atk_dataview_get_opaque(obj%voidptr)
        ! splicer end class.DataView.method.get_opaque
    end function dataview_get_opaque
    
    function dataview_get_buffer(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        type(databuffer) :: rv
        ! splicer begin class.DataView.method.get_buffer
        rv%voidptr = atk_dataview_get_buffer(obj%voidptr)
        ! splicer end class.DataView.method.get_buffer
    end function dataview_get_buffer
    
    function dataview_get_data_pointer(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        type(C_PTR) :: rv
        ! splicer begin class.DataView.method.get_data_pointer
        rv = atk_dataview_get_data_pointer(obj%voidptr)
        ! splicer end class.DataView.method.get_data_pointer
    end function dataview_get_data_pointer
    
    function dataview_get_owning_group(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        type(datagroup) :: rv
        ! splicer begin class.DataView.method.get_owning_group
        rv%voidptr = atk_dataview_get_owning_group(obj%voidptr)
        ! splicer end class.DataView.method.get_owning_group
    end function dataview_get_owning_group
    
    function dataview_get_type_id(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: rv
        ! splicer begin class.DataView.method.get_type_id
        rv = atk_dataview_get_type_id(obj%voidptr)
        ! splicer end class.DataView.method.get_type_id
    end function dataview_get_type_id
    
    function dataview_get_total_bytes(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.DataView.method.get_total_bytes
        rv = atk_dataview_get_total_bytes(obj%voidptr)
        ! splicer end class.DataView.method.get_total_bytes
    end function dataview_get_total_bytes
    
    function dataview_get_number_of_elements(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_SIZE_T) :: rv
        ! splicer begin class.DataView.method.get_number_of_elements
        rv = atk_dataview_get_number_of_elements(obj%voidptr)
        ! splicer end class.DataView.method.get_number_of_elements
    end function dataview_get_number_of_elements
    
    ! splicer begin class.DataView.additional_functions
    ! splicer end class.DataView.additional_functions
    
    function is_name_valid(name) result(rv)
        use iso_c_binding
        implicit none
        character(*) :: name
        logical :: rv
        ! splicer begin is_name_valid
        rv = name .ne. " "
        ! splicer end is_name_valid
    end function is_name_valid

end module sidre_mod

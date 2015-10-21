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
    interface c_loc_allocatable
    
       subroutine atk_c_loc_allocatable_int_scalar(variable, addr)
         use iso_c_binding
         integer(C_INT), allocatable, intent(IN) :: variable
         type(C_PTR), intent(OUT) :: addr
       end subroutine atk_c_loc_allocatable_int_scalar
    
       subroutine atk_c_loc_allocatable_int_1d(variable, addr)
         use iso_c_binding
         integer(C_INT), allocatable, intent(IN) :: variable(:)
         type(C_PTR), intent(OUT) :: addr
       end subroutine atk_c_loc_allocatable_int_1d
    
       subroutine atk_c_loc_allocatable_long_scalar(variable, addr)
         use iso_c_binding
         integer(C_LONG), allocatable, intent(IN) :: variable
         type(C_PTR), intent(OUT) :: addr
       end subroutine atk_c_loc_allocatable_long_scalar
    
       subroutine atk_c_loc_allocatable_long_1d(variable, addr)
         use iso_c_binding
         integer(C_LONG), allocatable, intent(IN) :: variable(:)
         type(C_PTR), intent(OUT) :: addr
       end subroutine atk_c_loc_allocatable_long_1d
    
       subroutine atk_c_loc_allocatable_float_scalar(variable, addr)
         use iso_c_binding
         real(C_FLOAT), allocatable, intent(IN) :: variable
         type(C_PTR), intent(OUT) :: addr
       end subroutine atk_c_loc_allocatable_float_scalar
    
       subroutine atk_c_loc_allocatable_float_1d(variable, addr)
         use iso_c_binding
         real(C_FLOAT), allocatable, intent(IN) :: variable(:)
         type(C_PTR), intent(OUT) :: addr
       end subroutine atk_c_loc_allocatable_float_1d
    
       subroutine atk_c_loc_allocatable_double_scalar(variable, addr)
         use iso_c_binding
         real(C_DOUBLE), allocatable, intent(IN) :: variable
         type(C_PTR), intent(OUT) :: addr
       end subroutine atk_c_loc_allocatable_double_scalar
    
       subroutine atk_c_loc_allocatable_double_1d(variable, addr)
         use iso_c_binding
         real(C_DOUBLE), allocatable, intent(IN) :: variable(:)
         type(C_PTR), intent(OUT) :: addr
       end subroutine atk_c_loc_allocatable_double_1d
    end interface c_loc_allocatable
    ! splicer end class.DataStore.module_top
    
    type datastore
        type(C_PTR) voidptr
        ! splicer begin class.DataStore.component_part
        ! splicer end class.DataStore.component_part
    contains
        procedure :: delete => datastore_delete
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
        procedure :: create_view_and_buffer_int => datagroup_create_view_and_buffer_int
        procedure :: create_view_and_buffer_long => datagroup_create_view_and_buffer_long
        procedure :: create_opaque_view => datagroup_create_opaque_view
        procedure :: create_view => datagroup_create_view
        procedure :: create_external_view_int => datagroup_create_external_view_int
        procedure :: create_external_view_long => datagroup_create_external_view_long
        procedure :: move_view => datagroup_move_view
        procedure :: copy_view => datagroup_copy_view
        procedure :: destroy_view_and_buffer => datagroup_destroy_view_and_buffer
        procedure :: get_view_from_name => datagroup_get_view_from_name
        procedure :: get_view_from_index => datagroup_get_view_from_index
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
        generic :: create_external_view => &
            ! splicer begin class.DataGroup.generic.create_external_view
            ! splicer end class.DataGroup.generic.create_external_view
            create_external_view_int,  &
            create_external_view_long
        generic :: create_view_and_buffer => &
            ! splicer begin class.DataGroup.generic.create_view_and_buffer
            ! splicer end class.DataGroup.generic.create_view_and_buffer
            create_view_and_buffer_simple,  &
            create_view_and_buffer_int,  &
            create_view_and_buffer_long
        generic :: get_view => &
            ! splicer begin class.DataGroup.generic.get_view
            ! splicer end class.DataGroup.generic.get_view
            get_view_from_name,  &
            get_view_from_index
        ! splicer begin class.DataGroup.type_bound_procedure_part
        procedure :: create_allocatable_view_int_scalar => datagroup_create_allocatable_view_int_scalar
        procedure :: create_allocatable_view_int_1d => datagroup_create_allocatable_view_int_1d
        procedure :: create_allocatable_view_long_scalar => datagroup_create_allocatable_view_long_scalar
        procedure :: create_allocatable_view_long_1d => datagroup_create_allocatable_view_long_1d
        procedure :: create_allocatable_view_float_scalar => datagroup_create_allocatable_view_float_scalar
        procedure :: create_allocatable_view_float_1d => datagroup_create_allocatable_view_float_1d
        procedure :: create_allocatable_view_double_scalar => datagroup_create_allocatable_view_double_scalar
        procedure :: create_allocatable_view_double_1d => datagroup_create_allocatable_view_double_1d
        generic :: create_allocatable_view => &
            create_allocatable_view_int_scalar,  &
            create_allocatable_view_int_1d,  &
            create_allocatable_view_long_scalar,  &
            create_allocatable_view_long_1d,  &
            create_allocatable_view_float_scalar,  &
            create_allocatable_view_float_1d,  &
            create_allocatable_view_double_scalar,  &
            create_allocatable_view_double_1d
        procedure :: create_array_view_int_scalar => datagroup_create_array_view_int_scalar
        procedure :: create_array_view_int_1d => datagroup_create_array_view_int_1d
        procedure :: create_array_view_long_scalar => datagroup_create_array_view_long_scalar
        procedure :: create_array_view_long_1d => datagroup_create_array_view_long_1d
        procedure :: create_array_view_float_scalar => datagroup_create_array_view_float_scalar
        procedure :: create_array_view_float_1d => datagroup_create_array_view_float_1d
        procedure :: create_array_view_double_scalar => datagroup_create_array_view_double_scalar
        procedure :: create_array_view_double_1d => datagroup_create_array_view_double_1d
        generic :: create_array_view => &
            create_array_view_int_scalar,  &
            create_array_view_int_1d,  &
            create_array_view_long_scalar,  &
            create_array_view_long_1d,  &
            create_array_view_float_scalar,  &
            create_array_view_float_1d,  &
            create_array_view_double_scalar,  &
            create_array_view_double_1d
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
        procedure :: declare_int => databuffer_declare_int
        procedure :: declare_long => databuffer_declare_long
        procedure :: allocate_existing => databuffer_allocate_existing
        procedure :: allocate_int => databuffer_allocate_int
        procedure :: allocate_long => databuffer_allocate_long
        procedure :: reallocate_int => databuffer_reallocate_int
        procedure :: reallocate_long => databuffer_reallocate_long
        procedure :: set_external_data => databuffer_set_external_data
        procedure :: is_external => databuffer_is_external
        procedure :: get_data => databuffer_get_data
        procedure :: get_total_bytes => databuffer_get_total_bytes
        procedure :: print => databuffer_print
        generic :: allocate => &
            ! splicer begin class.DataBuffer.generic.allocate
            ! splicer end class.DataBuffer.generic.allocate
            allocate_existing,  &
            allocate_int,  &
            allocate_long
        generic :: declare => &
            ! splicer begin class.DataBuffer.generic.declare
            ! splicer end class.DataBuffer.generic.declare
            declare_int,  &
            declare_long
        generic :: reallocate => &
            ! splicer begin class.DataBuffer.generic.reallocate
            ! splicer end class.DataBuffer.generic.reallocate
            reallocate_int,  &
            reallocate_long
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
        procedure :: declare_int => dataview_declare_int
        procedure :: declare_long => dataview_declare_long
        procedure :: allocate_simple => dataview_allocate_simple
        procedure :: allocate_int => dataview_allocate_int
        procedure :: allocate_long => dataview_allocate_long
        procedure :: reallocate_int => dataview_reallocate_int
        procedure :: reallocate_long => dataview_reallocate_long
        procedure :: has_buffer => dataview_has_buffer
        procedure :: is_opaque => dataview_is_opaque
        procedure :: get_name => dataview_get_name
        procedure :: get_opaque => dataview_get_opaque
        procedure :: get_buffer => dataview_get_buffer
        procedure :: get_data_pointer => dataview_get_data_pointer
        procedure :: set_value_int => dataview_set_value_int
        procedure :: set_value_long => dataview_set_value_long
        procedure :: set_value_float => dataview_set_value_float
        procedure :: set_value_double => dataview_set_value_double
        procedure :: get_value_int => dataview_get_value_int
        procedure :: get_value_long => dataview_get_value_long
        procedure :: get_value_float => dataview_get_value_float
        procedure :: get_value_double => dataview_get_value_double
        procedure :: get_owning_group => dataview_get_owning_group
        procedure :: get_type_id => dataview_get_type_id
        procedure :: get_total_bytes => dataview_get_total_bytes
        procedure :: get_number_of_elements => dataview_get_number_of_elements
        procedure :: print => dataview_print
        generic :: allocate => &
            ! splicer begin class.DataView.generic.allocate
            ! splicer end class.DataView.generic.allocate
            allocate_simple,  &
            allocate_int,  &
            allocate_long
        generic :: declare => &
            ! splicer begin class.DataView.generic.declare
            ! splicer end class.DataView.generic.declare
            declare_int,  &
            declare_long
        generic :: reallocate => &
            ! splicer begin class.DataView.generic.reallocate
            ! splicer end class.DataView.generic.reallocate
            reallocate_int,  &
            reallocate_long
        generic :: set_value => &
            ! splicer begin class.DataView.generic.set_value
            ! splicer end class.DataView.generic.set_value
            set_value_int,  &
            set_value_long,  &
            set_value_float,  &
            set_value_double
        ! splicer begin class.DataView.type_bound_procedure_part
        procedure :: get_value_int_scalar_ptr => dataview_get_value_int_scalar_ptr
        procedure :: get_value_int_1d_ptr => dataview_get_value_int_1d_ptr
        procedure :: get_value_long_scalar_ptr => dataview_get_value_long_scalar_ptr
        procedure :: get_value_long_1d_ptr => dataview_get_value_long_1d_ptr
        procedure :: get_value_float_scalar_ptr => dataview_get_value_float_scalar_ptr
        procedure :: get_value_float_1d_ptr => dataview_get_value_float_1d_ptr
        procedure :: get_value_double_scalar_ptr => dataview_get_value_double_scalar_ptr
        procedure :: get_value_double_1d_ptr => dataview_get_value_double_1d_ptr
        generic :: get_value => &
            get_value_int_scalar_ptr,  &
            get_value_int_1d_ptr,  &
            get_value_long_scalar_ptr,  &
            get_value_long_1d_ptr,  &
            get_value_float_scalar_ptr,  &
            get_value_float_1d_ptr,  &
            get_value_double_scalar_ptr,  &
            get_value_double_1d_ptr
        ! splicer end class.DataView.type_bound_procedure_part
    end type dataview
    
    
    interface operator (.eq.)
        module procedure datastore_eq
        module procedure datagroup_eq
        module procedure databuffer_eq
        module procedure dataview_eq
    end interface
    
    interface operator (.ne.)
        module procedure datastore_ne
        module procedure datagroup_ne
        module procedure databuffer_ne
        module procedure dataview_ne
    end interface
    
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
        
        function atk_datagroup_has_view_bufferify(self, name, Lname) result(rv) &
                bind(C, name="ATK_datagroup_has_view_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            logical(C_BOOL) :: rv
        end function atk_datagroup_has_view_bufferify
        
        function atk_datagroup_create_view_and_buffer_simple(self, name) result(rv) &
                bind(C, name="ATK_datagroup_create_view_and_buffer_simple")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_create_view_and_buffer_simple
        
        function atk_datagroup_create_view_and_buffer_bufferify(self, name, Lname) result(rv) &
                bind(C, name="ATK_datagroup_create_view_and_buffer_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            type(C_PTR) :: rv
        end function atk_datagroup_create_view_and_buffer_bufferify
        
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
        
        function atk_datagroup_create_opaque_view_bufferify(self, name, Lname, opaque_ptr) result(rv) &
                bind(C, name="ATK_datagroup_create_opaque_view_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            type(C_PTR), value, intent(IN) :: opaque_ptr
            type(C_PTR) :: rv
        end function atk_datagroup_create_opaque_view_bufferify
        
        function atk_datagroup_create_view(self, name, buff) result(rv) &
                bind(C, name="ATK_datagroup_create_view")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR), value, intent(IN) :: buff
            type(C_PTR) :: rv
        end function atk_datagroup_create_view
        
        function atk_datagroup_create_view_bufferify(self, name, Lname, buff) result(rv) &
                bind(C, name="ATK_datagroup_create_view_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            type(C_PTR), value, intent(IN) :: buff
            type(C_PTR) :: rv
        end function atk_datagroup_create_view_bufferify
        
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
        
        subroutine atk_datagroup_destroy_view_and_buffer_bufferify(self, name, Lname) &
                bind(C, name="ATK_datagroup_destroy_view_and_buffer_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
        end subroutine atk_datagroup_destroy_view_and_buffer_bufferify
        
        function atk_datagroup_get_view_from_name(self, name) result(rv) &
                bind(C, name="ATK_datagroup_get_view_from_name")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_get_view_from_name
        
        function atk_datagroup_get_view_bufferify(self, name, Lname) result(rv) &
                bind(C, name="ATK_datagroup_get_view_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            type(C_PTR) :: rv
        end function atk_datagroup_get_view_bufferify
        
        function atk_datagroup_get_view_from_index(self, idx) result(rv) &
                bind(C, name="ATK_datagroup_get_view_from_index")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            type(C_PTR) :: rv
        end function atk_datagroup_get_view_from_index
        
        pure function atk_datagroup_get_view_index(self, name) result(rv) &
                bind(C, name="ATK_datagroup_get_view_index")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT) :: rv
        end function atk_datagroup_get_view_index
        
        pure function atk_datagroup_get_view_index_bufferify(self, name, Lname) result(rv) &
                bind(C, name="ATK_datagroup_get_view_index_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            integer(C_INT) :: rv
        end function atk_datagroup_get_view_index_bufferify
        
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
        
        function atk_datagroup_has_group_bufferify(self, name, Lname) result(rv) &
                bind(C, name="ATK_datagroup_has_group_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            logical(C_BOOL) :: rv
        end function atk_datagroup_has_group_bufferify
        
        function atk_datagroup_create_group(self, name) result(rv) &
                bind(C, name="ATK_datagroup_create_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_create_group
        
        function atk_datagroup_create_group_bufferify(self, name, Lname) result(rv) &
                bind(C, name="ATK_datagroup_create_group_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            type(C_PTR) :: rv
        end function atk_datagroup_create_group_bufferify
        
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
        
        subroutine atk_datagroup_destroy_group_bufferify(self, name, Lname) &
                bind(C, name="ATK_datagroup_destroy_group_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
        end subroutine atk_datagroup_destroy_group_bufferify
        
        function atk_datagroup_get_group(self, name) result(rv) &
                bind(C, name="ATK_datagroup_get_group")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            type(C_PTR) :: rv
        end function atk_datagroup_get_group
        
        function atk_datagroup_get_group_bufferify(self, name, Lname) result(rv) &
                bind(C, name="ATK_datagroup_get_group_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            type(C_PTR) :: rv
        end function atk_datagroup_get_group_bufferify
        
        pure function atk_datagroup_get_group_index(self, name) result(rv) &
                bind(C, name="ATK_datagroup_get_group_index")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT) :: rv
        end function atk_datagroup_get_group_index
        
        pure function atk_datagroup_get_group_index_bufferify(self, name, Lname) result(rv) &
                bind(C, name="ATK_datagroup_get_group_index_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            integer(C_INT) :: rv
        end function atk_datagroup_get_group_index_bufferify
        
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
            character(kind=C_CHAR), intent(IN) :: obase(*)
            character(kind=C_CHAR), intent(IN) :: protocol(*)
        end subroutine atk_datagroup_save
        
        subroutine atk_datagroup_save_bufferify(self, obase, Lobase, protocol, Lprotocol) &
                bind(C, name="ATK_datagroup_save_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: obase(*)
            integer(C_INT), value, intent(IN) :: Lobase
            character(kind=C_CHAR), intent(IN) :: protocol(*)
            integer(C_INT), value, intent(IN) :: Lprotocol
        end subroutine atk_datagroup_save_bufferify
        
        subroutine atk_datagroup_load(self, obase, protocol) &
                bind(C, name="ATK_datagroup_load")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: obase(*)
            character(kind=C_CHAR), intent(IN) :: protocol(*)
        end subroutine atk_datagroup_load
        
        subroutine atk_datagroup_load_bufferify(self, obase, Lobase, protocol, Lprotocol) &
                bind(C, name="ATK_datagroup_load_bufferify")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: obase(*)
            integer(C_INT), value, intent(IN) :: Lobase
            character(kind=C_CHAR), intent(IN) :: protocol(*)
            integer(C_INT), value, intent(IN) :: Lprotocol
        end subroutine atk_datagroup_load_bufferify
        
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
        
        subroutine atk_databuffer_set_external_data(self, external_data) &
                bind(C, name="ATK_databuffer_set_external_data")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), value, intent(IN) :: external_data
        end subroutine atk_databuffer_set_external_data
        
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
        
        subroutine atk_databuffer_print(self) &
                bind(C, name="ATK_databuffer_print")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine atk_databuffer_print
        
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
        
        subroutine atk_dataview_allocate_simple(self) &
                bind(C, name="ATK_dataview_allocate_simple")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine atk_dataview_allocate_simple
        
        subroutine atk_dataview_allocate_from_type(self, type, len) &
                bind(C, name="ATK_dataview_allocate_from_type")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: len
        end subroutine atk_dataview_allocate_from_type
        
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
        
        subroutine atk_dataview_set_value_int(self, value) &
                bind(C, name="ATK_dataview_set_value_int")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: value
        end subroutine atk_dataview_set_value_int
        
        subroutine atk_dataview_set_value_long(self, value) &
                bind(C, name="ATK_dataview_set_value_long")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_LONG), value, intent(IN) :: value
        end subroutine atk_dataview_set_value_long
        
        subroutine atk_dataview_set_value_float(self, value) &
                bind(C, name="ATK_dataview_set_value_float")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_FLOAT), value, intent(IN) :: value
        end subroutine atk_dataview_set_value_float
        
        subroutine atk_dataview_set_value_double(self, value) &
                bind(C, name="ATK_dataview_set_value_double")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_DOUBLE), value, intent(IN) :: value
        end subroutine atk_dataview_set_value_double
        
        function atk_dataview_get_value_int(self) result(rv) &
                bind(C, name="ATK_dataview_get_value_int")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: rv
        end function atk_dataview_get_value_int
        
        function atk_dataview_get_value_long(self) result(rv) &
                bind(C, name="ATK_dataview_get_value_long")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_LONG) :: rv
        end function atk_dataview_get_value_long
        
        function atk_dataview_get_value_float(self) result(rv) &
                bind(C, name="ATK_dataview_get_value_float")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_FLOAT) :: rv
        end function atk_dataview_get_value_float
        
        function atk_dataview_get_value_double(self) result(rv) &
                bind(C, name="ATK_dataview_get_value_double")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_DOUBLE) :: rv
        end function atk_dataview_get_value_double
        
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
        
        subroutine atk_dataview_print(self) &
                bind(C, name="ATK_dataview_print")
            use iso_c_binding
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine atk_dataview_print
        
        ! splicer begin class.DataView.additional_interfaces
        ! splicer end class.DataView.additional_interfaces
        
        function atk_is_name_valid(name) result(rv) &
                bind(C, name="ATK_is_name_valid")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            logical(C_BOOL) :: rv
        end function atk_is_name_valid
        
        ! splicer begin additional_interfaces
        ! splicer end additional_interfaces
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
        class(datastore) :: obj
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
        rv%voidptr = atk_datastore_get_buffer(  &
            obj%voidptr,  &
            idx)
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
        call atk_datastore_destroy_buffer(  &
            obj%voidptr,  &
            id)
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
    
    subroutine datagroup_get_name(obj, name)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*), intent(OUT) :: name
        type(C_PTR) :: rv
        ! splicer begin class.DataGroup.method.get_name
        rv = atk_datagroup_get_name(obj%voidptr)
        call FccCopyPtr(name, len(name), rv)
        ! splicer end class.DataGroup.method.get_name
    end subroutine datagroup_get_name
    
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
        rv = atk_datagroup_has_view_bufferify(  &
            obj%voidptr,  &
            name,  &
            len_trim(name))
        ! splicer end class.DataGroup.method.has_view
    end function datagroup_has_view
    
    function datagroup_create_view_and_buffer_simple(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.create_view_and_buffer_simple
        rv%voidptr = atk_datagroup_create_view_and_buffer_bufferify(  &
            obj%voidptr,  &
            name,  &
            len_trim(name))
        ! splicer end class.DataGroup.method.create_view_and_buffer_simple
    end function datagroup_create_view_and_buffer_simple
    
    function datagroup_create_view_and_buffer_int(obj, name, type, len) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        integer(C_INT) :: type
        integer(C_INT) :: len
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.create_view_and_buffer_int
        rv%voidptr = atk_datagroup_create_view_and_buffer_from_type(  &
            obj%voidptr,  &
            trim(name) // C_NULL_CHAR,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataGroup.method.create_view_and_buffer_int
    end function datagroup_create_view_and_buffer_int
    
    function datagroup_create_view_and_buffer_long(obj, name, type, len) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        integer(C_INT) :: type
        integer(C_LONG) :: len
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.create_view_and_buffer_long
        rv%voidptr = atk_datagroup_create_view_and_buffer_from_type(  &
            obj%voidptr,  &
            trim(name) // C_NULL_CHAR,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataGroup.method.create_view_and_buffer_long
    end function datagroup_create_view_and_buffer_long
    
    function datagroup_create_opaque_view(obj, name, opaque_ptr) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(C_PTR) :: opaque_ptr
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.create_opaque_view
        rv%voidptr = atk_datagroup_create_opaque_view_bufferify(  &
            obj%voidptr,  &
            name,  &
            len_trim(name),  &
            opaque_ptr)
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
        rv%voidptr = atk_datagroup_create_view_bufferify(  &
            obj%voidptr,  &
            name,  &
            len_trim(name),  &
            buff%voidptr)
        ! splicer end class.DataGroup.method.create_view
    end function datagroup_create_view
    
    function datagroup_create_external_view_int(obj, name, external_data, type, len) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(C_PTR) :: external_data
        integer(C_INT) :: type
        integer(C_INT) :: len
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.create_external_view_int
        rv%voidptr = atk_datagroup_create_external_view(  &
            obj%voidptr,  &
            trim(name) // C_NULL_CHAR,  &
            external_data,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataGroup.method.create_external_view_int
    end function datagroup_create_external_view_int
    
    function datagroup_create_external_view_long(obj, name, external_data, type, len) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(C_PTR) :: external_data
        integer(C_INT) :: type
        integer(C_LONG) :: len
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.create_external_view_long
        rv%voidptr = atk_datagroup_create_external_view(  &
            obj%voidptr,  &
            trim(name) // C_NULL_CHAR,  &
            external_data,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataGroup.method.create_external_view_long
    end function datagroup_create_external_view_long
    
    function datagroup_move_view(obj, view) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(dataview) :: view
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.move_view
        rv%voidptr = atk_datagroup_move_view(  &
            obj%voidptr,  &
            view%voidptr)
        ! splicer end class.DataGroup.method.move_view
    end function datagroup_move_view
    
    function datagroup_copy_view(obj, view) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(dataview) :: view
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.copy_view
        rv%voidptr = atk_datagroup_copy_view(  &
            obj%voidptr,  &
            view%voidptr)
        ! splicer end class.DataGroup.method.copy_view
    end function datagroup_copy_view
    
    subroutine datagroup_destroy_view_and_buffer(obj, name)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        ! splicer begin class.DataGroup.method.destroy_view_and_buffer
        call atk_datagroup_destroy_view_and_buffer_bufferify(  &
            obj%voidptr,  &
            name,  &
            len_trim(name))
        ! splicer end class.DataGroup.method.destroy_view_and_buffer
    end subroutine datagroup_destroy_view_and_buffer
    
    function datagroup_get_view_from_name(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.get_view_from_name
        rv%voidptr = atk_datagroup_get_view_bufferify(  &
            obj%voidptr,  &
            name,  &
            len_trim(name))
        ! splicer end class.DataGroup.method.get_view_from_name
    end function datagroup_get_view_from_name
    
    function datagroup_get_view_from_index(obj, idx) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        integer(C_INT) :: idx
        type(dataview) :: rv
        ! splicer begin class.DataGroup.method.get_view_from_index
        rv%voidptr = atk_datagroup_get_view_from_index(  &
            obj%voidptr,  &
            idx)
        ! splicer end class.DataGroup.method.get_view_from_index
    end function datagroup_get_view_from_index
    
    function datagroup_get_view_index(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        integer(C_INT) :: rv
        ! splicer begin class.DataGroup.method.get_view_index
        rv = atk_datagroup_get_view_index_bufferify(  &
            obj%voidptr,  &
            name,  &
            len_trim(name))
        ! splicer end class.DataGroup.method.get_view_index
    end function datagroup_get_view_index
    
    subroutine datagroup_get_view_name(obj, idx, name)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        integer(C_INT) :: idx
        character(*), intent(OUT) :: name
        type(C_PTR) :: rv
        ! splicer begin class.DataGroup.method.get_view_name
        rv = atk_datagroup_get_view_name(  &
            obj%voidptr,  &
            idx)
        call FccCopyPtr(name, len(name), rv)
        ! splicer end class.DataGroup.method.get_view_name
    end subroutine datagroup_get_view_name
    
    function datagroup_has_group(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        logical :: rv
        ! splicer begin class.DataGroup.method.has_group
        rv = atk_datagroup_has_group_bufferify(  &
            obj%voidptr,  &
            name,  &
            len_trim(name))
        ! splicer end class.DataGroup.method.has_group
    end function datagroup_has_group
    
    function datagroup_create_group(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(datagroup) :: rv
        ! splicer begin class.DataGroup.method.create_group
        rv%voidptr = atk_datagroup_create_group_bufferify(  &
            obj%voidptr,  &
            name,  &
            len_trim(name))
        ! splicer end class.DataGroup.method.create_group
    end function datagroup_create_group
    
    function datagroup_move_group(obj, grp) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        type(datagroup) :: grp
        type(datagroup) :: rv
        ! splicer begin class.DataGroup.method.move_group
        rv%voidptr = atk_datagroup_move_group(  &
            obj%voidptr,  &
            grp%voidptr)
        ! splicer end class.DataGroup.method.move_group
    end function datagroup_move_group
    
    subroutine datagroup_destroy_group(obj, name)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        ! splicer begin class.DataGroup.method.destroy_group
        call atk_datagroup_destroy_group_bufferify(  &
            obj%voidptr,  &
            name,  &
            len_trim(name))
        ! splicer end class.DataGroup.method.destroy_group
    end subroutine datagroup_destroy_group
    
    function datagroup_get_group(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        type(datagroup) :: rv
        ! splicer begin class.DataGroup.method.get_group
        rv%voidptr = atk_datagroup_get_group_bufferify(  &
            obj%voidptr,  &
            name,  &
            len_trim(name))
        ! splicer end class.DataGroup.method.get_group
    end function datagroup_get_group
    
    function datagroup_get_group_index(obj, name) result(rv)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: name
        integer(C_INT) :: rv
        ! splicer begin class.DataGroup.method.get_group_index
        rv = atk_datagroup_get_group_index_bufferify(  &
            obj%voidptr,  &
            name,  &
            len_trim(name))
        ! splicer end class.DataGroup.method.get_group_index
    end function datagroup_get_group_index
    
    subroutine datagroup_get_group_name(obj, idx, name)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        integer(C_INT) :: idx
        character(*), intent(OUT) :: name
        type(C_PTR) :: rv
        ! splicer begin class.DataGroup.method.get_group_name
        rv = atk_datagroup_get_group_name(  &
            obj%voidptr,  &
            idx)
        call FccCopyPtr(name, len(name), rv)
        ! splicer end class.DataGroup.method.get_group_name
    end subroutine datagroup_get_group_name
    
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
        call atk_datagroup_save_bufferify(  &
            obj%voidptr,  &
            obase,  &
            len_trim(obase),  &
            protocol,  &
            len_trim(protocol))
        ! splicer end class.DataGroup.method.save
    end subroutine datagroup_save
    
    subroutine datagroup_load(obj, obase, protocol)
        use iso_c_binding
        implicit none
        class(datagroup) :: obj
        character(*) :: obase
        character(*) :: protocol
        ! splicer begin class.DataGroup.method.load
        call atk_datagroup_load_bufferify(  &
            obj%voidptr,  &
            obase,  &
            len_trim(obase),  &
            protocol,  &
            len_trim(protocol))
        ! splicer end class.DataGroup.method.load
    end subroutine datagroup_load
    
    ! splicer begin class.DataGroup.additional_functions
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_allocatable_view_int_scalar(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_fortran_allocatable_view(group, name, lname, addr, itype, rank) &
              bind(C,name="ATK_create_fortran_allocatable_view") &
              result(rv)
           use iso_c_binding
           type(C_PTR), value, intent(IN)    :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN) :: lname
           type(C_PTR), value                :: addr
           integer(C_INT), value, intent(IN) :: itype
           integer(C_INT), value, intent(IN) :: rank
           type(C_PTR) rv
           end function ATK_create_fortran_allocatable_view
        end interface
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        integer(C_INT), allocatable, intent(IN) :: value
        integer(C_INT) :: lname
        type(dataview) :: rv
        type(C_PTR) :: addr
        integer(C_INT), parameter :: rank = 0
        integer(C_INT), parameter :: itype = ATK_C_INT_T
    
        lname = len_trim(name)
        call c_loc_allocatable(value, addr)
        rv%voidptr = ATK_create_fortran_allocatable_view(group%voidptr, name, lname, addr, itype, rank)
    end function datagroup_create_allocatable_view_int_scalar
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_allocatable_view_int_1d(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_fortran_allocatable_view(group, name, lname, addr, itype, rank) &
              bind(C,name="ATK_create_fortran_allocatable_view") &
              result(rv)
           use iso_c_binding
           type(C_PTR), value, intent(IN)    :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN) :: lname
           type(C_PTR), value                :: addr
           integer(C_INT), value, intent(IN) :: itype
           integer(C_INT), value, intent(IN) :: rank
           type(C_PTR) rv
           end function ATK_create_fortran_allocatable_view
        end interface
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        integer(C_INT), allocatable, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(dataview) :: rv
        type(C_PTR) :: addr
        integer(C_INT), parameter :: rank = 1
        integer(C_INT), parameter :: itype = ATK_C_INT_T
    
        lname = len_trim(name)
        call c_loc_allocatable(value, addr)
        rv%voidptr = ATK_create_fortran_allocatable_view(group%voidptr, name, lname, addr, itype, rank)
    end function datagroup_create_allocatable_view_int_1d
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_allocatable_view_long_scalar(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_fortran_allocatable_view(group, name, lname, addr, itype, rank) &
              bind(C,name="ATK_create_fortran_allocatable_view") &
              result(rv)
           use iso_c_binding
           type(C_PTR), value, intent(IN)    :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN) :: lname
           type(C_PTR), value                :: addr
           integer(C_INT), value, intent(IN) :: itype
           integer(C_INT), value, intent(IN) :: rank
           type(C_PTR) rv
           end function ATK_create_fortran_allocatable_view
        end interface
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        integer(C_LONG), allocatable, intent(IN) :: value
        integer(C_INT) :: lname
        type(dataview) :: rv
        type(C_PTR) :: addr
        integer(C_INT), parameter :: rank = 0
        integer(C_INT), parameter :: itype = ATK_C_LONG_T
    
        lname = len_trim(name)
        call c_loc_allocatable(value, addr)
        rv%voidptr = ATK_create_fortran_allocatable_view(group%voidptr, name, lname, addr, itype, rank)
    end function datagroup_create_allocatable_view_long_scalar
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_allocatable_view_long_1d(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_fortran_allocatable_view(group, name, lname, addr, itype, rank) &
              bind(C,name="ATK_create_fortran_allocatable_view") &
              result(rv)
           use iso_c_binding
           type(C_PTR), value, intent(IN)    :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN) :: lname
           type(C_PTR), value                :: addr
           integer(C_INT), value, intent(IN) :: itype
           integer(C_INT), value, intent(IN) :: rank
           type(C_PTR) rv
           end function ATK_create_fortran_allocatable_view
        end interface
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        integer(C_LONG), allocatable, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(dataview) :: rv
        type(C_PTR) :: addr
        integer(C_INT), parameter :: rank = 1
        integer(C_INT), parameter :: itype = ATK_C_LONG_T
    
        lname = len_trim(name)
        call c_loc_allocatable(value, addr)
        rv%voidptr = ATK_create_fortran_allocatable_view(group%voidptr, name, lname, addr, itype, rank)
    end function datagroup_create_allocatable_view_long_1d
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_allocatable_view_float_scalar(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_fortran_allocatable_view(group, name, lname, addr, itype, rank) &
              bind(C,name="ATK_create_fortran_allocatable_view") &
              result(rv)
           use iso_c_binding
           type(C_PTR), value, intent(IN)    :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN) :: lname
           type(C_PTR), value                :: addr
           integer(C_INT), value, intent(IN) :: itype
           integer(C_INT), value, intent(IN) :: rank
           type(C_PTR) rv
           end function ATK_create_fortran_allocatable_view
        end interface
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        real(C_FLOAT), allocatable, intent(IN) :: value
        integer(C_INT) :: lname
        type(dataview) :: rv
        type(C_PTR) :: addr
        integer(C_INT), parameter :: rank = 0
        integer(C_INT), parameter :: itype = ATK_C_FLOAT_T
    
        lname = len_trim(name)
        call c_loc_allocatable(value, addr)
        rv%voidptr = ATK_create_fortran_allocatable_view(group%voidptr, name, lname, addr, itype, rank)
    end function datagroup_create_allocatable_view_float_scalar
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_allocatable_view_float_1d(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_fortran_allocatable_view(group, name, lname, addr, itype, rank) &
              bind(C,name="ATK_create_fortran_allocatable_view") &
              result(rv)
           use iso_c_binding
           type(C_PTR), value, intent(IN)    :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN) :: lname
           type(C_PTR), value                :: addr
           integer(C_INT), value, intent(IN) :: itype
           integer(C_INT), value, intent(IN) :: rank
           type(C_PTR) rv
           end function ATK_create_fortran_allocatable_view
        end interface
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        real(C_FLOAT), allocatable, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(dataview) :: rv
        type(C_PTR) :: addr
        integer(C_INT), parameter :: rank = 1
        integer(C_INT), parameter :: itype = ATK_C_FLOAT_T
    
        lname = len_trim(name)
        call c_loc_allocatable(value, addr)
        rv%voidptr = ATK_create_fortran_allocatable_view(group%voidptr, name, lname, addr, itype, rank)
    end function datagroup_create_allocatable_view_float_1d
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_allocatable_view_double_scalar(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_fortran_allocatable_view(group, name, lname, addr, itype, rank) &
              bind(C,name="ATK_create_fortran_allocatable_view") &
              result(rv)
           use iso_c_binding
           type(C_PTR), value, intent(IN)    :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN) :: lname
           type(C_PTR), value                :: addr
           integer(C_INT), value, intent(IN) :: itype
           integer(C_INT), value, intent(IN) :: rank
           type(C_PTR) rv
           end function ATK_create_fortran_allocatable_view
        end interface
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        real(C_DOUBLE), allocatable, intent(IN) :: value
        integer(C_INT) :: lname
        type(dataview) :: rv
        type(C_PTR) :: addr
        integer(C_INT), parameter :: rank = 0
        integer(C_INT), parameter :: itype = ATK_C_DOUBLE_T
    
        lname = len_trim(name)
        call c_loc_allocatable(value, addr)
        rv%voidptr = ATK_create_fortran_allocatable_view(group%voidptr, name, lname, addr, itype, rank)
    end function datagroup_create_allocatable_view_double_scalar
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_allocatable_view_double_1d(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_fortran_allocatable_view(group, name, lname, addr, itype, rank) &
              bind(C,name="ATK_create_fortran_allocatable_view") &
              result(rv)
           use iso_c_binding
           type(C_PTR), value, intent(IN)    :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN) :: lname
           type(C_PTR), value                :: addr
           integer(C_INT), value, intent(IN) :: itype
           integer(C_INT), value, intent(IN) :: rank
           type(C_PTR) rv
           end function ATK_create_fortran_allocatable_view
        end interface
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        real(C_DOUBLE), allocatable, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(dataview) :: rv
        type(C_PTR) :: addr
        integer(C_INT), parameter :: rank = 1
        integer(C_INT), parameter :: itype = ATK_C_DOUBLE_T
    
        lname = len_trim(name)
        call c_loc_allocatable(value, addr)
        rv%voidptr = ATK_create_fortran_allocatable_view(group%voidptr, name, lname, addr, itype, rank)
    end function datagroup_create_allocatable_view_double_1d
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_array_view_int_scalar(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_array_view(group, name, lname, addr, type, nitems) result(rv) bind(C,name="ATK_create_array_view")
           use iso_c_binding
           type(C_PTR), value, intent(IN)     :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN)  :: lname
           type(C_PTR), value,     intent(IN) :: addr
           integer(C_INT), value, intent(IN)  :: type
           integer(C_LONG), value, intent(IN) :: nitems
           type(C_PTR) rv
           end function ATK_create_array_view
        end interface
        external :: ATK_C_LOC
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        integer(C_INT), target, intent(IN) :: value
        integer(C_INT) :: lname
        type(dataview) :: rv
        integer(C_LONG) :: nitems
        integer(C_INT), parameter :: type = ATK_C_INT_T
        type(C_PTR) addr
    
        lname = len_trim(name)
        nitems = 1_C_LONG
        call ATK_C_LOC(value, addr)
        rv%voidptr = ATK_create_array_view(group%voidptr, name, lname, addr, type, nitems)
    end function datagroup_create_array_view_int_scalar
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_array_view_int_1d(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_array_view(group, name, lname, addr, type, nitems) result(rv) bind(C,name="ATK_create_array_view")
           use iso_c_binding
           type(C_PTR), value, intent(IN)     :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN)  :: lname
           type(C_PTR), value,     intent(IN) :: addr
           integer(C_INT), value, intent(IN)  :: type
           integer(C_LONG), value, intent(IN) :: nitems
           type(C_PTR) rv
           end function ATK_create_array_view
        end interface
        external :: ATK_C_LOC
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        integer(C_INT), target, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(dataview) :: rv
        integer(C_LONG) :: nitems
        integer(C_INT), parameter :: type = ATK_C_INT_T
        type(C_PTR) addr
    
        lname = len_trim(name)
        nitems = size(value, kind=1_C_LONG)
        call ATK_C_LOC(value, addr)
        rv%voidptr = ATK_create_array_view(group%voidptr, name, lname, addr, type, nitems)
    end function datagroup_create_array_view_int_1d
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_array_view_long_scalar(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_array_view(group, name, lname, addr, type, nitems) result(rv) bind(C,name="ATK_create_array_view")
           use iso_c_binding
           type(C_PTR), value, intent(IN)     :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN)  :: lname
           type(C_PTR), value,     intent(IN) :: addr
           integer(C_INT), value, intent(IN)  :: type
           integer(C_LONG), value, intent(IN) :: nitems
           type(C_PTR) rv
           end function ATK_create_array_view
        end interface
        external :: ATK_C_LOC
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        integer(C_LONG), target, intent(IN) :: value
        integer(C_INT) :: lname
        type(dataview) :: rv
        integer(C_LONG) :: nitems
        integer(C_INT), parameter :: type = ATK_C_LONG_T
        type(C_PTR) addr
    
        lname = len_trim(name)
        nitems = 1_C_LONG
        call ATK_C_LOC(value, addr)
        rv%voidptr = ATK_create_array_view(group%voidptr, name, lname, addr, type, nitems)
    end function datagroup_create_array_view_long_scalar
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_array_view_long_1d(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_array_view(group, name, lname, addr, type, nitems) result(rv) bind(C,name="ATK_create_array_view")
           use iso_c_binding
           type(C_PTR), value, intent(IN)     :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN)  :: lname
           type(C_PTR), value,     intent(IN) :: addr
           integer(C_INT), value, intent(IN)  :: type
           integer(C_LONG), value, intent(IN) :: nitems
           type(C_PTR) rv
           end function ATK_create_array_view
        end interface
        external :: ATK_C_LOC
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        integer(C_LONG), target, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(dataview) :: rv
        integer(C_LONG) :: nitems
        integer(C_INT), parameter :: type = ATK_C_LONG_T
        type(C_PTR) addr
    
        lname = len_trim(name)
        nitems = size(value, kind=1_C_LONG)
        call ATK_C_LOC(value, addr)
        rv%voidptr = ATK_create_array_view(group%voidptr, name, lname, addr, type, nitems)
    end function datagroup_create_array_view_long_1d
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_array_view_float_scalar(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_array_view(group, name, lname, addr, type, nitems) result(rv) bind(C,name="ATK_create_array_view")
           use iso_c_binding
           type(C_PTR), value, intent(IN)     :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN)  :: lname
           type(C_PTR), value,     intent(IN) :: addr
           integer(C_INT), value, intent(IN)  :: type
           integer(C_LONG), value, intent(IN) :: nitems
           type(C_PTR) rv
           end function ATK_create_array_view
        end interface
        external :: ATK_C_LOC
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        real(C_FLOAT), target, intent(IN) :: value
        integer(C_INT) :: lname
        type(dataview) :: rv
        integer(C_LONG) :: nitems
        integer(C_INT), parameter :: type = ATK_C_FLOAT_T
        type(C_PTR) addr
    
        lname = len_trim(name)
        nitems = 1_C_LONG
        call ATK_C_LOC(value, addr)
        rv%voidptr = ATK_create_array_view(group%voidptr, name, lname, addr, type, nitems)
    end function datagroup_create_array_view_float_scalar
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_array_view_float_1d(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_array_view(group, name, lname, addr, type, nitems) result(rv) bind(C,name="ATK_create_array_view")
           use iso_c_binding
           type(C_PTR), value, intent(IN)     :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN)  :: lname
           type(C_PTR), value,     intent(IN) :: addr
           integer(C_INT), value, intent(IN)  :: type
           integer(C_LONG), value, intent(IN) :: nitems
           type(C_PTR) rv
           end function ATK_create_array_view
        end interface
        external :: ATK_C_LOC
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        real(C_FLOAT), target, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(dataview) :: rv
        integer(C_LONG) :: nitems
        integer(C_INT), parameter :: type = ATK_C_FLOAT_T
        type(C_PTR) addr
    
        lname = len_trim(name)
        nitems = size(value, kind=1_C_LONG)
        call ATK_C_LOC(value, addr)
        rv%voidptr = ATK_create_array_view(group%voidptr, name, lname, addr, type, nitems)
    end function datagroup_create_array_view_float_1d
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_array_view_double_scalar(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_array_view(group, name, lname, addr, type, nitems) result(rv) bind(C,name="ATK_create_array_view")
           use iso_c_binding
           type(C_PTR), value, intent(IN)     :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN)  :: lname
           type(C_PTR), value,     intent(IN) :: addr
           integer(C_INT), value, intent(IN)  :: type
           integer(C_LONG), value, intent(IN) :: nitems
           type(C_PTR) rv
           end function ATK_create_array_view
        end interface
        external :: ATK_C_LOC
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        real(C_DOUBLE), target, intent(IN) :: value
        integer(C_INT) :: lname
        type(dataview) :: rv
        integer(C_LONG) :: nitems
        integer(C_INT), parameter :: type = ATK_C_DOUBLE_T
        type(C_PTR) addr
    
        lname = len_trim(name)
        nitems = 1_C_LONG
        call ATK_C_LOC(value, addr)
        rv%voidptr = ATK_create_array_view(group%voidptr, name, lname, addr, type, nitems)
    end function datagroup_create_array_view_double_scalar
    
    ! Generated by genfsidresplicer.py
    function datagroup_create_array_view_double_1d(group, name, value) result(rv)
        use iso_c_binding
        implicit none
    
        interface
           function ATK_create_array_view(group, name, lname, addr, type, nitems) result(rv) bind(C,name="ATK_create_array_view")
           use iso_c_binding
           type(C_PTR), value, intent(IN)     :: group
           character(kind=C_CHAR), intent(IN) :: name(*)
           integer(C_INT), value, intent(IN)  :: lname
           type(C_PTR), value,     intent(IN) :: addr
           integer(C_INT), value, intent(IN)  :: type
           integer(C_LONG), value, intent(IN) :: nitems
           type(C_PTR) rv
           end function ATK_create_array_view
        end interface
        external :: ATK_C_LOC
    
        class(datagroup), intent(IN) :: group
        character(*), intent(IN) :: name
        real(C_DOUBLE), target, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(dataview) :: rv
        integer(C_LONG) :: nitems
        integer(C_INT), parameter :: type = ATK_C_DOUBLE_T
        type(C_PTR) addr
    
        lname = len_trim(name)
        nitems = size(value, kind=1_C_LONG)
        call ATK_C_LOC(value, addr)
        rv%voidptr = ATK_create_array_view(group%voidptr, name, lname, addr, type, nitems)
    end function datagroup_create_array_view_double_1d
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
    
    subroutine databuffer_declare_int(obj, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_INT) :: len
        ! splicer begin class.DataBuffer.method.declare_int
        call atk_databuffer_declare(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataBuffer.method.declare_int
    end subroutine databuffer_declare_int
    
    subroutine databuffer_declare_long(obj, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataBuffer.method.declare_long
        call atk_databuffer_declare(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataBuffer.method.declare_long
    end subroutine databuffer_declare_long
    
    subroutine databuffer_allocate_existing(obj)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        ! splicer begin class.DataBuffer.method.allocate_existing
        call atk_databuffer_allocate_existing(obj%voidptr)
        ! splicer end class.DataBuffer.method.allocate_existing
    end subroutine databuffer_allocate_existing
    
    subroutine databuffer_allocate_int(obj, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_INT) :: len
        ! splicer begin class.DataBuffer.method.allocate_int
        call atk_databuffer_allocate_from_type(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataBuffer.method.allocate_int
    end subroutine databuffer_allocate_int
    
    subroutine databuffer_allocate_long(obj, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataBuffer.method.allocate_long
        call atk_databuffer_allocate_from_type(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataBuffer.method.allocate_long
    end subroutine databuffer_allocate_long
    
    subroutine databuffer_reallocate_int(obj, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_INT) :: len
        ! splicer begin class.DataBuffer.method.reallocate_int
        call atk_databuffer_reallocate(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataBuffer.method.reallocate_int
    end subroutine databuffer_reallocate_int
    
    subroutine databuffer_reallocate_long(obj, type, len)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataBuffer.method.reallocate_long
        call atk_databuffer_reallocate(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataBuffer.method.reallocate_long
    end subroutine databuffer_reallocate_long
    
    subroutine databuffer_set_external_data(obj, external_data)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        type(C_PTR) :: external_data
        ! splicer begin class.DataBuffer.method.set_external_data
        call atk_databuffer_set_external_data(  &
            obj%voidptr,  &
            external_data)
        ! splicer end class.DataBuffer.method.set_external_data
    end subroutine databuffer_set_external_data
    
    function databuffer_is_external(obj) result(rv)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        logical :: rv
        ! splicer begin class.DataBuffer.method.is_external
        rv = atk_databuffer_is_external(obj%voidptr)
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
    
    subroutine databuffer_print(obj)
        use iso_c_binding
        implicit none
        class(databuffer) :: obj
        ! splicer begin class.DataBuffer.method.print
        call atk_databuffer_print(obj%voidptr)
        ! splicer end class.DataBuffer.method.print
    end subroutine databuffer_print
    
    ! splicer begin class.DataBuffer.additional_functions
    ! splicer end class.DataBuffer.additional_functions
    
    subroutine dataview_declare_int(obj, type, len)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_INT) :: len
        ! splicer begin class.DataView.method.declare_int
        call atk_dataview_declare(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataView.method.declare_int
    end subroutine dataview_declare_int
    
    subroutine dataview_declare_long(obj, type, len)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataView.method.declare_long
        call atk_dataview_declare(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataView.method.declare_long
    end subroutine dataview_declare_long
    
    subroutine dataview_allocate_simple(obj)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        ! splicer begin class.DataView.method.allocate_simple
        call atk_dataview_allocate_simple(obj%voidptr)
        ! splicer end class.DataView.method.allocate_simple
    end subroutine dataview_allocate_simple
    
    subroutine dataview_allocate_int(obj, type, len)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_INT) :: len
        ! splicer begin class.DataView.method.allocate_int
        call atk_dataview_allocate_from_type(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataView.method.allocate_int
    end subroutine dataview_allocate_int
    
    subroutine dataview_allocate_long(obj, type, len)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataView.method.allocate_long
        call atk_dataview_allocate_from_type(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataView.method.allocate_long
    end subroutine dataview_allocate_long
    
    subroutine dataview_reallocate_int(obj, type, len)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_INT) :: len
        ! splicer begin class.DataView.method.reallocate_int
        call atk_dataview_reallocate(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataView.method.reallocate_int
    end subroutine dataview_reallocate_int
    
    subroutine dataview_reallocate_long(obj, type, len)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: type
        integer(C_LONG) :: len
        ! splicer begin class.DataView.method.reallocate_long
        call atk_dataview_reallocate(  &
            obj%voidptr,  &
            type,  &
            int(len, C_LONG))
        ! splicer end class.DataView.method.reallocate_long
    end subroutine dataview_reallocate_long
    
    function dataview_has_buffer(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        logical :: rv
        ! splicer begin class.DataView.method.has_buffer
        rv = atk_dataview_has_buffer(obj%voidptr)
        ! splicer end class.DataView.method.has_buffer
    end function dataview_has_buffer
    
    function dataview_is_opaque(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        logical :: rv
        ! splicer begin class.DataView.method.is_opaque
        rv = atk_dataview_is_opaque(obj%voidptr)
        ! splicer end class.DataView.method.is_opaque
    end function dataview_is_opaque
    
    subroutine dataview_get_name(obj, name)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        character(*), intent(OUT) :: name
        type(C_PTR) :: rv
        ! splicer begin class.DataView.method.get_name
        rv = atk_dataview_get_name(obj%voidptr)
        call FccCopyPtr(name, len(name), rv)
        ! splicer end class.DataView.method.get_name
    end subroutine dataview_get_name
    
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
    
    subroutine dataview_set_value_int(obj, value)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: value
        ! splicer begin class.DataView.method.set_value_int
        call atk_dataview_set_value_int(  &
            obj%voidptr,  &
            value)
        ! splicer end class.DataView.method.set_value_int
    end subroutine dataview_set_value_int
    
    subroutine dataview_set_value_long(obj, value)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_LONG) :: value
        ! splicer begin class.DataView.method.set_value_long
        call atk_dataview_set_value_long(  &
            obj%voidptr,  &
            value)
        ! splicer end class.DataView.method.set_value_long
    end subroutine dataview_set_value_long
    
    subroutine dataview_set_value_float(obj, value)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        real(C_FLOAT) :: value
        ! splicer begin class.DataView.method.set_value_float
        call atk_dataview_set_value_float(  &
            obj%voidptr,  &
            value)
        ! splicer end class.DataView.method.set_value_float
    end subroutine dataview_set_value_float
    
    subroutine dataview_set_value_double(obj, value)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        real(C_DOUBLE) :: value
        ! splicer begin class.DataView.method.set_value_double
        call atk_dataview_set_value_double(  &
            obj%voidptr,  &
            value)
        ! splicer end class.DataView.method.set_value_double
    end subroutine dataview_set_value_double
    
    function dataview_get_value_int(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_INT) :: rv
        ! splicer begin class.DataView.method.get_value_int
        rv = atk_dataview_get_value_int(obj%voidptr)
        ! splicer end class.DataView.method.get_value_int
    end function dataview_get_value_int
    
    function dataview_get_value_long(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        integer(C_LONG) :: rv
        ! splicer begin class.DataView.method.get_value_long
        rv = atk_dataview_get_value_long(obj%voidptr)
        ! splicer end class.DataView.method.get_value_long
    end function dataview_get_value_long
    
    function dataview_get_value_float(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        real(C_FLOAT) :: rv
        ! splicer begin class.DataView.method.get_value_float
        rv = atk_dataview_get_value_float(obj%voidptr)
        ! splicer end class.DataView.method.get_value_float
    end function dataview_get_value_float
    
    function dataview_get_value_double(obj) result(rv)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        real(C_DOUBLE) :: rv
        ! splicer begin class.DataView.method.get_value_double
        rv = atk_dataview_get_value_double(obj%voidptr)
        ! splicer end class.DataView.method.get_value_double
    end function dataview_get_value_double
    
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
    
    subroutine dataview_print(obj)
        use iso_c_binding
        implicit none
        class(dataview) :: obj
        ! splicer begin class.DataView.method.print
        call atk_dataview_print(obj%voidptr)
        ! splicer end class.DataView.method.print
    end subroutine dataview_print
    
    ! splicer begin class.DataView.additional_functions
    
    ! Generated by genfsidresplicer.py
    subroutine dataview_get_value_int_scalar_ptr(view, value)
        use iso_c_binding
        implicit none
        class(dataview), intent(IN) :: view
        integer(C_INT), pointer, intent(OUT) :: value
        type(C_PTR) cptr
    
        cptr = view%get_data_pointer()
        call c_f_pointer(cptr, value)
    end subroutine dataview_get_value_int_scalar_ptr
    
    ! Generated by genfsidresplicer.py
    subroutine dataview_get_value_int_1d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(dataview), intent(IN) :: view
        integer(C_INT), pointer, intent(OUT) :: value(:)
        type(C_PTR) cptr
        integer(C_SIZE_T) nelems
    
        cptr = view%get_data_pointer()
        nelems = view%get_number_of_elements()
        call c_f_pointer(cptr, value, [ nelems ])
    end subroutine dataview_get_value_int_1d_ptr
    
    ! Generated by genfsidresplicer.py
    subroutine dataview_get_value_long_scalar_ptr(view, value)
        use iso_c_binding
        implicit none
        class(dataview), intent(IN) :: view
        integer(C_LONG), pointer, intent(OUT) :: value
        type(C_PTR) cptr
    
        cptr = view%get_data_pointer()
        call c_f_pointer(cptr, value)
    end subroutine dataview_get_value_long_scalar_ptr
    
    ! Generated by genfsidresplicer.py
    subroutine dataview_get_value_long_1d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(dataview), intent(IN) :: view
        integer(C_LONG), pointer, intent(OUT) :: value(:)
        type(C_PTR) cptr
        integer(C_SIZE_T) nelems
    
        cptr = view%get_data_pointer()
        nelems = view%get_number_of_elements()
        call c_f_pointer(cptr, value, [ nelems ])
    end subroutine dataview_get_value_long_1d_ptr
    
    ! Generated by genfsidresplicer.py
    subroutine dataview_get_value_float_scalar_ptr(view, value)
        use iso_c_binding
        implicit none
        class(dataview), intent(IN) :: view
        real(C_FLOAT), pointer, intent(OUT) :: value
        type(C_PTR) cptr
    
        cptr = view%get_data_pointer()
        call c_f_pointer(cptr, value)
    end subroutine dataview_get_value_float_scalar_ptr
    
    ! Generated by genfsidresplicer.py
    subroutine dataview_get_value_float_1d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(dataview), intent(IN) :: view
        real(C_FLOAT), pointer, intent(OUT) :: value(:)
        type(C_PTR) cptr
        integer(C_SIZE_T) nelems
    
        cptr = view%get_data_pointer()
        nelems = view%get_number_of_elements()
        call c_f_pointer(cptr, value, [ nelems ])
    end subroutine dataview_get_value_float_1d_ptr
    
    ! Generated by genfsidresplicer.py
    subroutine dataview_get_value_double_scalar_ptr(view, value)
        use iso_c_binding
        implicit none
        class(dataview), intent(IN) :: view
        real(C_DOUBLE), pointer, intent(OUT) :: value
        type(C_PTR) cptr
    
        cptr = view%get_data_pointer()
        call c_f_pointer(cptr, value)
    end subroutine dataview_get_value_double_scalar_ptr
    
    ! Generated by genfsidresplicer.py
    subroutine dataview_get_value_double_1d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(dataview), intent(IN) :: view
        real(C_DOUBLE), pointer, intent(OUT) :: value(:)
        type(C_PTR) cptr
        integer(C_SIZE_T) nelems
    
        cptr = view%get_data_pointer()
        nelems = view%get_number_of_elements()
        call c_f_pointer(cptr, value, [ nelems ])
    end subroutine dataview_get_value_double_1d_ptr
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
    
    ! splicer begin additional_functions
    ! splicer end additional_functions
    
    function datastore_eq(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(datastore), intent(IN) ::a,b
        logical :: rv
        if (c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function datastore_eq
    
    function datastore_ne(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(datastore), intent(IN) ::a,b
        logical :: rv
        if (.not. c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function datastore_ne
    
    function datagroup_eq(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(datagroup), intent(IN) ::a,b
        logical :: rv
        if (c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function datagroup_eq
    
    function datagroup_ne(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(datagroup), intent(IN) ::a,b
        logical :: rv
        if (.not. c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function datagroup_ne
    
    function databuffer_eq(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(databuffer), intent(IN) ::a,b
        logical :: rv
        if (c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function databuffer_eq
    
    function databuffer_ne(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(databuffer), intent(IN) ::a,b
        logical :: rv
        if (.not. c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function databuffer_ne
    
    function dataview_eq(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(dataview), intent(IN) ::a,b
        logical :: rv
        if (c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function dataview_eq
    
    function dataview_ne(a,b) result (rv)
        use iso_c_binding, only: c_associated
        implicit none
        type(dataview), intent(IN) ::a,b
        logical :: rv
        if (.not. c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function dataview_ne

end module sidre_mod

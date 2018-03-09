! wrapfsidre.F
! This is generated code, do not edit
!
! Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
!
! Produced at the Lawrence Livermore National Laboratory
!
! LLNL-CODE-741217
!
! All rights reserved.
!
! This file is part of Axom.
!
! For details about use and distribution, please read axom/LICENSE.
!
!>
!! \file wrapfsidre.F
!! \brief Shroud generated wrapper for Sidre library
!<
! splicer begin file_top
#include "axom/config.hpp"
! splicer end file_top
module axom_sidre
    use iso_c_binding, only : C_PTR
    ! splicer begin module_use
    ! map conduit type names to sidre type names
    use conduit, only : &
        SIDRE_NO_TYPE_ID    => CONDUIT_EMPTY_ID, &
        SIDRE_INT8_ID       => CONDUIT_INT8_ID, &
        SIDRE_INT16_ID      => CONDUIT_INT16_ID, &
        SIDRE_INT32_ID      => CONDUIT_INT32_ID, &
        SIDRE_INT64_ID      => CONDUIT_INT64_ID, &
        SIDRE_UINT8_ID      => CONDUIT_UINT8_ID, &
        SIDRE_UINT16_ID     => CONDUIT_UINT16_ID, &
        SIDRE_UINT32_ID     => CONDUIT_UINT32_ID, &
        SIDRE_UINT64_ID     => CONDUIT_UINT64_ID, &
        SIDRE_FLOAT32_ID    => CONDUIT_FLOAT32_ID, &
        SIDRE_FLOAT64_ID    => CONDUIT_FLOAT64_ID, &
        SIDRE_CHAR8_STR_ID  => CONDUIT_CHAR8_STR_ID, &
        SIDRE_INT_ID        => CONDUIT_INT_ID, &
        SIDRE_UINT_ID       => CONDUIT_UINT_ID, &
        SIDRE_LONG_ID       => CONDUIT_LONG_ID, &
        SIDRE_ULONG_ID      => CONDUIT_ULONG_ID, &
        SIDRE_FLOAT_ID      => CONDUIT_FLOAT_ID, &
        SIDRE_DOUBLE_ID     => CONDUIT_DOUBLE_ID
    use, intrinsic :: iso_c_binding, only : C_LONG
    ! splicer end module_use
    implicit none

    ! splicer begin module_top
    integer, parameter :: MAXNAMESIZE = 128

    integer, parameter :: SIDRE_LENGTH = C_LONG

    integer, parameter :: invalid_index = -1
    ! splicer end module_top

    ! splicer begin class.DataStore.module_top
    ! splicer end class.DataStore.module_top

    type SidreDataStore
        type(C_PTR), private :: voidptr
        ! splicer begin class.DataStore.component_part
        ! splicer end class.DataStore.component_part
    contains
        procedure :: delete => datastore_delete
        procedure :: get_root => datastore_get_root
        procedure :: get_buffer => datastore_get_buffer
        procedure :: create_buffer_empty => datastore_create_buffer_empty
        procedure :: create_buffer_from_type_int => datastore_create_buffer_from_type_int
        procedure :: create_buffer_from_type_long => datastore_create_buffer_from_type_long
        procedure :: destroy_buffer => datastore_destroy_buffer
        procedure :: get_num_buffers => datastore_get_num_buffers
        procedure :: print => datastore_print
        procedure :: get_instance => datastore_get_instance
        procedure :: set_instance => datastore_set_instance
        procedure :: associated => datastore_associated
        generic :: create_buffer => &
            ! splicer begin class.DataStore.generic.create_buffer
            ! splicer end class.DataStore.generic.create_buffer
            create_buffer_empty,  &
            create_buffer_from_type_int,  &
            create_buffer_from_type_long
        ! splicer begin class.DataStore.type_bound_procedure_part
        ! splicer end class.DataStore.type_bound_procedure_part
    end type SidreDataStore

    ! splicer begin class.Group.module_top
    ! splicer end class.Group.module_top

    type SidreGroup
        type(C_PTR), private :: voidptr
        ! splicer begin class.Group.component_part
        ! splicer end class.Group.component_part
    contains
        procedure :: get_index => group_get_index
        procedure :: get_name => group_get_name
        procedure :: get_path => group_get_path
        procedure :: get_path_name => group_get_path_name
        procedure :: get_parent => group_get_parent
        procedure :: get_data_store => group_get_data_store
        procedure :: get_num_views => group_get_num_views
        procedure :: get_num_groups => group_get_num_groups
        procedure :: has_view => group_has_view
        procedure :: has_child_view => group_has_child_view
        procedure :: rename => group_rename
        procedure :: get_view_from_name => group_get_view_from_name
        procedure :: get_view_from_index => group_get_view_from_index
        procedure :: get_view_index => group_get_view_index
        procedure :: get_view_name => group_get_view_name
        procedure :: get_first_valid_view_index => group_get_first_valid_view_index
        procedure :: get_next_valid_view_index => group_get_next_valid_view_index
        procedure :: create_view_and_allocate_nelems_int => group_create_view_and_allocate_nelems_int
        procedure :: create_view_and_allocate_nelems_long => group_create_view_and_allocate_nelems_long
        procedure :: create_view_and_allocate_shape => group_create_view_and_allocate_shape
        procedure :: create_view_scalar_int => group_create_view_scalar_int
        procedure :: create_view_scalar_long => group_create_view_scalar_long
        procedure :: create_view_scalar_float => group_create_view_scalar_float
        procedure :: create_view_scalar_double => group_create_view_scalar_double
        procedure :: create_view_string => group_create_view_string
        procedure :: create_view_empty => group_create_view_empty
        procedure :: create_view_from_type_int => group_create_view_from_type_int
        procedure :: create_view_from_type_long => group_create_view_from_type_long
        procedure :: create_view_from_type_and_buffer_int => group_create_view_from_type_and_buffer_int
        procedure :: create_view_from_type_and_buffer_long => group_create_view_from_type_and_buffer_long
        procedure :: create_view_from_type_external_int => group_create_view_from_type_external_int
        procedure :: create_view_from_type_external_long => group_create_view_from_type_external_long
        procedure :: create_view_from_shape => group_create_view_from_shape
        procedure :: create_view_from_shape_and_buffer => group_create_view_from_shape_and_buffer
        procedure :: create_view_from_shape_external => group_create_view_from_shape_external
        procedure :: create_view_into_buffer => group_create_view_into_buffer
        procedure :: create_view_external => group_create_view_external
        procedure :: destroy_view => group_destroy_view
        procedure :: destroy_view_and_data_name => group_destroy_view_and_data_name
        procedure :: destroy_view_and_data_index => group_destroy_view_and_data_index
        procedure :: move_view => group_move_view
        procedure :: copy_view => group_copy_view
        procedure :: has_group => group_has_group
        procedure :: has_child_group => group_has_child_group
        procedure :: get_group_from_name => group_get_group_from_name
        procedure :: get_group_from_index => group_get_group_from_index
        procedure :: get_group_index => group_get_group_index
        procedure :: get_group_name => group_get_group_name
        procedure :: get_first_valid_group_index => group_get_first_valid_group_index
        procedure :: get_next_valid_group_index => group_get_next_valid_group_index
        procedure :: create_group => group_create_group
        procedure :: destroy_group_name => group_destroy_group_name
        procedure :: destroy_group_index => group_destroy_group_index
        procedure :: move_group => group_move_group
        procedure :: print => group_print
        procedure :: is_equivalent_to => group_is_equivalent_to
        procedure :: save => group_save
        procedure :: load_0 => group_load_0
        procedure :: load_1 => group_load_1
        procedure :: load_external_data => group_load_external_data
        procedure :: get_instance => group_get_instance
        procedure :: set_instance => group_set_instance
        procedure :: associated => group_associated
        generic :: create_view => &
            ! splicer begin class.Group.generic.create_view
            ! splicer end class.Group.generic.create_view
            create_view_empty,  &
            create_view_from_type_int,  &
            create_view_from_type_long,  &
            create_view_from_type_and_buffer_int,  &
            create_view_from_type_and_buffer_long,  &
            create_view_from_type_external_int,  &
            create_view_from_type_external_long,  &
            create_view_from_shape,  &
            create_view_from_shape_and_buffer,  &
            create_view_from_shape_external,  &
            create_view_into_buffer,  &
            create_view_external
        generic :: create_view_and_allocate => &
            ! splicer begin class.Group.generic.create_view_and_allocate
            ! splicer end class.Group.generic.create_view_and_allocate
            create_view_and_allocate_nelems_int,  &
            create_view_and_allocate_nelems_long,  &
            create_view_and_allocate_shape
        generic :: create_view_scalar => &
            ! splicer begin class.Group.generic.create_view_scalar
            ! splicer end class.Group.generic.create_view_scalar
            create_view_scalar_int,  &
            create_view_scalar_long,  &
            create_view_scalar_float,  &
            create_view_scalar_double
        generic :: destroy_group => &
            ! splicer begin class.Group.generic.destroy_group
            ! splicer end class.Group.generic.destroy_group
            destroy_group_name,  &
            destroy_group_index
        generic :: destroy_view_and_data => &
            ! splicer begin class.Group.generic.destroy_view_and_data
            ! splicer end class.Group.generic.destroy_view_and_data
            destroy_view_and_data_name,  &
            destroy_view_and_data_index
        generic :: get_group => &
            ! splicer begin class.Group.generic.get_group
            ! splicer end class.Group.generic.get_group
            get_group_from_name,  &
            get_group_from_index
        generic :: get_view => &
            ! splicer begin class.Group.generic.get_view
            ! splicer end class.Group.generic.get_view
            get_view_from_name,  &
            get_view_from_index
        generic :: load => &
            ! splicer begin class.Group.generic.load
            ! splicer end class.Group.generic.load
            load_0,  &
            load_1
        ! splicer begin class.Group.type_bound_procedure_part
        procedure :: get_scalar_int => group_get_scalar_int
        procedure :: get_scalar_long => group_get_scalar_long
        procedure :: get_scalar_float => group_get_scalar_float
        procedure :: get_scalar_double => group_get_scalar_double
        generic :: get_scalar => &
            get_scalar_int,  &
            get_scalar_long,  &
            get_scalar_float,  &
            get_scalar_double
        procedure :: set_scalar_int => group_set_scalar_int
        procedure :: set_scalar_long => group_set_scalar_long
        procedure :: set_scalar_float => group_set_scalar_float
        procedure :: set_scalar_double => group_set_scalar_double
        generic :: set_scalar => &
            set_scalar_int,  &
            set_scalar_long,  &
            set_scalar_float,  &
            set_scalar_double
        procedure :: create_array_view_int_scalar => group_create_array_view_int_scalar
        procedure :: create_array_view_int_1d => group_create_array_view_int_1d
        procedure :: create_array_view_int_2d => group_create_array_view_int_2d
        procedure :: create_array_view_int_3d => group_create_array_view_int_3d
        procedure :: create_array_view_int_4d => group_create_array_view_int_4d
        procedure :: create_array_view_long_scalar => group_create_array_view_long_scalar
        procedure :: create_array_view_long_1d => group_create_array_view_long_1d
        procedure :: create_array_view_long_2d => group_create_array_view_long_2d
        procedure :: create_array_view_long_3d => group_create_array_view_long_3d
        procedure :: create_array_view_long_4d => group_create_array_view_long_4d
        procedure :: create_array_view_float_scalar => group_create_array_view_float_scalar
        procedure :: create_array_view_float_1d => group_create_array_view_float_1d
        procedure :: create_array_view_float_2d => group_create_array_view_float_2d
        procedure :: create_array_view_float_3d => group_create_array_view_float_3d
        procedure :: create_array_view_float_4d => group_create_array_view_float_4d
        procedure :: create_array_view_double_scalar => group_create_array_view_double_scalar
        procedure :: create_array_view_double_1d => group_create_array_view_double_1d
        procedure :: create_array_view_double_2d => group_create_array_view_double_2d
        procedure :: create_array_view_double_3d => group_create_array_view_double_3d
        procedure :: create_array_view_double_4d => group_create_array_view_double_4d
        generic :: create_array_view => &
            create_array_view_int_scalar,  &
            create_array_view_int_1d,  &
            create_array_view_int_2d,  &
            create_array_view_int_3d,  &
            create_array_view_int_4d,  &
            create_array_view_long_scalar,  &
            create_array_view_long_1d,  &
            create_array_view_long_2d,  &
            create_array_view_long_3d,  &
            create_array_view_long_4d,  &
            create_array_view_float_scalar,  &
            create_array_view_float_1d,  &
            create_array_view_float_2d,  &
            create_array_view_float_3d,  &
            create_array_view_float_4d,  &
            create_array_view_double_scalar,  &
            create_array_view_double_1d,  &
            create_array_view_double_2d,  &
            create_array_view_double_3d,  &
            create_array_view_double_4d
        procedure :: set_array_data_ptr_int_scalar => group_set_array_data_ptr_int_scalar
        procedure :: set_array_data_ptr_int_1d => group_set_array_data_ptr_int_1d
        procedure :: set_array_data_ptr_int_2d => group_set_array_data_ptr_int_2d
        procedure :: set_array_data_ptr_int_3d => group_set_array_data_ptr_int_3d
        procedure :: set_array_data_ptr_int_4d => group_set_array_data_ptr_int_4d
        procedure :: set_array_data_ptr_long_scalar => group_set_array_data_ptr_long_scalar
        procedure :: set_array_data_ptr_long_1d => group_set_array_data_ptr_long_1d
        procedure :: set_array_data_ptr_long_2d => group_set_array_data_ptr_long_2d
        procedure :: set_array_data_ptr_long_3d => group_set_array_data_ptr_long_3d
        procedure :: set_array_data_ptr_long_4d => group_set_array_data_ptr_long_4d
        procedure :: set_array_data_ptr_float_scalar => group_set_array_data_ptr_float_scalar
        procedure :: set_array_data_ptr_float_1d => group_set_array_data_ptr_float_1d
        procedure :: set_array_data_ptr_float_2d => group_set_array_data_ptr_float_2d
        procedure :: set_array_data_ptr_float_3d => group_set_array_data_ptr_float_3d
        procedure :: set_array_data_ptr_float_4d => group_set_array_data_ptr_float_4d
        procedure :: set_array_data_ptr_double_scalar => group_set_array_data_ptr_double_scalar
        procedure :: set_array_data_ptr_double_1d => group_set_array_data_ptr_double_1d
        procedure :: set_array_data_ptr_double_2d => group_set_array_data_ptr_double_2d
        procedure :: set_array_data_ptr_double_3d => group_set_array_data_ptr_double_3d
        procedure :: set_array_data_ptr_double_4d => group_set_array_data_ptr_double_4d
        generic :: set_array_data_ptr => &
            set_array_data_ptr_int_scalar,  &
            set_array_data_ptr_int_1d,  &
            set_array_data_ptr_int_2d,  &
            set_array_data_ptr_int_3d,  &
            set_array_data_ptr_int_4d,  &
            set_array_data_ptr_long_scalar,  &
            set_array_data_ptr_long_1d,  &
            set_array_data_ptr_long_2d,  &
            set_array_data_ptr_long_3d,  &
            set_array_data_ptr_long_4d,  &
            set_array_data_ptr_float_scalar,  &
            set_array_data_ptr_float_1d,  &
            set_array_data_ptr_float_2d,  &
            set_array_data_ptr_float_3d,  &
            set_array_data_ptr_float_4d,  &
            set_array_data_ptr_double_scalar,  &
            set_array_data_ptr_double_1d,  &
            set_array_data_ptr_double_2d,  &
            set_array_data_ptr_double_3d,  &
            set_array_data_ptr_double_4d
        procedure :: get_string => group_get_string
        procedure :: set_string => group_set_string
        ! splicer end class.Group.type_bound_procedure_part
    end type SidreGroup

    ! splicer begin class.Buffer.module_top
    ! splicer end class.Buffer.module_top

    type SidreBuffer
        type(C_PTR), private :: voidptr
        ! splicer begin class.Buffer.component_part
        ! splicer end class.Buffer.component_part
    contains
        procedure :: get_index => buffer_get_index
        procedure :: get_num_views => buffer_get_num_views
        procedure :: describe_int => buffer_describe_int
        procedure :: describe_long => buffer_describe_long
        procedure :: allocate_existing => buffer_allocate_existing
        procedure :: allocate_from_type_int => buffer_allocate_from_type_int
        procedure :: allocate_from_type_long => buffer_allocate_from_type_long
        procedure :: reallocate_int => buffer_reallocate_int
        procedure :: reallocate_long => buffer_reallocate_long
        procedure :: get_void_ptr => buffer_get_void_ptr
        procedure :: get_type_id => buffer_get_type_id
        procedure :: get_num_elements => buffer_get_num_elements
        procedure :: get_total_bytes => buffer_get_total_bytes
        procedure :: get_bytes_per_element => buffer_get_bytes_per_element
        procedure :: print => buffer_print
        procedure :: get_instance => buffer_get_instance
        procedure :: set_instance => buffer_set_instance
        procedure :: associated => buffer_associated
        generic :: allocate => &
            ! splicer begin class.Buffer.generic.allocate
            ! splicer end class.Buffer.generic.allocate
            allocate_existing,  &
            allocate_from_type_int,  &
            allocate_from_type_long
        generic :: describe => &
            ! splicer begin class.Buffer.generic.describe
            ! splicer end class.Buffer.generic.describe
            describe_int,  &
            describe_long
        generic :: reallocate => &
            ! splicer begin class.Buffer.generic.reallocate
            ! splicer end class.Buffer.generic.reallocate
            reallocate_int,  &
            reallocate_long
        ! splicer begin class.Buffer.type_bound_procedure_part
        ! splicer end class.Buffer.type_bound_procedure_part
    end type SidreBuffer

    ! splicer begin class.View.module_top
    ! splicer end class.View.module_top

    type SidreView
        type(C_PTR), private :: voidptr
        ! splicer begin class.View.component_part
        ! splicer end class.View.component_part
    contains
        procedure :: allocate_simple => view_allocate_simple
        procedure :: allocate_from_type_int => view_allocate_from_type_int
        procedure :: allocate_from_type_long => view_allocate_from_type_long
        procedure :: reallocate_int => view_reallocate_int
        procedure :: reallocate_long => view_reallocate_long
        procedure :: apply_0 => view_apply_0
        procedure :: attach_buffer_only => view_attach_buffer_only
        procedure :: attach_buffer_type_int => view_attach_buffer_type_int
        procedure :: attach_buffer_type_long => view_attach_buffer_type_long
        procedure :: attach_buffer_shape => view_attach_buffer_shape
        procedure :: apply_nelems => view_apply_nelems
        procedure :: apply_nelems_offset => view_apply_nelems_offset
        procedure :: apply_nelems_offset_stride => view_apply_nelems_offset_stride
        procedure :: apply_type_nelems => view_apply_type_nelems
        procedure :: apply_type_nelems_offset => view_apply_type_nelems_offset
        procedure :: apply_type_nelems_offset_stride => view_apply_type_nelems_offset_stride
        procedure :: apply_type_shape => view_apply_type_shape
        procedure :: has_buffer => view_has_buffer
        procedure :: is_external => view_is_external
        procedure :: is_allocated => view_is_allocated
        procedure :: is_applied => view_is_applied
        procedure :: is_described => view_is_described
        procedure :: is_empty => view_is_empty
        procedure :: is_opaque => view_is_opaque
        procedure :: is_scalar => view_is_scalar
        procedure :: is_string => view_is_string
        procedure :: get_index => view_get_index
        procedure :: get_name => view_get_name
        procedure :: get_path => view_get_path
        procedure :: get_path_name => view_get_path_name
        procedure :: get_buffer => view_get_buffer
        procedure :: get_void_ptr => view_get_void_ptr
        procedure :: set_scalar_int => view_set_scalar_int
        procedure :: set_scalar_long => view_set_scalar_long
        procedure :: set_scalar_float => view_set_scalar_float
        procedure :: set_scalar_double => view_set_scalar_double
        procedure :: set_external_data_ptr_only => view_set_external_data_ptr_only
        procedure :: set_external_data_ptr_type_int => view_set_external_data_ptr_type_int
        procedure :: set_external_data_ptr_type_long => view_set_external_data_ptr_type_long
        procedure :: set_string => view_set_string
        procedure :: set_external_data_ptr_shape => view_set_external_data_ptr_shape
        procedure :: get_string => view_get_string
        procedure :: get_data_int => view_get_data_int
        procedure :: get_data_long => view_get_data_long
        procedure :: get_data_float => view_get_data_float
        procedure :: get_data_double => view_get_data_double
        procedure :: get_owning_group => view_get_owning_group
        procedure :: get_type_id => view_get_type_id
        procedure :: get_total_bytes => view_get_total_bytes
        procedure :: get_bytes_per_element => view_get_bytes_per_element
        procedure :: get_num_elements => view_get_num_elements
        procedure :: get_offset => view_get_offset
        procedure :: get_stride => view_get_stride
        procedure :: get_num_dimensions => view_get_num_dimensions
        procedure :: get_shape => view_get_shape
        procedure :: rename => view_rename
        procedure :: print => view_print
        procedure :: get_instance => view_get_instance
        procedure :: set_instance => view_set_instance
        procedure :: associated => view_associated
        generic :: allocate => &
            ! splicer begin class.View.generic.allocate
            ! splicer end class.View.generic.allocate
            allocate_simple,  &
            allocate_from_type_int,  &
            allocate_from_type_long
        generic :: apply => &
            ! splicer begin class.View.generic.apply
            ! splicer end class.View.generic.apply
            apply_0,  &
            apply_nelems,  &
            apply_nelems_offset,  &
            apply_nelems_offset_stride,  &
            apply_type_nelems,  &
            apply_type_nelems_offset,  &
            apply_type_nelems_offset_stride,  &
            apply_type_shape
        generic :: attach_buffer => &
            ! splicer begin class.View.generic.attach_buffer
            ! splicer end class.View.generic.attach_buffer
            attach_buffer_only,  &
            attach_buffer_type_int,  &
            attach_buffer_type_long,  &
            attach_buffer_shape
        generic :: reallocate => &
            ! splicer begin class.View.generic.reallocate
            ! splicer end class.View.generic.reallocate
            reallocate_int,  &
            reallocate_long
        generic :: set_external_data_ptr => &
            ! splicer begin class.View.generic.set_external_data_ptr
            ! splicer end class.View.generic.set_external_data_ptr
            set_external_data_ptr_only,  &
            set_external_data_ptr_type_int,  &
            set_external_data_ptr_type_long,  &
            set_external_data_ptr_shape
        generic :: set_scalar => &
            ! splicer begin class.View.generic.set_scalar
            ! splicer end class.View.generic.set_scalar
            set_scalar_int,  &
            set_scalar_long,  &
            set_scalar_float,  &
            set_scalar_double
        ! splicer begin class.View.type_bound_procedure_part
        procedure :: get_data_int_scalar_ptr => view_get_data_int_scalar_ptr
        procedure :: get_data_int_1d_ptr => view_get_data_int_1d_ptr
        procedure :: get_data_int_2d_ptr => view_get_data_int_2d_ptr
        procedure :: get_data_int_3d_ptr => view_get_data_int_3d_ptr
        procedure :: get_data_int_4d_ptr => view_get_data_int_4d_ptr
        procedure :: get_data_long_scalar_ptr => view_get_data_long_scalar_ptr
        procedure :: get_data_long_1d_ptr => view_get_data_long_1d_ptr
        procedure :: get_data_long_2d_ptr => view_get_data_long_2d_ptr
        procedure :: get_data_long_3d_ptr => view_get_data_long_3d_ptr
        procedure :: get_data_long_4d_ptr => view_get_data_long_4d_ptr
        procedure :: get_data_float_scalar_ptr => view_get_data_float_scalar_ptr
        procedure :: get_data_float_1d_ptr => view_get_data_float_1d_ptr
        procedure :: get_data_float_2d_ptr => view_get_data_float_2d_ptr
        procedure :: get_data_float_3d_ptr => view_get_data_float_3d_ptr
        procedure :: get_data_float_4d_ptr => view_get_data_float_4d_ptr
        procedure :: get_data_double_scalar_ptr => view_get_data_double_scalar_ptr
        procedure :: get_data_double_1d_ptr => view_get_data_double_1d_ptr
        procedure :: get_data_double_2d_ptr => view_get_data_double_2d_ptr
        procedure :: get_data_double_3d_ptr => view_get_data_double_3d_ptr
        procedure :: get_data_double_4d_ptr => view_get_data_double_4d_ptr
        generic :: get_data => &
            get_data_int_scalar_ptr,  &
            get_data_int_1d_ptr,  &
            get_data_int_2d_ptr,  &
            get_data_int_3d_ptr,  &
            get_data_int_4d_ptr,  &
            get_data_long_scalar_ptr,  &
            get_data_long_1d_ptr,  &
            get_data_long_2d_ptr,  &
            get_data_long_3d_ptr,  &
            get_data_long_4d_ptr,  &
            get_data_float_scalar_ptr,  &
            get_data_float_1d_ptr,  &
            get_data_float_2d_ptr,  &
            get_data_float_3d_ptr,  &
            get_data_float_4d_ptr,  &
            get_data_double_scalar_ptr,  &
            get_data_double_1d_ptr,  &
            get_data_double_2d_ptr,  &
            get_data_double_3d_ptr,  &
            get_data_double_4d_ptr
        procedure :: set_array_data_ptr_int_scalar => view_set_array_data_ptr_int_scalar
        procedure :: set_array_data_ptr_int_1d => view_set_array_data_ptr_int_1d
        procedure :: set_array_data_ptr_int_2d => view_set_array_data_ptr_int_2d
        procedure :: set_array_data_ptr_int_3d => view_set_array_data_ptr_int_3d
        procedure :: set_array_data_ptr_int_4d => view_set_array_data_ptr_int_4d
        procedure :: set_array_data_ptr_long_scalar => view_set_array_data_ptr_long_scalar
        procedure :: set_array_data_ptr_long_1d => view_set_array_data_ptr_long_1d
        procedure :: set_array_data_ptr_long_2d => view_set_array_data_ptr_long_2d
        procedure :: set_array_data_ptr_long_3d => view_set_array_data_ptr_long_3d
        procedure :: set_array_data_ptr_long_4d => view_set_array_data_ptr_long_4d
        procedure :: set_array_data_ptr_float_scalar => view_set_array_data_ptr_float_scalar
        procedure :: set_array_data_ptr_float_1d => view_set_array_data_ptr_float_1d
        procedure :: set_array_data_ptr_float_2d => view_set_array_data_ptr_float_2d
        procedure :: set_array_data_ptr_float_3d => view_set_array_data_ptr_float_3d
        procedure :: set_array_data_ptr_float_4d => view_set_array_data_ptr_float_4d
        procedure :: set_array_data_ptr_double_scalar => view_set_array_data_ptr_double_scalar
        procedure :: set_array_data_ptr_double_1d => view_set_array_data_ptr_double_1d
        procedure :: set_array_data_ptr_double_2d => view_set_array_data_ptr_double_2d
        procedure :: set_array_data_ptr_double_3d => view_set_array_data_ptr_double_3d
        procedure :: set_array_data_ptr_double_4d => view_set_array_data_ptr_double_4d
        generic :: set_array_data_ptr => &
            set_array_data_ptr_int_scalar,  &
            set_array_data_ptr_int_1d,  &
            set_array_data_ptr_int_2d,  &
            set_array_data_ptr_int_3d,  &
            set_array_data_ptr_int_4d,  &
            set_array_data_ptr_long_scalar,  &
            set_array_data_ptr_long_1d,  &
            set_array_data_ptr_long_2d,  &
            set_array_data_ptr_long_3d,  &
            set_array_data_ptr_long_4d,  &
            set_array_data_ptr_float_scalar,  &
            set_array_data_ptr_float_1d,  &
            set_array_data_ptr_float_2d,  &
            set_array_data_ptr_float_3d,  &
            set_array_data_ptr_float_4d,  &
            set_array_data_ptr_double_scalar,  &
            set_array_data_ptr_double_1d,  &
            set_array_data_ptr_double_2d,  &
            set_array_data_ptr_double_3d,  &
            set_array_data_ptr_double_4d
        ! splicer end class.View.type_bound_procedure_part
    end type SidreView

    interface operator (.eq.)
        module procedure datastore_eq
        module procedure group_eq
        module procedure buffer_eq
        module procedure view_eq
    end interface

    interface operator (.ne.)
        module procedure datastore_ne
        module procedure group_ne
        module procedure buffer_ne
        module procedure view_ne
    end interface

    interface

        function c_datastore_new() &
                result(SHT_rv) &
                bind(C, name="SIDRE_datastore_new")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR) :: SHT_rv
        end function c_datastore_new

        subroutine c_datastore_delete(self) &
                bind(C, name="SIDRE_datastore_delete")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine c_datastore_delete

        function c_datastore_get_root(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_datastore_get_root")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: SHT_rv
        end function c_datastore_get_root

        function c_datastore_get_buffer(self, idx) &
                result(SHT_rv) &
                bind(C, name="SIDRE_datastore_get_buffer")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            type(C_PTR) :: SHT_rv
        end function c_datastore_get_buffer

        function c_datastore_create_buffer_empty(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_datastore_create_buffer_empty")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: SHT_rv
        end function c_datastore_create_buffer_empty

        function c_datastore_create_buffer_from_type(self, type, &
                num_elems) &
                result(SHT_rv) &
                bind(C, name="SIDRE_datastore_create_buffer_from_type")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            type(C_PTR) :: SHT_rv
        end function c_datastore_create_buffer_from_type

        subroutine c_datastore_destroy_buffer(self, id) &
                bind(C, name="SIDRE_datastore_destroy_buffer")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: id
        end subroutine c_datastore_destroy_buffer

        pure function c_datastore_get_num_buffers(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_datastore_get_num_buffers")
            use iso_c_binding, only : C_PTR, C_SIZE_T
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: SHT_rv
        end function c_datastore_get_num_buffers

        subroutine c_datastore_print(self) &
                bind(C, name="SIDRE_datastore_print")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine c_datastore_print

        ! splicer begin class.DataStore.additional_interfaces
        ! splicer end class.DataStore.additional_interfaces

        function c_group_get_index(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_index")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: SHT_rv
        end function c_group_get_index

        pure function c_group_get_name(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_name")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) SHT_rv
        end function c_group_get_name

        subroutine c_group_get_name_bufferify(self, SHF_rv, NSHF_rv) &
                bind(C, name="SIDRE_group_get_name_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(OUT) :: SHF_rv(*)
            integer(C_INT), value, intent(IN) :: NSHF_rv
        end subroutine c_group_get_name_bufferify

        subroutine c_group_get_path_bufferify(self, SHF_rv, NSHF_rv) &
                bind(C, name="SIDRE_group_get_path_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(OUT) :: SHF_rv(*)
            integer(C_INT), value, intent(IN) :: NSHF_rv
        end subroutine c_group_get_path_bufferify

        subroutine c_group_get_path_name_bufferify(self, SHF_rv, &
                NSHF_rv) &
                bind(C, name="SIDRE_group_get_path_name_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(OUT) :: SHF_rv(*)
            integer(C_INT), value, intent(IN) :: NSHF_rv
        end subroutine c_group_get_path_name_bufferify

        pure function c_group_get_parent(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_parent")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: SHT_rv
        end function c_group_get_parent

        pure function c_group_get_data_store(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_data_store")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: SHT_rv
        end function c_group_get_data_store

        pure function c_group_get_num_views(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_num_views")
            use iso_c_binding, only : C_PTR, C_SIZE_T
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: SHT_rv
        end function c_group_get_num_views

        pure function c_group_get_num_groups(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_num_groups")
            use iso_c_binding, only : C_PTR, C_SIZE_T
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: SHT_rv
        end function c_group_get_num_groups

        pure function c_group_has_view(self, path) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_has_view")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            logical(C_BOOL) :: SHT_rv
        end function c_group_has_view

        pure function c_group_has_view_bufferify(self, path, Lpath) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_has_view_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            logical(C_BOOL) :: SHT_rv
        end function c_group_has_view_bufferify

        pure function c_group_has_child_view(self, name) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_has_child_view")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            logical(C_BOOL) :: SHT_rv
        end function c_group_has_child_view

        pure function c_group_has_child_view_bufferify(self, name, &
                Lname) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_has_child_view_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            logical(C_BOOL) :: SHT_rv
        end function c_group_has_child_view_bufferify

        function c_group_rename(self, new_name) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_rename")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: new_name(*)
            logical(C_BOOL) :: SHT_rv
        end function c_group_rename

        function c_group_rename_bufferify(self, new_name, Lnew_name) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_rename_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: new_name(*)
            integer(C_INT), value, intent(IN) :: Lnew_name
            logical(C_BOOL) :: SHT_rv
        end function c_group_rename_bufferify

        function c_group_get_view_from_name(self, path) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_view_from_name")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            type(C_PTR) :: SHT_rv
        end function c_group_get_view_from_name

        function c_group_get_view_from_name_bufferify(self, path, Lpath) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_view_from_name_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            type(C_PTR) :: SHT_rv
        end function c_group_get_view_from_name_bufferify

        function c_group_get_view_from_index(self, idx) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_view_from_index")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            type(C_PTR) :: SHT_rv
        end function c_group_get_view_from_index

        pure function c_group_get_view_index(self, name) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_view_index")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT) :: SHT_rv
        end function c_group_get_view_index

        pure function c_group_get_view_index_bufferify(self, name, &
                Lname) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_view_index_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            integer(C_INT) :: SHT_rv
        end function c_group_get_view_index_bufferify

        pure function c_group_get_view_name(self, idx) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_view_name")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            type(C_PTR) SHT_rv
        end function c_group_get_view_name

        subroutine c_group_get_view_name_bufferify(self, idx, SHF_rv, &
                NSHF_rv) &
                bind(C, name="SIDRE_group_get_view_name_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            character(kind=C_CHAR), intent(OUT) :: SHF_rv(*)
            integer(C_INT), value, intent(IN) :: NSHF_rv
        end subroutine c_group_get_view_name_bufferify

        pure function c_group_get_first_valid_view_index(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_first_valid_view_index")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: SHT_rv
        end function c_group_get_first_valid_view_index

        pure function c_group_get_next_valid_view_index(self, idx) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_next_valid_view_index")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            integer(C_INT) :: SHT_rv
        end function c_group_get_next_valid_view_index

        function c_group_create_view_and_allocate_nelems(self, path, &
                type, num_elems) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_and_allocate_nelems")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_and_allocate_nelems

        function c_group_create_view_and_allocate_nelems_bufferify(self, &
                path, Lpath, type, num_elems) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_and_allocate_nelems_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_and_allocate_nelems_bufferify

        function c_group_create_view_and_allocate_shape(self, path, &
                type, ndims, shape) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_and_allocate_shape")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: type
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_LONG), intent(IN) :: shape(*)
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_and_allocate_shape

        function c_group_create_view_and_allocate_shape_bufferify(self, &
                path, Lpath, type, ndims, shape) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_and_allocate_shape_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            integer(C_INT), value, intent(IN) :: type
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_LONG), intent(IN) :: shape(*)
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_and_allocate_shape_bufferify

        function c_group_create_view_scalar_int(self, path, value) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_scalar_int")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: value
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_scalar_int

        function c_group_create_view_scalar_int_bufferify(self, path, &
                Lpath, value) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_scalar_int_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            integer(C_INT), value, intent(IN) :: value
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_scalar_int_bufferify

        function c_group_create_view_scalar_long(self, path, value) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_scalar_long")
            use iso_c_binding, only : C_CHAR, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_LONG), value, intent(IN) :: value
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_scalar_long

        function c_group_create_view_scalar_long_bufferify(self, path, &
                Lpath, value) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_scalar_long_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            integer(C_LONG), value, intent(IN) :: value
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_scalar_long_bufferify

        function c_group_create_view_scalar_float(self, path, value) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_scalar_float")
            use iso_c_binding, only : C_CHAR, C_FLOAT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            real(C_FLOAT), value, intent(IN) :: value
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_scalar_float

        function c_group_create_view_scalar_float_bufferify(self, path, &
                Lpath, value) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_scalar_float_bufferify")
            use iso_c_binding, only : C_CHAR, C_FLOAT, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            real(C_FLOAT), value, intent(IN) :: value
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_scalar_float_bufferify

        function c_group_create_view_scalar_double(self, path, value) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_scalar_double")
            use iso_c_binding, only : C_CHAR, C_DOUBLE, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            real(C_DOUBLE), value, intent(IN) :: value
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_scalar_double

        function c_group_create_view_scalar_double_bufferify(self, path, &
                Lpath, value) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_scalar_double_bufferify")
            use iso_c_binding, only : C_CHAR, C_DOUBLE, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            real(C_DOUBLE), value, intent(IN) :: value
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_scalar_double_bufferify

        function c_group_create_view_string(self, path, value) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_string")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            character(kind=C_CHAR), intent(IN) :: value(*)
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_string

        function c_group_create_view_string_bufferify(self, path, Lpath, &
                value, Lvalue) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_string_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            character(kind=C_CHAR), intent(IN) :: value(*)
            integer(C_INT), value, intent(IN) :: Lvalue
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_string_bufferify

        function c_group_create_view_empty(self, path) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_empty")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_empty

        function c_group_create_view_empty_bufferify(self, path, Lpath) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_empty_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_empty_bufferify

        function c_group_create_view_from_type(self, path, type, &
                num_elems) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_from_type")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_from_type

        function c_group_create_view_from_type_bufferify(self, path, &
                Lpath, type, num_elems) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_from_type_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_from_type_bufferify

        function c_group_create_view_from_type_and_buffer(self, path, &
                type, num_elems, buff) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_from_type_and_buffer")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            type(C_PTR), value, intent(IN) :: buff
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_from_type_and_buffer

        function c_group_create_view_from_type_and_buffer_bufferify( &
                self, path, Lpath, type, num_elems, buff) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_from_type_and_buffer_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            type(C_PTR), value, intent(IN) :: buff
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_from_type_and_buffer_bufferify

        function c_group_create_view_from_type_external(self, path, &
                type, num_elems, external_ptr) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_from_type_external")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            type(C_PTR), value, intent(IN) :: external_ptr
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_from_type_external

        function c_group_create_view_from_type_external_bufferify(self, &
                path, Lpath, type, num_elems, external_ptr) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_from_type_external_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            type(C_PTR), value, intent(IN) :: external_ptr
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_from_type_external_bufferify

        function c_group_create_view_from_shape(self, path, type, ndims, &
                shape) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_from_shape")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: type
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_LONG), intent(IN) :: shape(*)
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_from_shape

        function c_group_create_view_from_shape_bufferify(self, path, &
                Lpath, type, ndims, shape) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_from_shape_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            integer(C_INT), value, intent(IN) :: type
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_LONG), intent(IN) :: shape(*)
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_from_shape_bufferify

        function c_group_create_view_from_shape_and_buffer(self, path, &
                type, ndims, shape, buff) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_from_shape_and_buffer")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: type
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_LONG), intent(IN) :: shape(*)
            type(C_PTR), value, intent(IN) :: buff
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_from_shape_and_buffer

        function c_group_create_view_from_shape_and_buffer_bufferify( &
                self, path, Lpath, type, ndims, shape, buff) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_from_shape_and_buffer_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            integer(C_INT), value, intent(IN) :: type
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_LONG), intent(IN) :: shape(*)
            type(C_PTR), value, intent(IN) :: buff
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_from_shape_and_buffer_bufferify

        function c_group_create_view_from_shape_external(self, path, &
                type, ndims, shape, external_ptr) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_from_shape_external")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: type
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_LONG), intent(IN) :: shape(*)
            type(C_PTR), value, intent(IN) :: external_ptr
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_from_shape_external

        function c_group_create_view_from_shape_external_bufferify(self, &
                path, Lpath, type, ndims, shape, external_ptr) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_from_shape_external_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            integer(C_INT), value, intent(IN) :: type
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_LONG), intent(IN) :: shape(*)
            type(C_PTR), value, intent(IN) :: external_ptr
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_from_shape_external_bufferify

        function c_group_create_view_into_buffer(self, path, buff) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_into_buffer")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            type(C_PTR), value, intent(IN) :: buff
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_into_buffer

        function c_group_create_view_into_buffer_bufferify(self, path, &
                Lpath, buff) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_into_buffer_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            type(C_PTR), value, intent(IN) :: buff
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_into_buffer_bufferify

        function c_group_create_view_external(self, path, external_ptr) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_external")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            type(C_PTR), value, intent(IN) :: external_ptr
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_external

        function c_group_create_view_external_bufferify(self, path, &
                Lpath, external_ptr) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_view_external_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            type(C_PTR), value, intent(IN) :: external_ptr
            type(C_PTR) :: SHT_rv
        end function c_group_create_view_external_bufferify

        subroutine c_group_destroy_view(self, path) &
                bind(C, name="SIDRE_group_destroy_view")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
        end subroutine c_group_destroy_view

        subroutine c_group_destroy_view_bufferify(self, path, Lpath) &
                bind(C, name="SIDRE_group_destroy_view_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
        end subroutine c_group_destroy_view_bufferify

        subroutine c_group_destroy_view_and_data_name(self, path) &
                bind(C, name="SIDRE_group_destroy_view_and_data_name")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
        end subroutine c_group_destroy_view_and_data_name

        subroutine c_group_destroy_view_and_data_name_bufferify(self, &
                path, Lpath) &
                bind(C, name="SIDRE_group_destroy_view_and_data_name_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
        end subroutine c_group_destroy_view_and_data_name_bufferify

        subroutine c_group_destroy_view_and_data_index(self, idx) &
                bind(C, name="SIDRE_group_destroy_view_and_data_index")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
        end subroutine c_group_destroy_view_and_data_index

        function c_group_move_view(self, view) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_move_view")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), value, intent(IN) :: view
            type(C_PTR) :: SHT_rv
        end function c_group_move_view

        function c_group_copy_view(self, view) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_copy_view")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), value, intent(IN) :: view
            type(C_PTR) :: SHT_rv
        end function c_group_copy_view

        function c_group_has_group(self, path) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_has_group")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            logical(C_BOOL) :: SHT_rv
        end function c_group_has_group

        function c_group_has_group_bufferify(self, path, Lpath) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_has_group_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            logical(C_BOOL) :: SHT_rv
        end function c_group_has_group_bufferify

        function c_group_has_child_group(self, name) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_has_child_group")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            logical(C_BOOL) :: SHT_rv
        end function c_group_has_child_group

        function c_group_has_child_group_bufferify(self, name, Lname) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_has_child_group_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            logical(C_BOOL) :: SHT_rv
        end function c_group_has_child_group_bufferify

        function c_group_get_group_from_name(self, path) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_group_from_name")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            type(C_PTR) :: SHT_rv
        end function c_group_get_group_from_name

        function c_group_get_group_from_name_bufferify(self, path, &
                Lpath) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_group_from_name_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            type(C_PTR) :: SHT_rv
        end function c_group_get_group_from_name_bufferify

        function c_group_get_group_from_index(self, idx) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_group_from_index")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            type(C_PTR) :: SHT_rv
        end function c_group_get_group_from_index

        pure function c_group_get_group_index(self, name) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_group_index")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT) :: SHT_rv
        end function c_group_get_group_index

        pure function c_group_get_group_index_bufferify(self, name, &
                Lname) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_group_index_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN) :: Lname
            integer(C_INT) :: SHT_rv
        end function c_group_get_group_index_bufferify

        pure function c_group_get_group_name(self, idx) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_group_name")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            type(C_PTR) SHT_rv
        end function c_group_get_group_name

        subroutine c_group_get_group_name_bufferify(self, idx, SHF_rv, &
                NSHF_rv) &
                bind(C, name="SIDRE_group_get_group_name_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            character(kind=C_CHAR), intent(OUT) :: SHF_rv(*)
            integer(C_INT), value, intent(IN) :: NSHF_rv
        end subroutine c_group_get_group_name_bufferify

        pure function c_group_get_first_valid_group_index(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_first_valid_group_index")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: SHT_rv
        end function c_group_get_first_valid_group_index

        pure function c_group_get_next_valid_group_index(self, idx) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_get_next_valid_group_index")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
            integer(C_INT) :: SHT_rv
        end function c_group_get_next_valid_group_index

        function c_group_create_group(self, path) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_group")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            type(C_PTR) :: SHT_rv
        end function c_group_create_group

        function c_group_create_group_bufferify(self, path, Lpath) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_create_group_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
            type(C_PTR) :: SHT_rv
        end function c_group_create_group_bufferify

        subroutine c_group_destroy_group_name(self, path) &
                bind(C, name="SIDRE_group_destroy_group_name")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
        end subroutine c_group_destroy_group_name

        subroutine c_group_destroy_group_name_bufferify(self, path, &
                Lpath) &
                bind(C, name="SIDRE_group_destroy_group_name_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: path(*)
            integer(C_INT), value, intent(IN) :: Lpath
        end subroutine c_group_destroy_group_name_bufferify

        subroutine c_group_destroy_group_index(self, idx) &
                bind(C, name="SIDRE_group_destroy_group_index")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: idx
        end subroutine c_group_destroy_group_index

        function c_group_move_group(self, grp) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_move_group")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), value, intent(IN) :: grp
            type(C_PTR) :: SHT_rv
        end function c_group_move_group

        subroutine c_group_print(self) &
                bind(C, name="SIDRE_group_print")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine c_group_print

        pure function c_group_is_equivalent_to(self, other) &
                result(SHT_rv) &
                bind(C, name="SIDRE_group_is_equivalent_to")
            use iso_c_binding, only : C_BOOL, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), value, intent(IN) :: other
            logical(C_BOOL) :: SHT_rv
        end function c_group_is_equivalent_to

        subroutine c_group_save(self, file_path, protocol) &
                bind(C, name="SIDRE_group_save")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: file_path(*)
            character(kind=C_CHAR), intent(IN) :: protocol(*)
        end subroutine c_group_save

        subroutine c_group_save_bufferify(self, file_path, Lfile_path, &
                protocol, Lprotocol) &
                bind(C, name="SIDRE_group_save_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: file_path(*)
            integer(C_INT), value, intent(IN) :: Lfile_path
            character(kind=C_CHAR), intent(IN) :: protocol(*)
            integer(C_INT), value, intent(IN) :: Lprotocol
        end subroutine c_group_save_bufferify

        subroutine c_group_load_0(self, file_path, protocol) &
                bind(C, name="SIDRE_group_load_0")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: file_path(*)
            character(kind=C_CHAR), intent(IN) :: protocol(*)
        end subroutine c_group_load_0

        subroutine c_group_load_0_bufferify(self, file_path, Lfile_path, &
                protocol, Lprotocol) &
                bind(C, name="SIDRE_group_load_0_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: file_path(*)
            integer(C_INT), value, intent(IN) :: Lfile_path
            character(kind=C_CHAR), intent(IN) :: protocol(*)
            integer(C_INT), value, intent(IN) :: Lprotocol
        end subroutine c_group_load_0_bufferify

        subroutine c_group_load_1(self, file_path, protocol, &
                preserve_contents) &
                bind(C, name="SIDRE_group_load_1")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: file_path(*)
            character(kind=C_CHAR), intent(IN) :: protocol(*)
            logical(C_BOOL), value, intent(IN) :: preserve_contents
        end subroutine c_group_load_1

        subroutine c_group_load_1_bufferify(self, file_path, Lfile_path, &
                protocol, Lprotocol, preserve_contents) &
                bind(C, name="SIDRE_group_load_1_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: file_path(*)
            integer(C_INT), value, intent(IN) :: Lfile_path
            character(kind=C_CHAR), intent(IN) :: protocol(*)
            integer(C_INT), value, intent(IN) :: Lprotocol
            logical(C_BOOL), value, intent(IN) :: preserve_contents
        end subroutine c_group_load_1_bufferify

        subroutine c_group_load_external_data(self, file_path) &
                bind(C, name="SIDRE_group_load_external_data")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: file_path(*)
        end subroutine c_group_load_external_data

        subroutine c_group_load_external_data_bufferify(self, file_path, &
                Lfile_path) &
                bind(C, name="SIDRE_group_load_external_data_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: file_path(*)
            integer(C_INT), value, intent(IN) :: Lfile_path
        end subroutine c_group_load_external_data_bufferify

        ! splicer begin class.Group.additional_interfaces
        ! splicer end class.Group.additional_interfaces

        pure function c_buffer_get_index(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_buffer_get_index")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: SHT_rv
        end function c_buffer_get_index

        pure function c_buffer_get_num_views(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_buffer_get_num_views")
            use iso_c_binding, only : C_PTR, C_SIZE_T
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: SHT_rv
        end function c_buffer_get_num_views

        subroutine c_buffer_describe(self, type, num_elems) &
                bind(C, name="SIDRE_buffer_describe")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
        end subroutine c_buffer_describe

        subroutine c_buffer_allocate_existing(self) &
                bind(C, name="SIDRE_buffer_allocate_existing")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine c_buffer_allocate_existing

        subroutine c_buffer_allocate_from_type(self, type, num_elems) &
                bind(C, name="SIDRE_buffer_allocate_from_type")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
        end subroutine c_buffer_allocate_from_type

        subroutine c_buffer_reallocate(self, num_elems) &
                bind(C, name="SIDRE_buffer_reallocate")
            use iso_c_binding, only : C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_LONG), value, intent(IN) :: num_elems
        end subroutine c_buffer_reallocate

        function c_buffer_get_void_ptr(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_buffer_get_void_ptr")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: SHT_rv
        end function c_buffer_get_void_ptr

        pure function c_buffer_get_type_id(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_buffer_get_type_id")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: SHT_rv
        end function c_buffer_get_type_id

        pure function c_buffer_get_num_elements(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_buffer_get_num_elements")
            use iso_c_binding, only : C_PTR, C_SIZE_T
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: SHT_rv
        end function c_buffer_get_num_elements

        pure function c_buffer_get_total_bytes(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_buffer_get_total_bytes")
            use iso_c_binding, only : C_PTR, C_SIZE_T
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: SHT_rv
        end function c_buffer_get_total_bytes

        pure function c_buffer_get_bytes_per_element(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_buffer_get_bytes_per_element")
            use iso_c_binding, only : C_PTR, C_SIZE_T
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: SHT_rv
        end function c_buffer_get_bytes_per_element

        subroutine c_buffer_print(self) &
                bind(C, name="SIDRE_buffer_print")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine c_buffer_print

        ! splicer begin class.Buffer.additional_interfaces
        ! splicer end class.Buffer.additional_interfaces

        subroutine c_view_allocate_simple(self) &
                bind(C, name="SIDRE_view_allocate_simple")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine c_view_allocate_simple

        subroutine c_view_allocate_from_type(self, type, num_elems) &
                bind(C, name="SIDRE_view_allocate_from_type")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
        end subroutine c_view_allocate_from_type

        subroutine c_view_reallocate(self, num_elems) &
                bind(C, name="SIDRE_view_reallocate")
            use iso_c_binding, only : C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_LONG), value, intent(IN) :: num_elems
        end subroutine c_view_reallocate

        subroutine c_view_apply_0(self) &
                bind(C, name="SIDRE_view_apply_0")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine c_view_apply_0

        subroutine c_view_attach_buffer_only(self, buff) &
                bind(C, name="SIDRE_view_attach_buffer_only")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), value, intent(IN) :: buff
        end subroutine c_view_attach_buffer_only

        subroutine c_view_attach_buffer_type(self, type, num_elems, &
                buff) &
                bind(C, name="SIDRE_view_attach_buffer_type")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            type(C_PTR), value, intent(IN) :: buff
        end subroutine c_view_attach_buffer_type

        subroutine c_view_attach_buffer_shape(self, type, ndims, shape, &
                buff) &
                bind(C, name="SIDRE_view_attach_buffer_shape")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_LONG), intent(IN) :: shape(*)
            type(C_PTR), value, intent(IN) :: buff
        end subroutine c_view_attach_buffer_shape

        subroutine c_view_apply_nelems(self, num_elems) &
                bind(C, name="SIDRE_view_apply_nelems")
            use iso_c_binding, only : C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_LONG), value, intent(IN) :: num_elems
        end subroutine c_view_apply_nelems

        subroutine c_view_apply_nelems_offset(self, num_elems, offset) &
                bind(C, name="SIDRE_view_apply_nelems_offset")
            use iso_c_binding, only : C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_LONG), value, intent(IN) :: num_elems
            integer(C_LONG), value, intent(IN) :: offset
        end subroutine c_view_apply_nelems_offset

        subroutine c_view_apply_nelems_offset_stride(self, num_elems, &
                offset, stride) &
                bind(C, name="SIDRE_view_apply_nelems_offset_stride")
            use iso_c_binding, only : C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_LONG), value, intent(IN) :: num_elems
            integer(C_LONG), value, intent(IN) :: offset
            integer(C_LONG), value, intent(IN) :: stride
        end subroutine c_view_apply_nelems_offset_stride

        subroutine c_view_apply_type_nelems(self, type, num_elems) &
                bind(C, name="SIDRE_view_apply_type_nelems")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
        end subroutine c_view_apply_type_nelems

        subroutine c_view_apply_type_nelems_offset(self, type, &
                num_elems, offset) &
                bind(C, name="SIDRE_view_apply_type_nelems_offset")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            integer(C_LONG), value, intent(IN) :: offset
        end subroutine c_view_apply_type_nelems_offset

        subroutine c_view_apply_type_nelems_offset_stride(self, type, &
                num_elems, offset, stride) &
                bind(C, name="SIDRE_view_apply_type_nelems_offset_stride")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            integer(C_LONG), value, intent(IN) :: offset
            integer(C_LONG), value, intent(IN) :: stride
        end subroutine c_view_apply_type_nelems_offset_stride

        subroutine c_view_apply_type_shape(self, type, ndims, shape) &
                bind(C, name="SIDRE_view_apply_type_shape")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_LONG), intent(IN) :: shape(*)
        end subroutine c_view_apply_type_shape

        pure function c_view_has_buffer(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_has_buffer")
            use iso_c_binding, only : C_BOOL, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: SHT_rv
        end function c_view_has_buffer

        pure function c_view_is_external(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_is_external")
            use iso_c_binding, only : C_BOOL, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: SHT_rv
        end function c_view_is_external

        function c_view_is_allocated(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_is_allocated")
            use iso_c_binding, only : C_BOOL, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: SHT_rv
        end function c_view_is_allocated

        pure function c_view_is_applied(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_is_applied")
            use iso_c_binding, only : C_BOOL, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: SHT_rv
        end function c_view_is_applied

        pure function c_view_is_described(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_is_described")
            use iso_c_binding, only : C_BOOL, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: SHT_rv
        end function c_view_is_described

        pure function c_view_is_empty(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_is_empty")
            use iso_c_binding, only : C_BOOL, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: SHT_rv
        end function c_view_is_empty

        pure function c_view_is_opaque(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_is_opaque")
            use iso_c_binding, only : C_BOOL, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: SHT_rv
        end function c_view_is_opaque

        pure function c_view_is_scalar(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_is_scalar")
            use iso_c_binding, only : C_BOOL, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: SHT_rv
        end function c_view_is_scalar

        pure function c_view_is_string(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_is_string")
            use iso_c_binding, only : C_BOOL, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            logical(C_BOOL) :: SHT_rv
        end function c_view_is_string

        function c_view_get_index(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_index")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: SHT_rv
        end function c_view_get_index

        pure function c_view_get_name(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_name")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) SHT_rv
        end function c_view_get_name

        subroutine c_view_get_name_bufferify(self, SHF_rv, NSHF_rv) &
                bind(C, name="SIDRE_view_get_name_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(OUT) :: SHF_rv(*)
            integer(C_INT), value, intent(IN) :: NSHF_rv
        end subroutine c_view_get_name_bufferify

        subroutine c_view_get_path_bufferify(self, SHF_rv, NSHF_rv) &
                bind(C, name="SIDRE_view_get_path_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(OUT) :: SHF_rv(*)
            integer(C_INT), value, intent(IN) :: NSHF_rv
        end subroutine c_view_get_path_bufferify

        subroutine c_view_get_path_name_bufferify(self, SHF_rv, NSHF_rv) &
                bind(C, name="SIDRE_view_get_path_name_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(OUT) :: SHF_rv(*)
            integer(C_INT), value, intent(IN) :: NSHF_rv
        end subroutine c_view_get_path_name_bufferify

        function c_view_get_buffer(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_buffer")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: SHT_rv
        end function c_view_get_buffer

        pure function c_view_get_void_ptr(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_void_ptr")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: SHT_rv
        end function c_view_get_void_ptr

        subroutine c_view_set_scalar_int(self, value) &
                bind(C, name="SIDRE_view_set_scalar_int")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: value
        end subroutine c_view_set_scalar_int

        subroutine c_view_set_scalar_long(self, value) &
                bind(C, name="SIDRE_view_set_scalar_long")
            use iso_c_binding, only : C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_LONG), value, intent(IN) :: value
        end subroutine c_view_set_scalar_long

        subroutine c_view_set_scalar_float(self, value) &
                bind(C, name="SIDRE_view_set_scalar_float")
            use iso_c_binding, only : C_FLOAT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_FLOAT), value, intent(IN) :: value
        end subroutine c_view_set_scalar_float

        subroutine c_view_set_scalar_double(self, value) &
                bind(C, name="SIDRE_view_set_scalar_double")
            use iso_c_binding, only : C_DOUBLE, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_DOUBLE), value, intent(IN) :: value
        end subroutine c_view_set_scalar_double

        subroutine c_view_set_external_data_ptr_only(self, external_ptr) &
                bind(C, name="SIDRE_view_set_external_data_ptr_only")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR), value, intent(IN) :: external_ptr
        end subroutine c_view_set_external_data_ptr_only

        subroutine c_view_set_external_data_ptr_type(self, type, &
                num_elems, external_ptr) &
                bind(C, name="SIDRE_view_set_external_data_ptr_type")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_LONG), value, intent(IN) :: num_elems
            type(C_PTR), value, intent(IN) :: external_ptr
        end subroutine c_view_set_external_data_ptr_type

        subroutine c_view_set_string(self, value) &
                bind(C, name="SIDRE_view_set_string")
            use iso_c_binding, only : C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: value(*)
        end subroutine c_view_set_string

        subroutine c_view_set_string_bufferify(self, value, Lvalue) &
                bind(C, name="SIDRE_view_set_string_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: value(*)
            integer(C_INT), value, intent(IN) :: Lvalue
        end subroutine c_view_set_string_bufferify

        subroutine c_view_set_external_data_ptr_shape(self, type, ndims, &
                shape, external_ptr) &
                bind(C, name="SIDRE_view_set_external_data_ptr_shape")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: type
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_LONG), intent(IN) :: shape(*)
            type(C_PTR), value, intent(IN) :: external_ptr
        end subroutine c_view_set_external_data_ptr_shape

        function c_view_get_string(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_string")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) SHT_rv
        end function c_view_get_string

        subroutine c_view_get_string_bufferify(self, name, Nname) &
                bind(C, name="SIDRE_view_get_string_bufferify")
            use iso_c_binding, only : C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(OUT) :: name(*)
            integer(C_INT), value, intent(IN) :: Nname
        end subroutine c_view_get_string_bufferify

        function c_view_get_data_int(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_data_int")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: SHT_rv
        end function c_view_get_data_int

        function c_view_get_data_long(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_data_long")
            use iso_c_binding, only : C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_LONG) :: SHT_rv
        end function c_view_get_data_long

        function c_view_get_data_float(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_data_float")
            use iso_c_binding, only : C_FLOAT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_FLOAT) :: SHT_rv
        end function c_view_get_data_float

        function c_view_get_data_double(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_data_double")
            use iso_c_binding, only : C_DOUBLE, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            real(C_DOUBLE) :: SHT_rv
        end function c_view_get_data_double

        function c_view_get_owning_group(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_owning_group")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            type(C_PTR) :: SHT_rv
        end function c_view_get_owning_group

        pure function c_view_get_type_id(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_type_id")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: SHT_rv
        end function c_view_get_type_id

        pure function c_view_get_total_bytes(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_total_bytes")
            use iso_c_binding, only : C_PTR, C_SIZE_T
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: SHT_rv
        end function c_view_get_total_bytes

        pure function c_view_get_bytes_per_element(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_bytes_per_element")
            use iso_c_binding, only : C_PTR, C_SIZE_T
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: SHT_rv
        end function c_view_get_bytes_per_element

        pure function c_view_get_num_elements(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_num_elements")
            use iso_c_binding, only : C_PTR, C_SIZE_T
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: SHT_rv
        end function c_view_get_num_elements

        pure function c_view_get_offset(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_offset")
            use iso_c_binding, only : C_PTR, C_SIZE_T
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: SHT_rv
        end function c_view_get_offset

        pure function c_view_get_stride(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_stride")
            use iso_c_binding, only : C_PTR, C_SIZE_T
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_SIZE_T) :: SHT_rv
        end function c_view_get_stride

        pure function c_view_get_num_dimensions(self) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_num_dimensions")
            use iso_c_binding, only : C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT) :: SHT_rv
        end function c_view_get_num_dimensions

        function c_view_get_shape(self, ndims, shape) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_get_shape")
            use iso_c_binding, only : C_INT, C_LONG, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            integer(C_INT), value, intent(IN) :: ndims
            integer(C_LONG), intent(OUT) :: shape(*)
            integer(C_INT) :: SHT_rv
        end function c_view_get_shape

        function c_view_rename(self, new_name) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_rename")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: new_name(*)
            logical(C_BOOL) :: SHT_rv
        end function c_view_rename

        function c_view_rename_bufferify(self, new_name, Lnew_name) &
                result(SHT_rv) &
                bind(C, name="SIDRE_view_rename_bufferify")
            use iso_c_binding, only : C_BOOL, C_CHAR, C_INT, C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
            character(kind=C_CHAR), intent(IN) :: new_name(*)
            integer(C_INT), value, intent(IN) :: Lnew_name
            logical(C_BOOL) :: SHT_rv
        end function c_view_rename_bufferify

        subroutine c_view_print(self) &
                bind(C, name="SIDRE_view_print")
            use iso_c_binding, only : C_PTR
            implicit none
            type(C_PTR), value, intent(IN) :: self
        end subroutine c_view_print

        ! splicer begin class.View.additional_interfaces
        ! splicer end class.View.additional_interfaces

        function c_name_is_valid(name) &
                result(SHT_rv) &
                bind(C, name="SIDRE_name_is_valid")
            use iso_c_binding, only : C_BOOL, C_CHAR
            implicit none
            character(kind=C_CHAR), intent(IN) :: name(*)
            logical(C_BOOL) :: SHT_rv
        end function c_name_is_valid

        ! splicer begin additional_interfaces
        function SIDRE_create_array_view(group, name, lname, addr, type, rank, extents) &
              result(rv) bind(C,name="SIDRE_create_array_view")
            use iso_c_binding
            import SIDRE_LENGTH
            type(C_PTR), value, intent(IN)     :: group
            character(kind=C_CHAR), intent(IN) :: name(*)
            integer(C_INT), value, intent(IN)  :: lname
            type(C_PTR), value,     intent(IN) :: addr
            integer(C_INT), value, intent(IN)  :: type
            integer(C_INT), value, intent(IN)  :: rank
            integer(SIDRE_LENGTH), intent(IN)  :: extents(*)
            type(C_PTR) rv
        end function SIDRE_create_array_view
        ! splicer end additional_interfaces
    end interface

contains

    function datastore_new() &
            result(SHT_rv)
        type(SidreDataStore) :: SHT_rv
        ! splicer begin class.DataStore.method.new
        SHT_rv%voidptr = c_datastore_new()
        ! splicer end class.DataStore.method.new
    end function datastore_new

    subroutine datastore_delete(obj)
        use iso_c_binding, only : C_NULL_PTR
        class(SidreDataStore) :: obj
        ! splicer begin class.DataStore.method.delete
        call c_datastore_delete(obj%voidptr)
        obj%voidptr = C_NULL_PTR
        ! splicer end class.DataStore.method.delete
    end subroutine datastore_delete

    function datastore_get_root(obj) &
            result(SHT_rv)
        class(SidreDataStore) :: obj
        type(SidreGroup) :: SHT_rv
        ! splicer begin class.DataStore.method.get_root
        SHT_rv%voidptr = c_datastore_get_root(obj%voidptr)
        ! splicer end class.DataStore.method.get_root
    end function datastore_get_root

    function datastore_get_buffer(obj, idx) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreDataStore) :: obj
        integer(C_INT), value, intent(IN) :: idx
        type(SidreBuffer) :: SHT_rv
        ! splicer begin class.DataStore.method.get_buffer
        SHT_rv%voidptr = c_datastore_get_buffer(obj%voidptr, idx)
        ! splicer end class.DataStore.method.get_buffer
    end function datastore_get_buffer

    function datastore_create_buffer_empty(obj) &
            result(SHT_rv)
        class(SidreDataStore) :: obj
        type(SidreBuffer) :: SHT_rv
        ! splicer begin class.DataStore.method.create_buffer_empty
        SHT_rv%voidptr = c_datastore_create_buffer_empty(obj%voidptr)
        ! splicer end class.DataStore.method.create_buffer_empty
    end function datastore_create_buffer_empty

    function datastore_create_buffer_from_type_int(obj, type, num_elems) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreDataStore) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: num_elems
        type(SidreBuffer) :: SHT_rv
        ! splicer begin class.DataStore.method.create_buffer_from_type_int
        SHT_rv%voidptr = c_datastore_create_buffer_from_type(obj%voidptr, &
            type, int(num_elems, C_LONG))
        ! splicer end class.DataStore.method.create_buffer_from_type_int
    end function datastore_create_buffer_from_type_int

    function datastore_create_buffer_from_type_long(obj, type, &
            num_elems) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreDataStore) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        type(SidreBuffer) :: SHT_rv
        ! splicer begin class.DataStore.method.create_buffer_from_type_long
        SHT_rv%voidptr = c_datastore_create_buffer_from_type(obj%voidptr, &
            type, int(num_elems, C_LONG))
        ! splicer end class.DataStore.method.create_buffer_from_type_long
    end function datastore_create_buffer_from_type_long

    subroutine datastore_destroy_buffer(obj, id)
        use iso_c_binding, only : C_INT
        class(SidreDataStore) :: obj
        integer(C_INT), value, intent(IN) :: id
        ! splicer begin class.DataStore.method.destroy_buffer
        call c_datastore_destroy_buffer(obj%voidptr, id)
        ! splicer end class.DataStore.method.destroy_buffer
    end subroutine datastore_destroy_buffer

    function datastore_get_num_buffers(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_SIZE_T
        class(SidreDataStore) :: obj
        integer(C_SIZE_T) :: SHT_rv
        ! splicer begin class.DataStore.method.get_num_buffers
        SHT_rv = c_datastore_get_num_buffers(obj%voidptr)
        ! splicer end class.DataStore.method.get_num_buffers
    end function datastore_get_num_buffers

    subroutine datastore_print(obj)
        class(SidreDataStore) :: obj
        ! splicer begin class.DataStore.method.print
        call c_datastore_print(obj%voidptr)
        ! splicer end class.DataStore.method.print
    end subroutine datastore_print

    function datastore_get_instance(obj) result (voidptr)
        use iso_c_binding, only: C_PTR
        class(SidreDataStore), intent(IN) :: obj
        type(C_PTR) :: voidptr
        voidptr = obj%voidptr
    end function datastore_get_instance

    subroutine datastore_set_instance(obj, voidptr)
        use iso_c_binding, only: C_PTR
        class(SidreDataStore), intent(INOUT) :: obj
        type(C_PTR), intent(IN) :: voidptr
        obj%voidptr = voidptr
    end subroutine datastore_set_instance

    function datastore_associated(obj) result (rv)
        use iso_c_binding, only: c_associated
        class(SidreDataStore), intent(IN) :: obj
        logical rv
        rv = c_associated(obj%voidptr)
    end function datastore_associated

    ! splicer begin class.DataStore.additional_functions
    ! splicer end class.DataStore.additional_functions

    function group_get_index(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        integer(C_INT) :: SHT_rv
        ! splicer begin class.Group.method.get_index
        SHT_rv = c_group_get_index(obj%voidptr)
        ! splicer end class.Group.method.get_index
    end function group_get_index

    function group_get_name(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        class(SidreGroup) :: obj
        character(kind=C_CHAR, len=MAXNAMESIZE) :: SHT_rv
        ! splicer begin class.Group.method.get_name
        call c_group_get_name_bufferify(obj%voidptr, SHT_rv, &
            len(SHT_rv, kind=C_INT))
        ! splicer end class.Group.method.get_name
    end function group_get_name

    function group_get_path(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        class(SidreGroup) :: obj
        character(kind=C_CHAR, len=MAXNAMESIZE) :: SHT_rv
        ! splicer begin class.Group.method.get_path
        call c_group_get_path_bufferify(obj%voidptr, SHT_rv, &
            len(SHT_rv, kind=C_INT))
        ! splicer end class.Group.method.get_path
    end function group_get_path

    function group_get_path_name(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        class(SidreGroup) :: obj
        character(kind=C_CHAR, len=MAXNAMESIZE) :: SHT_rv
        ! splicer begin class.Group.method.get_path_name
        call c_group_get_path_name_bufferify(obj%voidptr, SHT_rv, &
            len(SHT_rv, kind=C_INT))
        ! splicer end class.Group.method.get_path_name
    end function group_get_path_name

    function group_get_parent(obj) &
            result(SHT_rv)
        class(SidreGroup) :: obj
        type(SidreGroup) :: SHT_rv
        ! splicer begin class.Group.method.get_parent
        SHT_rv%voidptr = c_group_get_parent(obj%voidptr)
        ! splicer end class.Group.method.get_parent
    end function group_get_parent

    function group_get_data_store(obj) &
            result(SHT_rv)
        class(SidreGroup) :: obj
        type(SidreDataStore) :: SHT_rv
        ! splicer begin class.Group.method.get_data_store
        SHT_rv%voidptr = c_group_get_data_store(obj%voidptr)
        ! splicer end class.Group.method.get_data_store
    end function group_get_data_store

    function group_get_num_views(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_SIZE_T
        class(SidreGroup) :: obj
        integer(C_SIZE_T) :: SHT_rv
        ! splicer begin class.Group.method.get_num_views
        SHT_rv = c_group_get_num_views(obj%voidptr)
        ! splicer end class.Group.method.get_num_views
    end function group_get_num_views

    function group_get_num_groups(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_SIZE_T
        class(SidreGroup) :: obj
        integer(C_SIZE_T) :: SHT_rv
        ! splicer begin class.Group.method.get_num_groups
        SHT_rv = c_group_get_num_groups(obj%voidptr)
        ! splicer end class.Group.method.get_num_groups
    end function group_get_num_groups

    function group_has_view(obj, path) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL, C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        logical :: SHT_rv
        ! splicer begin class.Group.method.has_view
        SHT_rv = c_group_has_view_bufferify(obj%voidptr, path, &
            len_trim(path, kind=C_INT))
        ! splicer end class.Group.method.has_view
    end function group_has_view

    function group_has_child_view(obj, name) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL, C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: name
        logical :: SHT_rv
        ! splicer begin class.Group.method.has_child_view
        SHT_rv = c_group_has_child_view_bufferify(obj%voidptr, name, &
            len_trim(name, kind=C_INT))
        ! splicer end class.Group.method.has_child_view
    end function group_has_child_view

    function group_rename(obj, new_name) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL, C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: new_name
        logical :: SHT_rv
        ! splicer begin class.Group.method.rename
        SHT_rv = c_group_rename_bufferify(obj%voidptr, new_name, &
            len_trim(new_name, kind=C_INT))
        ! splicer end class.Group.method.rename
    end function group_rename

    function group_get_view_from_name(obj, path) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.get_view_from_name
        SHT_rv%voidptr = c_group_get_view_from_name_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT))
        ! splicer end class.Group.method.get_view_from_name
    end function group_get_view_from_name

    function group_get_view_from_index(obj, idx) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        integer(C_INT), value, intent(IN) :: idx
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.get_view_from_index
        SHT_rv%voidptr = c_group_get_view_from_index(obj%voidptr, idx)
        ! splicer end class.Group.method.get_view_from_index
    end function group_get_view_from_index

    function group_get_view_index(obj, name) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: name
        integer(C_INT) :: SHT_rv
        ! splicer begin class.Group.method.get_view_index
        SHT_rv = c_group_get_view_index_bufferify(obj%voidptr, name, &
            len_trim(name, kind=C_INT))
        ! splicer end class.Group.method.get_view_index
    end function group_get_view_index

    function group_get_view_name(obj, idx) &
            result(SHT_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        class(SidreGroup) :: obj
        integer(C_INT), value, intent(IN) :: idx
        character(kind=C_CHAR, len=MAXNAMESIZE) :: SHT_rv
        ! splicer begin class.Group.method.get_view_name
        call c_group_get_view_name_bufferify(obj%voidptr, idx, SHT_rv, &
            len(SHT_rv, kind=C_INT))
        ! splicer end class.Group.method.get_view_name
    end function group_get_view_name

    function group_get_first_valid_view_index(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        integer(C_INT) :: SHT_rv
        ! splicer begin class.Group.method.get_first_valid_view_index
        SHT_rv = c_group_get_first_valid_view_index(obj%voidptr)
        ! splicer end class.Group.method.get_first_valid_view_index
    end function group_get_first_valid_view_index

    function group_get_next_valid_view_index(obj, idx) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        integer(C_INT), value, intent(IN) :: idx
        integer(C_INT) :: SHT_rv
        ! splicer begin class.Group.method.get_next_valid_view_index
        SHT_rv = c_group_get_next_valid_view_index(obj%voidptr, idx)
        ! splicer end class.Group.method.get_next_valid_view_index
    end function group_get_next_valid_view_index

    function group_create_view_and_allocate_nelems_int(obj, path, type, &
            num_elems) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: num_elems
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_and_allocate_nelems_int
        SHT_rv%voidptr = c_group_create_view_and_allocate_nelems_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), type, &
            int(num_elems, C_LONG))
        ! splicer end class.Group.method.create_view_and_allocate_nelems_int
    end function group_create_view_and_allocate_nelems_int

    function group_create_view_and_allocate_nelems_long(obj, path, type, &
            num_elems) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_and_allocate_nelems_long
        SHT_rv%voidptr = c_group_create_view_and_allocate_nelems_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), type, &
            int(num_elems, C_LONG))
        ! splicer end class.Group.method.create_view_and_allocate_nelems_long
    end function group_create_view_and_allocate_nelems_long

    function group_create_view_and_allocate_shape(obj, path, type, &
            ndims, shape) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: ndims
        integer(C_LONG), intent(IN) :: shape(*)
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_and_allocate_shape
        SHT_rv%voidptr = c_group_create_view_and_allocate_shape_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), type, ndims, shape)
        ! splicer end class.Group.method.create_view_and_allocate_shape
    end function group_create_view_and_allocate_shape

    function group_create_view_scalar_int(obj, path, value) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: value
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_scalar_int
        SHT_rv%voidptr = c_group_create_view_scalar_int_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), value)
        ! splicer end class.Group.method.create_view_scalar_int
    end function group_create_view_scalar_int

    function group_create_view_scalar_long(obj, path, value) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_LONG), value, intent(IN) :: value
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_scalar_long
        SHT_rv%voidptr = c_group_create_view_scalar_long_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), value)
        ! splicer end class.Group.method.create_view_scalar_long
    end function group_create_view_scalar_long

    function group_create_view_scalar_float(obj, path, value) &
            result(SHT_rv)
        use iso_c_binding, only : C_FLOAT, C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        real(C_FLOAT), value, intent(IN) :: value
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_scalar_float
        SHT_rv%voidptr = c_group_create_view_scalar_float_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), value)
        ! splicer end class.Group.method.create_view_scalar_float
    end function group_create_view_scalar_float

    function group_create_view_scalar_double(obj, path, value) &
            result(SHT_rv)
        use iso_c_binding, only : C_DOUBLE, C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        real(C_DOUBLE), value, intent(IN) :: value
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_scalar_double
        SHT_rv%voidptr = c_group_create_view_scalar_double_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), value)
        ! splicer end class.Group.method.create_view_scalar_double
    end function group_create_view_scalar_double

    function group_create_view_string(obj, path, value) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        character(*), intent(IN) :: value
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_string
        SHT_rv%voidptr = c_group_create_view_string_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), value, &
            len_trim(value, kind=C_INT))
        ! splicer end class.Group.method.create_view_string
    end function group_create_view_string

    function group_create_view_empty(obj, path) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_empty
        SHT_rv%voidptr = c_group_create_view_empty_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT))
        ! splicer end class.Group.method.create_view_empty
    end function group_create_view_empty

    function group_create_view_from_type_int(obj, path, type, num_elems) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: num_elems
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_from_type_int
        SHT_rv%voidptr = c_group_create_view_from_type_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), type, &
            int(num_elems, C_LONG))
        ! splicer end class.Group.method.create_view_from_type_int
    end function group_create_view_from_type_int

    function group_create_view_from_type_long(obj, path, type, &
            num_elems) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_from_type_long
        SHT_rv%voidptr = c_group_create_view_from_type_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), type, &
            int(num_elems, C_LONG))
        ! splicer end class.Group.method.create_view_from_type_long
    end function group_create_view_from_type_long

    function group_create_view_from_type_and_buffer_int(obj, path, type, &
            num_elems, buff) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: num_elems
        type(SidreBuffer), value, intent(IN) :: buff
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_from_type_and_buffer_int
        SHT_rv%voidptr = c_group_create_view_from_type_and_buffer_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), type, &
            int(num_elems, C_LONG), buff%get_instance())
        ! splicer end class.Group.method.create_view_from_type_and_buffer_int
    end function group_create_view_from_type_and_buffer_int

    function group_create_view_from_type_and_buffer_long(obj, path, &
            type, num_elems, buff) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        type(SidreBuffer), value, intent(IN) :: buff
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_from_type_and_buffer_long
        SHT_rv%voidptr = c_group_create_view_from_type_and_buffer_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), type, &
            int(num_elems, C_LONG), buff%get_instance())
        ! splicer end class.Group.method.create_view_from_type_and_buffer_long
    end function group_create_view_from_type_and_buffer_long

    function group_create_view_from_type_external_int(obj, path, type, &
            num_elems, external_ptr) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG, C_PTR
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: num_elems
        type(C_PTR), value, intent(IN) :: external_ptr
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_from_type_external_int
        SHT_rv%voidptr = c_group_create_view_from_type_external_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), type, &
            int(num_elems, C_LONG), external_ptr)
        ! splicer end class.Group.method.create_view_from_type_external_int
    end function group_create_view_from_type_external_int

    function group_create_view_from_type_external_long(obj, path, type, &
            num_elems, external_ptr) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG, C_PTR
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        type(C_PTR), value, intent(IN) :: external_ptr
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_from_type_external_long
        SHT_rv%voidptr = c_group_create_view_from_type_external_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), type, &
            int(num_elems, C_LONG), external_ptr)
        ! splicer end class.Group.method.create_view_from_type_external_long
    end function group_create_view_from_type_external_long

    function group_create_view_from_shape(obj, path, type, ndims, shape) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: ndims
        integer(C_LONG), intent(IN) :: shape(*)
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_from_shape
        SHT_rv%voidptr = c_group_create_view_from_shape_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), type, ndims, shape)
        ! splicer end class.Group.method.create_view_from_shape
    end function group_create_view_from_shape

    function group_create_view_from_shape_and_buffer(obj, path, type, &
            ndims, shape, buff) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: ndims
        integer(C_LONG), intent(IN) :: shape(*)
        type(SidreBuffer), value, intent(IN) :: buff
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_from_shape_and_buffer
        SHT_rv%voidptr = c_group_create_view_from_shape_and_buffer_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), type, ndims, shape, &
            buff%get_instance())
        ! splicer end class.Group.method.create_view_from_shape_and_buffer
    end function group_create_view_from_shape_and_buffer

    function group_create_view_from_shape_external(obj, path, type, &
            ndims, shape, external_ptr) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG, C_PTR
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: ndims
        integer(C_LONG), intent(IN) :: shape(*)
        type(C_PTR), value, intent(IN) :: external_ptr
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_from_shape_external
        SHT_rv%voidptr = c_group_create_view_from_shape_external_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), type, ndims, shape, &
            external_ptr)
        ! splicer end class.Group.method.create_view_from_shape_external
    end function group_create_view_from_shape_external

    function group_create_view_into_buffer(obj, path, buff) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        type(SidreBuffer), value, intent(IN) :: buff
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_into_buffer
        SHT_rv%voidptr = c_group_create_view_into_buffer_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), buff%get_instance())
        ! splicer end class.Group.method.create_view_into_buffer
    end function group_create_view_into_buffer

    function group_create_view_external(obj, path, external_ptr) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_PTR
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        type(C_PTR), value, intent(IN) :: external_ptr
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.create_view_external
        SHT_rv%voidptr = c_group_create_view_external_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT), external_ptr)
        ! splicer end class.Group.method.create_view_external
    end function group_create_view_external

    subroutine group_destroy_view(obj, path)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        ! splicer begin class.Group.method.destroy_view
        call c_group_destroy_view_bufferify(obj%voidptr, path, &
            len_trim(path, kind=C_INT))
        ! splicer end class.Group.method.destroy_view
    end subroutine group_destroy_view

    subroutine group_destroy_view_and_data_name(obj, path)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        ! splicer begin class.Group.method.destroy_view_and_data_name
        call c_group_destroy_view_and_data_name_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT))
        ! splicer end class.Group.method.destroy_view_and_data_name
    end subroutine group_destroy_view_and_data_name

    subroutine group_destroy_view_and_data_index(obj, idx)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        integer(C_INT), value, intent(IN) :: idx
        ! splicer begin class.Group.method.destroy_view_and_data_index
        call c_group_destroy_view_and_data_index(obj%voidptr, idx)
        ! splicer end class.Group.method.destroy_view_and_data_index
    end subroutine group_destroy_view_and_data_index

    function group_move_view(obj, view) &
            result(SHT_rv)
        class(SidreGroup) :: obj
        type(SidreView), value, intent(IN) :: view
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.move_view
        SHT_rv%voidptr = c_group_move_view(obj%voidptr, &
            view%get_instance())
        ! splicer end class.Group.method.move_view
    end function group_move_view

    function group_copy_view(obj, view) &
            result(SHT_rv)
        class(SidreGroup) :: obj
        type(SidreView), value, intent(IN) :: view
        type(SidreView) :: SHT_rv
        ! splicer begin class.Group.method.copy_view
        SHT_rv%voidptr = c_group_copy_view(obj%voidptr, &
            view%get_instance())
        ! splicer end class.Group.method.copy_view
    end function group_copy_view

    function group_has_group(obj, path) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL, C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        logical :: SHT_rv
        ! splicer begin class.Group.method.has_group
        SHT_rv = c_group_has_group_bufferify(obj%voidptr, path, &
            len_trim(path, kind=C_INT))
        ! splicer end class.Group.method.has_group
    end function group_has_group

    function group_has_child_group(obj, name) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL, C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: name
        logical :: SHT_rv
        ! splicer begin class.Group.method.has_child_group
        SHT_rv = c_group_has_child_group_bufferify(obj%voidptr, name, &
            len_trim(name, kind=C_INT))
        ! splicer end class.Group.method.has_child_group
    end function group_has_child_group

    function group_get_group_from_name(obj, path) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        type(SidreGroup) :: SHT_rv
        ! splicer begin class.Group.method.get_group_from_name
        SHT_rv%voidptr = c_group_get_group_from_name_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT))
        ! splicer end class.Group.method.get_group_from_name
    end function group_get_group_from_name

    function group_get_group_from_index(obj, idx) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        integer(C_INT), value, intent(IN) :: idx
        type(SidreGroup) :: SHT_rv
        ! splicer begin class.Group.method.get_group_from_index
        SHT_rv%voidptr = c_group_get_group_from_index(obj%voidptr, idx)
        ! splicer end class.Group.method.get_group_from_index
    end function group_get_group_from_index

    function group_get_group_index(obj, name) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: name
        integer(C_INT) :: SHT_rv
        ! splicer begin class.Group.method.get_group_index
        SHT_rv = c_group_get_group_index_bufferify(obj%voidptr, name, &
            len_trim(name, kind=C_INT))
        ! splicer end class.Group.method.get_group_index
    end function group_get_group_index

    function group_get_group_name(obj, idx) &
            result(SHT_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        class(SidreGroup) :: obj
        integer(C_INT), value, intent(IN) :: idx
        character(kind=C_CHAR, len=MAXNAMESIZE) :: SHT_rv
        ! splicer begin class.Group.method.get_group_name
        call c_group_get_group_name_bufferify(obj%voidptr, idx, SHT_rv, &
            len(SHT_rv, kind=C_INT))
        ! splicer end class.Group.method.get_group_name
    end function group_get_group_name

    function group_get_first_valid_group_index(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        integer(C_INT) :: SHT_rv
        ! splicer begin class.Group.method.get_first_valid_group_index
        SHT_rv = c_group_get_first_valid_group_index(obj%voidptr)
        ! splicer end class.Group.method.get_first_valid_group_index
    end function group_get_first_valid_group_index

    function group_get_next_valid_group_index(obj, idx) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        integer(C_INT), value, intent(IN) :: idx
        integer(C_INT) :: SHT_rv
        ! splicer begin class.Group.method.get_next_valid_group_index
        SHT_rv = c_group_get_next_valid_group_index(obj%voidptr, idx)
        ! splicer end class.Group.method.get_next_valid_group_index
    end function group_get_next_valid_group_index

    function group_create_group(obj, path) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        type(SidreGroup) :: SHT_rv
        ! splicer begin class.Group.method.create_group
        SHT_rv%voidptr = c_group_create_group_bufferify(obj%voidptr, &
            path, len_trim(path, kind=C_INT))
        ! splicer end class.Group.method.create_group
    end function group_create_group

    subroutine group_destroy_group_name(obj, path)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: path
        ! splicer begin class.Group.method.destroy_group_name
        call c_group_destroy_group_name_bufferify(obj%voidptr, path, &
            len_trim(path, kind=C_INT))
        ! splicer end class.Group.method.destroy_group_name
    end subroutine group_destroy_group_name

    subroutine group_destroy_group_index(obj, idx)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        integer(C_INT), value, intent(IN) :: idx
        ! splicer begin class.Group.method.destroy_group_index
        call c_group_destroy_group_index(obj%voidptr, idx)
        ! splicer end class.Group.method.destroy_group_index
    end subroutine group_destroy_group_index

    function group_move_group(obj, grp) &
            result(SHT_rv)
        class(SidreGroup) :: obj
        type(SidreGroup), value, intent(IN) :: grp
        type(SidreGroup) :: SHT_rv
        ! splicer begin class.Group.method.move_group
        SHT_rv%voidptr = c_group_move_group(obj%voidptr, &
            grp%get_instance())
        ! splicer end class.Group.method.move_group
    end function group_move_group

    subroutine group_print(obj)
        class(SidreGroup) :: obj
        ! splicer begin class.Group.method.print
        call c_group_print(obj%voidptr)
        ! splicer end class.Group.method.print
    end subroutine group_print

    function group_is_equivalent_to(obj, other) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        class(SidreGroup) :: obj
        type(SidreGroup), value, intent(IN) :: other
        logical :: SHT_rv
        ! splicer begin class.Group.method.is_equivalent_to
        SHT_rv = c_group_is_equivalent_to(obj%voidptr, &
            other%get_instance())
        ! splicer end class.Group.method.is_equivalent_to
    end function group_is_equivalent_to

    subroutine group_save(obj, file_path, protocol)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: file_path
        character(*), intent(IN) :: protocol
        ! splicer begin class.Group.method.save
        call c_group_save_bufferify(obj%voidptr, file_path, &
            len_trim(file_path, kind=C_INT), protocol, &
            len_trim(protocol, kind=C_INT))
        ! splicer end class.Group.method.save
    end subroutine group_save

    subroutine group_load_0(obj, file_path, protocol)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: file_path
        character(*), intent(IN) :: protocol
        ! splicer begin class.Group.method.load_0
        call c_group_load_0_bufferify(obj%voidptr, file_path, &
            len_trim(file_path, kind=C_INT), protocol, &
            len_trim(protocol, kind=C_INT))
        ! splicer end class.Group.method.load_0
    end subroutine group_load_0

    subroutine group_load_1(obj, file_path, protocol, preserve_contents)
        use iso_c_binding, only : C_BOOL, C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: file_path
        character(*), intent(IN) :: protocol
        logical, value, intent(IN) :: preserve_contents
        logical(C_BOOL) SH_preserve_contents
        SH_preserve_contents = preserve_contents  ! coerce to C_BOOL
        ! splicer begin class.Group.method.load_1
        call c_group_load_1_bufferify(obj%voidptr, file_path, &
            len_trim(file_path, kind=C_INT), protocol, &
            len_trim(protocol, kind=C_INT), SH_preserve_contents)
        ! splicer end class.Group.method.load_1
    end subroutine group_load_1

    subroutine group_load_external_data(obj, file_path)
        use iso_c_binding, only : C_INT
        class(SidreGroup) :: obj
        character(*), intent(IN) :: file_path
        ! splicer begin class.Group.method.load_external_data
        call c_group_load_external_data_bufferify(obj%voidptr, &
            file_path, len_trim(file_path, kind=C_INT))
        ! splicer end class.Group.method.load_external_data
    end subroutine group_load_external_data

    function group_get_instance(obj) result (voidptr)
        use iso_c_binding, only: C_PTR
        class(SidreGroup), intent(IN) :: obj
        type(C_PTR) :: voidptr
        voidptr = obj%voidptr
    end function group_get_instance

    subroutine group_set_instance(obj, voidptr)
        use iso_c_binding, only: C_PTR
        class(SidreGroup), intent(INOUT) :: obj
        type(C_PTR), intent(IN) :: voidptr
        obj%voidptr = voidptr
    end subroutine group_set_instance

    function group_associated(obj) result (rv)
        use iso_c_binding, only: c_associated
        class(SidreGroup), intent(IN) :: obj
        logical rv
        rv = c_associated(obj%voidptr)
    end function group_associated

    ! splicer begin class.Group.additional_functions

    ! Generated by genfsidresplicer.py
    subroutine group_get_scalar_int(grp, name, value)
        use iso_c_binding
        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_INT), intent(OUT) :: value
        integer(C_INT) :: lname
        type(C_PTR) view

        lname = len_trim(name)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        value = c_view_get_data_int(view)
    end subroutine group_get_scalar_int

    ! Generated by genfsidresplicer.py
    subroutine group_get_scalar_long(grp, name, value)
        use iso_c_binding
        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_LONG), intent(OUT) :: value
        integer(C_INT) :: lname
        type(C_PTR) view

        lname = len_trim(name)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        value = c_view_get_data_long(view)
    end subroutine group_get_scalar_long

    ! Generated by genfsidresplicer.py
    subroutine group_get_scalar_float(grp, name, value)
        use iso_c_binding
        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_FLOAT), intent(OUT) :: value
        integer(C_INT) :: lname
        type(C_PTR) view

        lname = len_trim(name)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        value = c_view_get_data_float(view)
    end subroutine group_get_scalar_float

    ! Generated by genfsidresplicer.py
    subroutine group_get_scalar_double(grp, name, value)
        use iso_c_binding
        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_DOUBLE), intent(OUT) :: value
        integer(C_INT) :: lname
        type(C_PTR) view

        lname = len_trim(name)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        value = c_view_get_data_double(view)
    end subroutine group_get_scalar_double

    ! Generated by genfsidresplicer.py
    subroutine group_set_scalar_int(grp, name, value)
        use iso_c_binding
        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_INT), intent(IN) :: value
        integer(C_INT) :: lname
        type(C_PTR) view

        lname = len_trim(name)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        call c_view_set_scalar_int(view, value)
    end subroutine group_set_scalar_int

    ! Generated by genfsidresplicer.py
    subroutine group_set_scalar_long(grp, name, value)
        use iso_c_binding
        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_LONG), intent(IN) :: value
        integer(C_INT) :: lname
        type(C_PTR) view

        lname = len_trim(name)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        call c_view_set_scalar_long(view, value)
    end subroutine group_set_scalar_long

    ! Generated by genfsidresplicer.py
    subroutine group_set_scalar_float(grp, name, value)
        use iso_c_binding
        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_FLOAT), intent(IN) :: value
        integer(C_INT) :: lname
        type(C_PTR) view

        lname = len_trim(name)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        call c_view_set_scalar_float(view, value)
    end subroutine group_set_scalar_float

    ! Generated by genfsidresplicer.py
    subroutine group_set_scalar_double(grp, name, value)
        use iso_c_binding
        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_DOUBLE), intent(IN) :: value
        integer(C_INT) :: lname
        type(C_PTR) view

        lname = len_trim(name)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        call c_view_set_scalar_double(view, value)
    end subroutine group_set_scalar_double

    ! Generated by genfsidresplicer.py
    function group_create_array_view_int_scalar(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_INT), target, intent(IN) :: value
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(1)
        integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value, addr)
#endif
        extents(1) = 1_SIDRE_LENGTH
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 0, extents)
    end function group_create_array_view_int_scalar

    ! Generated by genfsidresplicer.py
    function group_create_array_view_int_1d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_INT), target, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(1)
        integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 1, extents)
    end function group_create_array_view_int_1d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_int_2d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_INT), target, intent(IN) :: value(:,:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(2)
        integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 2, extents)
    end function group_create_array_view_int_2d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_int_3d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_INT), target, intent(IN) :: value(:,:,:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(3)
        integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 3, extents)
    end function group_create_array_view_int_3d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_int_4d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_INT), target, intent(IN) :: value(:,:,:,:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(4)
        integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3),lbound(value,4)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 4, extents)
    end function group_create_array_view_int_4d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_long_scalar(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_LONG), target, intent(IN) :: value
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(1)
        integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value, addr)
#endif
        extents(1) = 1_SIDRE_LENGTH
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 0, extents)
    end function group_create_array_view_long_scalar

    ! Generated by genfsidresplicer.py
    function group_create_array_view_long_1d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_LONG), target, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(1)
        integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 1, extents)
    end function group_create_array_view_long_1d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_long_2d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_LONG), target, intent(IN) :: value(:,:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(2)
        integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 2, extents)
    end function group_create_array_view_long_2d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_long_3d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_LONG), target, intent(IN) :: value(:,:,:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(3)
        integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 3, extents)
    end function group_create_array_view_long_3d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_long_4d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        integer(C_LONG), target, intent(IN) :: value(:,:,:,:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(4)
        integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3),lbound(value,4)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 4, extents)
    end function group_create_array_view_long_4d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_float_scalar(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_FLOAT), target, intent(IN) :: value
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(1)
        integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value, addr)
#endif
        extents(1) = 1_SIDRE_LENGTH
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 0, extents)
    end function group_create_array_view_float_scalar

    ! Generated by genfsidresplicer.py
    function group_create_array_view_float_1d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_FLOAT), target, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(1)
        integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 1, extents)
    end function group_create_array_view_float_1d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_float_2d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_FLOAT), target, intent(IN) :: value(:,:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(2)
        integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 2, extents)
    end function group_create_array_view_float_2d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_float_3d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_FLOAT), target, intent(IN) :: value(:,:,:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(3)
        integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 3, extents)
    end function group_create_array_view_float_3d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_float_4d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_FLOAT), target, intent(IN) :: value(:,:,:,:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(4)
        integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3),lbound(value,4)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 4, extents)
    end function group_create_array_view_float_4d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_double_scalar(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_DOUBLE), target, intent(IN) :: value
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(1)
        integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value, addr)
#endif
        extents(1) = 1_SIDRE_LENGTH
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 0, extents)
    end function group_create_array_view_double_scalar

    ! Generated by genfsidresplicer.py
    function group_create_array_view_double_1d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_DOUBLE), target, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(1)
        integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 1, extents)
    end function group_create_array_view_double_1d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_double_2d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_DOUBLE), target, intent(IN) :: value(:,:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(2)
        integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 2, extents)
    end function group_create_array_view_double_2d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_double_3d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_DOUBLE), target, intent(IN) :: value(:,:,:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(3)
        integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 3, extents)
    end function group_create_array_view_double_3d

    ! Generated by genfsidresplicer.py
    function group_create_array_view_double_4d(grp, name, value) result(rv)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        real(C_DOUBLE), target, intent(IN) :: value(:,:,:,:)
        integer(C_INT) :: lname
        type(SidreView) :: rv
        integer(SIDRE_LENGTH) :: extents(4)
        integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

        lname = len_trim(name)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3),lbound(value,4)), addr)
#endif
        extents = shape(value, kind=SIDRE_LENGTH)
        rv%voidptr = c_group_create_view_external_bufferify( &
            grp%voidptr, name, lname, addr)
        call c_view_apply_type_shape(rv%voidptr, type, 4, extents)
    end function group_create_array_view_double_4d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_int_scalar(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        integer(C_INT), target, intent(IN) :: value
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents(1) = 1_SIDRE_LENGTH
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value, addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 0, extents)
        endif
    end subroutine group_set_array_data_ptr_int_scalar

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_int_1d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        integer(C_INT), target, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 1, extents)
        endif
    end subroutine group_set_array_data_ptr_int_1d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_int_2d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        integer(C_INT), target, intent(IN) :: value(:,:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(2)
    !    integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 2, extents)
        endif
    end subroutine group_set_array_data_ptr_int_2d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_int_3d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        integer(C_INT), target, intent(IN) :: value(:,:,:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(3)
    !    integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 3, extents)
        endif
    end subroutine group_set_array_data_ptr_int_3d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_int_4d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        integer(C_INT), target, intent(IN) :: value(:,:,:,:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(4)
    !    integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3),lbound(value,4)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 4, extents)
        endif
    end subroutine group_set_array_data_ptr_int_4d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_long_scalar(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        integer(C_LONG), target, intent(IN) :: value
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents(1) = 1_SIDRE_LENGTH
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value, addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 0, extents)
        endif
    end subroutine group_set_array_data_ptr_long_scalar

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_long_1d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        integer(C_LONG), target, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 1, extents)
        endif
    end subroutine group_set_array_data_ptr_long_1d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_long_2d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        integer(C_LONG), target, intent(IN) :: value(:,:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(2)
    !    integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 2, extents)
        endif
    end subroutine group_set_array_data_ptr_long_2d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_long_3d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        integer(C_LONG), target, intent(IN) :: value(:,:,:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(3)
    !    integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 3, extents)
        endif
    end subroutine group_set_array_data_ptr_long_3d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_long_4d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        integer(C_LONG), target, intent(IN) :: value(:,:,:,:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(4)
    !    integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3),lbound(value,4)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 4, extents)
        endif
    end subroutine group_set_array_data_ptr_long_4d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_float_scalar(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        real(C_FLOAT), target, intent(IN) :: value
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents(1) = 1_SIDRE_LENGTH
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value, addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 0, extents)
        endif
    end subroutine group_set_array_data_ptr_float_scalar

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_float_1d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        real(C_FLOAT), target, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 1, extents)
        endif
    end subroutine group_set_array_data_ptr_float_1d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_float_2d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        real(C_FLOAT), target, intent(IN) :: value(:,:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(2)
    !    integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 2, extents)
        endif
    end subroutine group_set_array_data_ptr_float_2d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_float_3d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        real(C_FLOAT), target, intent(IN) :: value(:,:,:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(3)
    !    integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 3, extents)
        endif
    end subroutine group_set_array_data_ptr_float_3d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_float_4d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        real(C_FLOAT), target, intent(IN) :: value(:,:,:,:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(4)
    !    integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3),lbound(value,4)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 4, extents)
        endif
    end subroutine group_set_array_data_ptr_float_4d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_double_scalar(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        real(C_DOUBLE), target, intent(IN) :: value
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents(1) = 1_SIDRE_LENGTH
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value, addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 0, extents)
        endif
    end subroutine group_set_array_data_ptr_double_scalar

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_double_1d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        real(C_DOUBLE), target, intent(IN) :: value(:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 1, extents)
        endif
    end subroutine group_set_array_data_ptr_double_1d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_double_2d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        real(C_DOUBLE), target, intent(IN) :: value(:,:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(2)
    !    integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 2, extents)
        endif
    end subroutine group_set_array_data_ptr_double_2d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_double_3d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        real(C_DOUBLE), target, intent(IN) :: value(:,:,:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(3)
    !    integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 3, extents)
        endif
    end subroutine group_set_array_data_ptr_double_3d

    ! Generated by genfsidresplicer.py
    ! This function does nothing if view name does not exist in group.
    subroutine group_set_array_data_ptr_double_4d(grp, name, value)
        use iso_c_binding
        implicit none

        class(SidreGroup), intent(IN) :: grp
        character(len=*), intent(IN) :: name
        real(C_DOUBLE), target, intent(IN) :: value(:,:,:,:)
        integer(C_INT) :: lname
        type(C_ptr) view
    !    integer(SIDRE_LENGTH) :: extents(4)
    !    integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

        lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        if (c_associated(view)) then
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
            addr = c_loc(value)
#else
            call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3),lbound(value,4)), addr)
#endif
            call c_view_set_external_data_ptr_only(view, addr)
    !        call c_view_apply_type_shape(rv%voidptr, type, 4, extents)
        endif
    end subroutine group_set_array_data_ptr_double_4d

    subroutine group_get_string(grp, name, value)
        use iso_c_binding
        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        character(*), intent(OUT) :: value
        integer(C_INT) :: lname
        type(C_PTR) view

        lname = len_trim(name)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        call c_view_get_string_bufferify(view, value, len(value, kind=C_INT))
    end subroutine group_get_string

    subroutine group_set_string(grp, name, value)
        use iso_c_binding
        class(SidreGroup), intent(IN) :: grp
        character(*), intent(IN) :: name
        character(*), intent(IN) :: value
        integer(C_INT) :: lname
        type(C_PTR) view

        lname = len_trim(name)
        view = c_group_get_view_from_name_bufferify(grp%voidptr, name, lname)
        call c_view_set_string_bufferify(view, value, len_trim(value, kind=C_INT))
    end subroutine group_set_string

    ! splicer end class.Group.additional_functions

    function buffer_get_index(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreBuffer) :: obj
        integer(C_INT) :: SHT_rv
        ! splicer begin class.Buffer.method.get_index
        SHT_rv = c_buffer_get_index(obj%voidptr)
        ! splicer end class.Buffer.method.get_index
    end function buffer_get_index

    function buffer_get_num_views(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_SIZE_T
        class(SidreBuffer) :: obj
        integer(C_SIZE_T) :: SHT_rv
        ! splicer begin class.Buffer.method.get_num_views
        SHT_rv = c_buffer_get_num_views(obj%voidptr)
        ! splicer end class.Buffer.method.get_num_views
    end function buffer_get_num_views

    subroutine buffer_describe_int(obj, type, num_elems)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreBuffer) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: num_elems
        ! splicer begin class.Buffer.method.describe_int
        call c_buffer_describe(obj%voidptr, type, &
            int(num_elems, C_LONG))
        ! splicer end class.Buffer.method.describe_int
    end subroutine buffer_describe_int

    subroutine buffer_describe_long(obj, type, num_elems)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreBuffer) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        ! splicer begin class.Buffer.method.describe_long
        call c_buffer_describe(obj%voidptr, type, &
            int(num_elems, C_LONG))
        ! splicer end class.Buffer.method.describe_long
    end subroutine buffer_describe_long

    subroutine buffer_allocate_existing(obj)
        class(SidreBuffer) :: obj
        ! splicer begin class.Buffer.method.allocate_existing
        call c_buffer_allocate_existing(obj%voidptr)
        ! splicer end class.Buffer.method.allocate_existing
    end subroutine buffer_allocate_existing

    subroutine buffer_allocate_from_type_int(obj, type, num_elems)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreBuffer) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: num_elems
        ! splicer begin class.Buffer.method.allocate_from_type_int
        call c_buffer_allocate_from_type(obj%voidptr, type, &
            int(num_elems, C_LONG))
        ! splicer end class.Buffer.method.allocate_from_type_int
    end subroutine buffer_allocate_from_type_int

    subroutine buffer_allocate_from_type_long(obj, type, num_elems)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreBuffer) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        ! splicer begin class.Buffer.method.allocate_from_type_long
        call c_buffer_allocate_from_type(obj%voidptr, type, &
            int(num_elems, C_LONG))
        ! splicer end class.Buffer.method.allocate_from_type_long
    end subroutine buffer_allocate_from_type_long

    subroutine buffer_reallocate_int(obj, num_elems)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreBuffer) :: obj
        integer(C_INT), value, intent(IN) :: num_elems
        ! splicer begin class.Buffer.method.reallocate_int
        call c_buffer_reallocate(obj%voidptr, int(num_elems, C_LONG))
        ! splicer end class.Buffer.method.reallocate_int
    end subroutine buffer_reallocate_int

    subroutine buffer_reallocate_long(obj, num_elems)
        use iso_c_binding, only : C_LONG
        class(SidreBuffer) :: obj
        integer(C_LONG), value, intent(IN) :: num_elems
        ! splicer begin class.Buffer.method.reallocate_long
        call c_buffer_reallocate(obj%voidptr, int(num_elems, C_LONG))
        ! splicer end class.Buffer.method.reallocate_long
    end subroutine buffer_reallocate_long

    function buffer_get_void_ptr(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_PTR
        class(SidreBuffer) :: obj
        type(C_PTR) :: SHT_rv
        ! splicer begin class.Buffer.method.get_void_ptr
        SHT_rv = c_buffer_get_void_ptr(obj%voidptr)
        ! splicer end class.Buffer.method.get_void_ptr
    end function buffer_get_void_ptr

    function buffer_get_type_id(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreBuffer) :: obj
        integer(C_INT) :: SHT_rv
        ! splicer begin class.Buffer.method.get_type_id
        SHT_rv = c_buffer_get_type_id(obj%voidptr)
        ! splicer end class.Buffer.method.get_type_id
    end function buffer_get_type_id

    function buffer_get_num_elements(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_SIZE_T
        class(SidreBuffer) :: obj
        integer(C_SIZE_T) :: SHT_rv
        ! splicer begin class.Buffer.method.get_num_elements
        SHT_rv = c_buffer_get_num_elements(obj%voidptr)
        ! splicer end class.Buffer.method.get_num_elements
    end function buffer_get_num_elements

    function buffer_get_total_bytes(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_SIZE_T
        class(SidreBuffer) :: obj
        integer(C_SIZE_T) :: SHT_rv
        ! splicer begin class.Buffer.method.get_total_bytes
        SHT_rv = c_buffer_get_total_bytes(obj%voidptr)
        ! splicer end class.Buffer.method.get_total_bytes
    end function buffer_get_total_bytes

    function buffer_get_bytes_per_element(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_SIZE_T
        class(SidreBuffer) :: obj
        integer(C_SIZE_T) :: SHT_rv
        ! splicer begin class.Buffer.method.get_bytes_per_element
        SHT_rv = c_buffer_get_bytes_per_element(obj%voidptr)
        ! splicer end class.Buffer.method.get_bytes_per_element
    end function buffer_get_bytes_per_element

    subroutine buffer_print(obj)
        class(SidreBuffer) :: obj
        ! splicer begin class.Buffer.method.print
        call c_buffer_print(obj%voidptr)
        ! splicer end class.Buffer.method.print
    end subroutine buffer_print

    function buffer_get_instance(obj) result (voidptr)
        use iso_c_binding, only: C_PTR
        class(SidreBuffer), intent(IN) :: obj
        type(C_PTR) :: voidptr
        voidptr = obj%voidptr
    end function buffer_get_instance

    subroutine buffer_set_instance(obj, voidptr)
        use iso_c_binding, only: C_PTR
        class(SidreBuffer), intent(INOUT) :: obj
        type(C_PTR), intent(IN) :: voidptr
        obj%voidptr = voidptr
    end subroutine buffer_set_instance

    function buffer_associated(obj) result (rv)
        use iso_c_binding, only: c_associated
        class(SidreBuffer), intent(IN) :: obj
        logical rv
        rv = c_associated(obj%voidptr)
    end function buffer_associated

    ! splicer begin class.Buffer.additional_functions
    ! splicer end class.Buffer.additional_functions

    subroutine view_allocate_simple(obj)
        class(SidreView) :: obj
        ! splicer begin class.View.method.allocate_simple
        call c_view_allocate_simple(obj%voidptr)
        ! splicer end class.View.method.allocate_simple
    end subroutine view_allocate_simple

    subroutine view_allocate_from_type_int(obj, type, num_elems)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: num_elems
        ! splicer begin class.View.method.allocate_from_type_int
        call c_view_allocate_from_type(obj%voidptr, type, &
            int(num_elems, C_LONG))
        ! splicer end class.View.method.allocate_from_type_int
    end subroutine view_allocate_from_type_int

    subroutine view_allocate_from_type_long(obj, type, num_elems)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        ! splicer begin class.View.method.allocate_from_type_long
        call c_view_allocate_from_type(obj%voidptr, type, &
            int(num_elems, C_LONG))
        ! splicer end class.View.method.allocate_from_type_long
    end subroutine view_allocate_from_type_long

    subroutine view_reallocate_int(obj, num_elems)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: num_elems
        ! splicer begin class.View.method.reallocate_int
        call c_view_reallocate(obj%voidptr, int(num_elems, C_LONG))
        ! splicer end class.View.method.reallocate_int
    end subroutine view_reallocate_int

    subroutine view_reallocate_long(obj, num_elems)
        use iso_c_binding, only : C_LONG
        class(SidreView) :: obj
        integer(C_LONG), value, intent(IN) :: num_elems
        ! splicer begin class.View.method.reallocate_long
        call c_view_reallocate(obj%voidptr, int(num_elems, C_LONG))
        ! splicer end class.View.method.reallocate_long
    end subroutine view_reallocate_long

    subroutine view_apply_0(obj)
        class(SidreView) :: obj
        ! splicer begin class.View.method.apply_0
        call c_view_apply_0(obj%voidptr)
        ! splicer end class.View.method.apply_0
    end subroutine view_apply_0

    subroutine view_attach_buffer_only(obj, buff)
        class(SidreView) :: obj
        type(SidreBuffer), value, intent(IN) :: buff
        ! splicer begin class.View.method.attach_buffer_only
        call c_view_attach_buffer_only(obj%voidptr, buff%get_instance())
        ! splicer end class.View.method.attach_buffer_only
    end subroutine view_attach_buffer_only

    subroutine view_attach_buffer_type_int(obj, type, num_elems, buff)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: num_elems
        type(SidreBuffer), value, intent(IN) :: buff
        ! splicer begin class.View.method.attach_buffer_type_int
        call c_view_attach_buffer_type(obj%voidptr, type, &
            int(num_elems, C_LONG), buff%get_instance())
        ! splicer end class.View.method.attach_buffer_type_int
    end subroutine view_attach_buffer_type_int

    subroutine view_attach_buffer_type_long(obj, type, num_elems, buff)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        type(SidreBuffer), value, intent(IN) :: buff
        ! splicer begin class.View.method.attach_buffer_type_long
        call c_view_attach_buffer_type(obj%voidptr, type, &
            int(num_elems, C_LONG), buff%get_instance())
        ! splicer end class.View.method.attach_buffer_type_long
    end subroutine view_attach_buffer_type_long

    subroutine view_attach_buffer_shape(obj, type, ndims, shape, buff)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: ndims
        integer(C_LONG), intent(IN) :: shape(*)
        type(SidreBuffer), value, intent(IN) :: buff
        ! splicer begin class.View.method.attach_buffer_shape
        call c_view_attach_buffer_shape(obj%voidptr, type, ndims, shape, &
            buff%get_instance())
        ! splicer end class.View.method.attach_buffer_shape
    end subroutine view_attach_buffer_shape

    subroutine view_apply_nelems(obj, num_elems)
        use iso_c_binding, only : C_LONG
        class(SidreView) :: obj
        integer(C_LONG), value, intent(IN) :: num_elems
        ! splicer begin class.View.method.apply_nelems
        call c_view_apply_nelems(obj%voidptr, num_elems)
        ! splicer end class.View.method.apply_nelems
    end subroutine view_apply_nelems

    subroutine view_apply_nelems_offset(obj, num_elems, offset)
        use iso_c_binding, only : C_LONG
        class(SidreView) :: obj
        integer(C_LONG), value, intent(IN) :: num_elems
        integer(C_LONG), value, intent(IN) :: offset
        ! splicer begin class.View.method.apply_nelems_offset
        call c_view_apply_nelems_offset(obj%voidptr, num_elems, offset)
        ! splicer end class.View.method.apply_nelems_offset
    end subroutine view_apply_nelems_offset

    subroutine view_apply_nelems_offset_stride(obj, num_elems, offset, &
            stride)
        use iso_c_binding, only : C_LONG
        class(SidreView) :: obj
        integer(C_LONG), value, intent(IN) :: num_elems
        integer(C_LONG), value, intent(IN) :: offset
        integer(C_LONG), value, intent(IN) :: stride
        ! splicer begin class.View.method.apply_nelems_offset_stride
        call c_view_apply_nelems_offset_stride(obj%voidptr, num_elems, &
            offset, stride)
        ! splicer end class.View.method.apply_nelems_offset_stride
    end subroutine view_apply_nelems_offset_stride

    subroutine view_apply_type_nelems(obj, type, num_elems)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        ! splicer begin class.View.method.apply_type_nelems
        call c_view_apply_type_nelems(obj%voidptr, type, num_elems)
        ! splicer end class.View.method.apply_type_nelems
    end subroutine view_apply_type_nelems

    subroutine view_apply_type_nelems_offset(obj, type, num_elems, &
            offset)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        integer(C_LONG), value, intent(IN) :: offset
        ! splicer begin class.View.method.apply_type_nelems_offset
        call c_view_apply_type_nelems_offset(obj%voidptr, type, &
            num_elems, offset)
        ! splicer end class.View.method.apply_type_nelems_offset
    end subroutine view_apply_type_nelems_offset

    subroutine view_apply_type_nelems_offset_stride(obj, type, &
            num_elems, offset, stride)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        integer(C_LONG), value, intent(IN) :: offset
        integer(C_LONG), value, intent(IN) :: stride
        ! splicer begin class.View.method.apply_type_nelems_offset_stride
        call c_view_apply_type_nelems_offset_stride(obj%voidptr, type, &
            num_elems, offset, stride)
        ! splicer end class.View.method.apply_type_nelems_offset_stride
    end subroutine view_apply_type_nelems_offset_stride

    subroutine view_apply_type_shape(obj, type, ndims, shape)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: ndims
        integer(C_LONG), intent(IN) :: shape(*)
        ! splicer begin class.View.method.apply_type_shape
        call c_view_apply_type_shape(obj%voidptr, type, ndims, shape)
        ! splicer end class.View.method.apply_type_shape
    end subroutine view_apply_type_shape

    function view_has_buffer(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        class(SidreView) :: obj
        logical :: SHT_rv
        ! splicer begin class.View.method.has_buffer
        SHT_rv = c_view_has_buffer(obj%voidptr)
        ! splicer end class.View.method.has_buffer
    end function view_has_buffer

    function view_is_external(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        class(SidreView) :: obj
        logical :: SHT_rv
        ! splicer begin class.View.method.is_external
        SHT_rv = c_view_is_external(obj%voidptr)
        ! splicer end class.View.method.is_external
    end function view_is_external

    function view_is_allocated(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        class(SidreView) :: obj
        logical :: SHT_rv
        ! splicer begin class.View.method.is_allocated
        SHT_rv = c_view_is_allocated(obj%voidptr)
        ! splicer end class.View.method.is_allocated
    end function view_is_allocated

    function view_is_applied(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        class(SidreView) :: obj
        logical :: SHT_rv
        ! splicer begin class.View.method.is_applied
        SHT_rv = c_view_is_applied(obj%voidptr)
        ! splicer end class.View.method.is_applied
    end function view_is_applied

    function view_is_described(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        class(SidreView) :: obj
        logical :: SHT_rv
        ! splicer begin class.View.method.is_described
        SHT_rv = c_view_is_described(obj%voidptr)
        ! splicer end class.View.method.is_described
    end function view_is_described

    function view_is_empty(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        class(SidreView) :: obj
        logical :: SHT_rv
        ! splicer begin class.View.method.is_empty
        SHT_rv = c_view_is_empty(obj%voidptr)
        ! splicer end class.View.method.is_empty
    end function view_is_empty

    function view_is_opaque(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        class(SidreView) :: obj
        logical :: SHT_rv
        ! splicer begin class.View.method.is_opaque
        SHT_rv = c_view_is_opaque(obj%voidptr)
        ! splicer end class.View.method.is_opaque
    end function view_is_opaque

    function view_is_scalar(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        class(SidreView) :: obj
        logical :: SHT_rv
        ! splicer begin class.View.method.is_scalar
        SHT_rv = c_view_is_scalar(obj%voidptr)
        ! splicer end class.View.method.is_scalar
    end function view_is_scalar

    function view_is_string(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        class(SidreView) :: obj
        logical :: SHT_rv
        ! splicer begin class.View.method.is_string
        SHT_rv = c_view_is_string(obj%voidptr)
        ! splicer end class.View.method.is_string
    end function view_is_string

    function view_get_index(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreView) :: obj
        integer(C_INT) :: SHT_rv
        ! splicer begin class.View.method.get_index
        SHT_rv = c_view_get_index(obj%voidptr)
        ! splicer end class.View.method.get_index
    end function view_get_index

    function view_get_name(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        class(SidreView) :: obj
        character(kind=C_CHAR, len=MAXNAMESIZE) :: SHT_rv
        ! splicer begin class.View.method.get_name
        call c_view_get_name_bufferify(obj%voidptr, SHT_rv, &
            len(SHT_rv, kind=C_INT))
        ! splicer end class.View.method.get_name
    end function view_get_name

    function view_get_path(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        class(SidreView) :: obj
        character(kind=C_CHAR, len=MAXNAMESIZE) :: SHT_rv
        ! splicer begin class.View.method.get_path
        call c_view_get_path_bufferify(obj%voidptr, SHT_rv, &
            len(SHT_rv, kind=C_INT))
        ! splicer end class.View.method.get_path
    end function view_get_path

    function view_get_path_name(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_CHAR, C_INT
        class(SidreView) :: obj
        character(kind=C_CHAR, len=MAXNAMESIZE) :: SHT_rv
        ! splicer begin class.View.method.get_path_name
        call c_view_get_path_name_bufferify(obj%voidptr, SHT_rv, &
            len(SHT_rv, kind=C_INT))
        ! splicer end class.View.method.get_path_name
    end function view_get_path_name

    function view_get_buffer(obj) &
            result(SHT_rv)
        class(SidreView) :: obj
        type(SidreBuffer) :: SHT_rv
        ! splicer begin class.View.method.get_buffer
        SHT_rv%voidptr = c_view_get_buffer(obj%voidptr)
        ! splicer end class.View.method.get_buffer
    end function view_get_buffer

    function view_get_void_ptr(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_PTR
        class(SidreView) :: obj
        type(C_PTR) :: SHT_rv
        ! splicer begin class.View.method.get_void_ptr
        SHT_rv = c_view_get_void_ptr(obj%voidptr)
        ! splicer end class.View.method.get_void_ptr
    end function view_get_void_ptr

    subroutine view_set_scalar_int(obj, value)
        use iso_c_binding, only : C_INT
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: value
        ! splicer begin class.View.method.set_scalar_int
        call c_view_set_scalar_int(obj%voidptr, value)
        ! splicer end class.View.method.set_scalar_int
    end subroutine view_set_scalar_int

    subroutine view_set_scalar_long(obj, value)
        use iso_c_binding, only : C_LONG
        class(SidreView) :: obj
        integer(C_LONG), value, intent(IN) :: value
        ! splicer begin class.View.method.set_scalar_long
        call c_view_set_scalar_long(obj%voidptr, value)
        ! splicer end class.View.method.set_scalar_long
    end subroutine view_set_scalar_long

    subroutine view_set_scalar_float(obj, value)
        use iso_c_binding, only : C_FLOAT
        class(SidreView) :: obj
        real(C_FLOAT), value, intent(IN) :: value
        ! splicer begin class.View.method.set_scalar_float
        call c_view_set_scalar_float(obj%voidptr, value)
        ! splicer end class.View.method.set_scalar_float
    end subroutine view_set_scalar_float

    subroutine view_set_scalar_double(obj, value)
        use iso_c_binding, only : C_DOUBLE
        class(SidreView) :: obj
        real(C_DOUBLE), value, intent(IN) :: value
        ! splicer begin class.View.method.set_scalar_double
        call c_view_set_scalar_double(obj%voidptr, value)
        ! splicer end class.View.method.set_scalar_double
    end subroutine view_set_scalar_double

    subroutine view_set_external_data_ptr_only(obj, external_ptr)
        use iso_c_binding, only : C_PTR
        class(SidreView) :: obj
        type(C_PTR), value, intent(IN) :: external_ptr
        ! splicer begin class.View.method.set_external_data_ptr_only
        call c_view_set_external_data_ptr_only(obj%voidptr, &
            external_ptr)
        ! splicer end class.View.method.set_external_data_ptr_only
    end subroutine view_set_external_data_ptr_only

    subroutine view_set_external_data_ptr_type_int(obj, type, num_elems, &
            external_ptr)
        use iso_c_binding, only : C_INT, C_LONG, C_PTR
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: num_elems
        type(C_PTR), value, intent(IN) :: external_ptr
        ! splicer begin class.View.method.set_external_data_ptr_type_int
        call c_view_set_external_data_ptr_type(obj%voidptr, type, &
            int(num_elems, C_LONG), external_ptr)
        ! splicer end class.View.method.set_external_data_ptr_type_int
    end subroutine view_set_external_data_ptr_type_int

    subroutine view_set_external_data_ptr_type_long(obj, type, &
            num_elems, external_ptr)
        use iso_c_binding, only : C_INT, C_LONG, C_PTR
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_LONG), value, intent(IN) :: num_elems
        type(C_PTR), value, intent(IN) :: external_ptr
        ! splicer begin class.View.method.set_external_data_ptr_type_long
        call c_view_set_external_data_ptr_type(obj%voidptr, type, &
            int(num_elems, C_LONG), external_ptr)
        ! splicer end class.View.method.set_external_data_ptr_type_long
    end subroutine view_set_external_data_ptr_type_long

    subroutine view_set_string(obj, value)
        use iso_c_binding, only : C_INT
        class(SidreView) :: obj
        character(*), intent(IN) :: value
        ! splicer begin class.View.method.set_string
        call c_view_set_string_bufferify(obj%voidptr, value, &
            len_trim(value, kind=C_INT))
        ! splicer end class.View.method.set_string
    end subroutine view_set_string

    subroutine view_set_external_data_ptr_shape(obj, type, ndims, shape, &
            external_ptr)
        use iso_c_binding, only : C_INT, C_LONG, C_PTR
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: type
        integer(C_INT), value, intent(IN) :: ndims
        integer(C_LONG), intent(IN) :: shape(*)
        type(C_PTR), value, intent(IN) :: external_ptr
        ! splicer begin class.View.method.set_external_data_ptr_shape
        call c_view_set_external_data_ptr_shape(obj%voidptr, type, &
            ndims, shape, external_ptr)
        ! splicer end class.View.method.set_external_data_ptr_shape
    end subroutine view_set_external_data_ptr_shape

    subroutine view_get_string(obj, name)
        use iso_c_binding, only : C_INT
        class(SidreView) :: obj
        character(*), intent(OUT) :: name
        ! splicer begin class.View.method.get_string
        call c_view_get_string_bufferify(obj%voidptr, name, &
            len(name, kind=C_INT))
        ! splicer end class.View.method.get_string
    end subroutine view_get_string

    function view_get_data_int(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreView) :: obj
        integer(C_INT) :: SHT_rv
        ! splicer begin class.View.method.get_data_int
        SHT_rv = c_view_get_data_int(obj%voidptr)
        ! splicer end class.View.method.get_data_int
    end function view_get_data_int

    function view_get_data_long(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_LONG
        class(SidreView) :: obj
        integer(C_LONG) :: SHT_rv
        ! splicer begin class.View.method.get_data_long
        SHT_rv = c_view_get_data_long(obj%voidptr)
        ! splicer end class.View.method.get_data_long
    end function view_get_data_long

    function view_get_data_float(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_FLOAT
        class(SidreView) :: obj
        real(C_FLOAT) :: SHT_rv
        ! splicer begin class.View.method.get_data_float
        SHT_rv = c_view_get_data_float(obj%voidptr)
        ! splicer end class.View.method.get_data_float
    end function view_get_data_float

    function view_get_data_double(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_DOUBLE
        class(SidreView) :: obj
        real(C_DOUBLE) :: SHT_rv
        ! splicer begin class.View.method.get_data_double
        SHT_rv = c_view_get_data_double(obj%voidptr)
        ! splicer end class.View.method.get_data_double
    end function view_get_data_double

    function view_get_owning_group(obj) &
            result(SHT_rv)
        class(SidreView) :: obj
        type(SidreGroup) :: SHT_rv
        ! splicer begin class.View.method.get_owning_group
        SHT_rv%voidptr = c_view_get_owning_group(obj%voidptr)
        ! splicer end class.View.method.get_owning_group
    end function view_get_owning_group

    function view_get_type_id(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreView) :: obj
        integer(C_INT) :: SHT_rv
        ! splicer begin class.View.method.get_type_id
        SHT_rv = c_view_get_type_id(obj%voidptr)
        ! splicer end class.View.method.get_type_id
    end function view_get_type_id

    function view_get_total_bytes(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_SIZE_T
        class(SidreView) :: obj
        integer(C_SIZE_T) :: SHT_rv
        ! splicer begin class.View.method.get_total_bytes
        SHT_rv = c_view_get_total_bytes(obj%voidptr)
        ! splicer end class.View.method.get_total_bytes
    end function view_get_total_bytes

    function view_get_bytes_per_element(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_SIZE_T
        class(SidreView) :: obj
        integer(C_SIZE_T) :: SHT_rv
        ! splicer begin class.View.method.get_bytes_per_element
        SHT_rv = c_view_get_bytes_per_element(obj%voidptr)
        ! splicer end class.View.method.get_bytes_per_element
    end function view_get_bytes_per_element

    function view_get_num_elements(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_SIZE_T
        class(SidreView) :: obj
        integer(C_SIZE_T) :: SHT_rv
        ! splicer begin class.View.method.get_num_elements
        SHT_rv = c_view_get_num_elements(obj%voidptr)
        ! splicer end class.View.method.get_num_elements
    end function view_get_num_elements

    function view_get_offset(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_SIZE_T
        class(SidreView) :: obj
        integer(C_SIZE_T) :: SHT_rv
        ! splicer begin class.View.method.get_offset
        SHT_rv = c_view_get_offset(obj%voidptr)
        ! splicer end class.View.method.get_offset
    end function view_get_offset

    function view_get_stride(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_SIZE_T
        class(SidreView) :: obj
        integer(C_SIZE_T) :: SHT_rv
        ! splicer begin class.View.method.get_stride
        SHT_rv = c_view_get_stride(obj%voidptr)
        ! splicer end class.View.method.get_stride
    end function view_get_stride

    function view_get_num_dimensions(obj) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT
        class(SidreView) :: obj
        integer(C_INT) :: SHT_rv
        ! splicer begin class.View.method.get_num_dimensions
        SHT_rv = c_view_get_num_dimensions(obj%voidptr)
        ! splicer end class.View.method.get_num_dimensions
    end function view_get_num_dimensions

    function view_get_shape(obj, ndims, shape) &
            result(SHT_rv)
        use iso_c_binding, only : C_INT, C_LONG
        class(SidreView) :: obj
        integer(C_INT), value, intent(IN) :: ndims
        integer(C_LONG), intent(OUT) :: shape(*)
        integer(C_INT) :: SHT_rv
        ! splicer begin class.View.method.get_shape
        SHT_rv = c_view_get_shape(obj%voidptr, ndims, shape)
        ! splicer end class.View.method.get_shape
    end function view_get_shape

    function view_rename(obj, new_name) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL, C_INT
        class(SidreView) :: obj
        character(*), intent(IN) :: new_name
        logical :: SHT_rv
        ! splicer begin class.View.method.rename
        SHT_rv = c_view_rename_bufferify(obj%voidptr, new_name, &
            len_trim(new_name, kind=C_INT))
        ! splicer end class.View.method.rename
    end function view_rename

    subroutine view_print(obj)
        class(SidreView) :: obj
        ! splicer begin class.View.method.print
        call c_view_print(obj%voidptr)
        ! splicer end class.View.method.print
    end subroutine view_print

    function view_get_instance(obj) result (voidptr)
        use iso_c_binding, only: C_PTR
        class(SidreView), intent(IN) :: obj
        type(C_PTR) :: voidptr
        voidptr = obj%voidptr
    end function view_get_instance

    subroutine view_set_instance(obj, voidptr)
        use iso_c_binding, only: C_PTR
        class(SidreView), intent(INOUT) :: obj
        type(C_PTR), intent(IN) :: voidptr
        obj%voidptr = voidptr
    end subroutine view_set_instance

    function view_associated(obj) result (rv)
        use iso_c_binding, only: c_associated
        class(SidreView), intent(IN) :: obj
        logical rv
        rv = c_associated(obj%voidptr)
    end function view_associated

    ! splicer begin class.View.additional_functions

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_int_scalar_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        integer(C_INT), pointer, intent(OUT) :: value
        integer(C_INT), pointer :: tmp(:)
        type(C_PTR) cptr
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          call c_f_pointer(cptr, value)
        else
          nullify(value)
        endif
    end subroutine view_get_data_int_scalar_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_int_1d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        integer(C_INT), pointer, intent(OUT) :: value(:)
        integer(C_INT), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(1)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(1, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_int_1d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_int_2d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        integer(C_INT), pointer, intent(OUT) :: value(:,:)
        integer(C_INT), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(2)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(2, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_int_2d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_int_3d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        integer(C_INT), pointer, intent(OUT) :: value(:,:,:)
        integer(C_INT), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(3)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(3, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_int_3d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_int_4d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        integer(C_INT), pointer, intent(OUT) :: value(:,:,:,:)
        integer(C_INT), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(4)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(4, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_int_4d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_long_scalar_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        integer(C_LONG), pointer, intent(OUT) :: value
        integer(C_LONG), pointer :: tmp(:)
        type(C_PTR) cptr
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          call c_f_pointer(cptr, value)
        else
          nullify(value)
        endif
    end subroutine view_get_data_long_scalar_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_long_1d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        integer(C_LONG), pointer, intent(OUT) :: value(:)
        integer(C_LONG), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(1)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(1, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_long_1d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_long_2d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        integer(C_LONG), pointer, intent(OUT) :: value(:,:)
        integer(C_LONG), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(2)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(2, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_long_2d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_long_3d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        integer(C_LONG), pointer, intent(OUT) :: value(:,:,:)
        integer(C_LONG), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(3)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(3, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_long_3d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_long_4d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        integer(C_LONG), pointer, intent(OUT) :: value(:,:,:,:)
        integer(C_LONG), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(4)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(4, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_long_4d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_float_scalar_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        real(C_FLOAT), pointer, intent(OUT) :: value
        real(C_FLOAT), pointer :: tmp(:)
        type(C_PTR) cptr
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          call c_f_pointer(cptr, value)
        else
          nullify(value)
        endif
    end subroutine view_get_data_float_scalar_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_float_1d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        real(C_FLOAT), pointer, intent(OUT) :: value(:)
        real(C_FLOAT), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(1)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(1, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_float_1d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_float_2d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        real(C_FLOAT), pointer, intent(OUT) :: value(:,:)
        real(C_FLOAT), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(2)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(2, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_float_2d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_float_3d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        real(C_FLOAT), pointer, intent(OUT) :: value(:,:,:)
        real(C_FLOAT), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(3)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(3, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_float_3d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_float_4d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        real(C_FLOAT), pointer, intent(OUT) :: value(:,:,:,:)
        real(C_FLOAT), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(4)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(4, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_float_4d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_double_scalar_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        real(C_DOUBLE), pointer, intent(OUT) :: value
        real(C_DOUBLE), pointer :: tmp(:)
        type(C_PTR) cptr
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          call c_f_pointer(cptr, value)
        else
          nullify(value)
        endif
    end subroutine view_get_data_double_scalar_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_double_1d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        real(C_DOUBLE), pointer, intent(OUT) :: value(:)
        real(C_DOUBLE), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(1)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(1, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_double_1d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_double_2d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        real(C_DOUBLE), pointer, intent(OUT) :: value(:,:)
        real(C_DOUBLE), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(2)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(2, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_double_2d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_double_3d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        real(C_DOUBLE), pointer, intent(OUT) :: value(:,:,:)
        real(C_DOUBLE), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(3)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(3, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_double_3d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_get_data_double_4d_ptr(view, value)
        use iso_c_binding
        implicit none
        class(SidreView), intent(IN) :: view
        real(C_DOUBLE), pointer, intent(OUT) :: value(:,:,:,:)
        real(C_DOUBLE), pointer :: tmp(:)
        type(C_PTR) cptr
        integer rank
        integer(SIDRE_LENGTH) extents(4)
        integer(SIDRE_LENGTH) :: offset

        cptr = view%get_void_ptr()
        if (c_associated(cptr)) then
          offset = view%get_offset()
          if (offset > 0) then
            call c_f_pointer(cptr, tmp, [offset+1])   ! +1 to convert 0-based offset to 1-based index
            cptr = c_loc(tmp(offset+1))               ! Emulates pointer arithmetic
          endif

          rank = view%get_shape(4, extents)
          call c_f_pointer(cptr, value, extents)
        else
          nullify(value)
        endif
    end subroutine view_get_data_double_4d_ptr

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_int_scalar(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        integer(C_INT), target, intent(IN) :: value
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents(1) = 1_SIDRE_LENGTH
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value, addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 0, extents)
    end subroutine view_set_array_data_ptr_int_scalar

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_int_1d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        integer(C_INT), target, intent(IN) :: value(:)
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 1, extents)
    end subroutine view_set_array_data_ptr_int_1d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_int_2d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        integer(C_INT), target, intent(IN) :: value(:,:)
    !    integer(SIDRE_LENGTH) :: extents(2)
    !    integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 2, extents)
    end subroutine view_set_array_data_ptr_int_2d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_int_3d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        integer(C_INT), target, intent(IN) :: value(:,:,:)
    !    integer(SIDRE_LENGTH) :: extents(3)
    !    integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 3, extents)
    end subroutine view_set_array_data_ptr_int_3d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_int_4d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        integer(C_INT), target, intent(IN) :: value(:,:,:,:)
    !    integer(SIDRE_LENGTH) :: extents(4)
    !    integer(C_INT), parameter :: type = SIDRE_INT_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3),lbound(value,4)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 4, extents)
    end subroutine view_set_array_data_ptr_int_4d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_long_scalar(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        integer(C_LONG), target, intent(IN) :: value
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents(1) = 1_SIDRE_LENGTH
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value, addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 0, extents)
    end subroutine view_set_array_data_ptr_long_scalar

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_long_1d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        integer(C_LONG), target, intent(IN) :: value(:)
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 1, extents)
    end subroutine view_set_array_data_ptr_long_1d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_long_2d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        integer(C_LONG), target, intent(IN) :: value(:,:)
    !    integer(SIDRE_LENGTH) :: extents(2)
    !    integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 2, extents)
    end subroutine view_set_array_data_ptr_long_2d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_long_3d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        integer(C_LONG), target, intent(IN) :: value(:,:,:)
    !    integer(SIDRE_LENGTH) :: extents(3)
    !    integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 3, extents)
    end subroutine view_set_array_data_ptr_long_3d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_long_4d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        integer(C_LONG), target, intent(IN) :: value(:,:,:,:)
    !    integer(SIDRE_LENGTH) :: extents(4)
    !    integer(C_INT), parameter :: type = SIDRE_LONG_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3),lbound(value,4)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 4, extents)
    end subroutine view_set_array_data_ptr_long_4d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_float_scalar(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        real(C_FLOAT), target, intent(IN) :: value
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents(1) = 1_SIDRE_LENGTH
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value, addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 0, extents)
    end subroutine view_set_array_data_ptr_float_scalar

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_float_1d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        real(C_FLOAT), target, intent(IN) :: value(:)
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 1, extents)
    end subroutine view_set_array_data_ptr_float_1d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_float_2d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        real(C_FLOAT), target, intent(IN) :: value(:,:)
    !    integer(SIDRE_LENGTH) :: extents(2)
    !    integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 2, extents)
    end subroutine view_set_array_data_ptr_float_2d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_float_3d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        real(C_FLOAT), target, intent(IN) :: value(:,:,:)
    !    integer(SIDRE_LENGTH) :: extents(3)
    !    integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 3, extents)
    end subroutine view_set_array_data_ptr_float_3d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_float_4d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        real(C_FLOAT), target, intent(IN) :: value(:,:,:,:)
    !    integer(SIDRE_LENGTH) :: extents(4)
    !    integer(C_INT), parameter :: type = SIDRE_FLOAT_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3),lbound(value,4)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 4, extents)
    end subroutine view_set_array_data_ptr_float_4d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_double_scalar(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        real(C_DOUBLE), target, intent(IN) :: value
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents(1) = 1_SIDRE_LENGTH
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value, addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 0, extents)
    end subroutine view_set_array_data_ptr_double_scalar

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_double_1d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        real(C_DOUBLE), target, intent(IN) :: value(:)
    !    integer(SIDRE_LENGTH) :: extents(1)
    !    integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 1, extents)
    end subroutine view_set_array_data_ptr_double_1d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_double_2d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        real(C_DOUBLE), target, intent(IN) :: value(:,:)
    !    integer(SIDRE_LENGTH) :: extents(2)
    !    integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 2, extents)
    end subroutine view_set_array_data_ptr_double_2d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_double_3d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        real(C_DOUBLE), target, intent(IN) :: value(:,:,:)
    !    integer(SIDRE_LENGTH) :: extents(3)
    !    integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 3, extents)
    end subroutine view_set_array_data_ptr_double_3d

    ! Generated by genfsidresplicer.py
    subroutine view_set_array_data_ptr_double_4d(view, value)
        use iso_c_binding
        implicit none

        class(SidreView), intent(IN) :: view
        real(C_DOUBLE), target, intent(IN) :: value(:,:,:,:)
    !    integer(SIDRE_LENGTH) :: extents(4)
    !    integer(C_INT), parameter :: type = SIDRE_DOUBLE_ID
        type(C_PTR) addr

    !    lname = len_trim(name)
    !    extents = shape(value, kind=SIDRE_LENGTH)
#ifdef USE_C_LOC_WITH_ASSUMED_SHAPE
        addr = c_loc(value)
#else
        call SIDRE_C_LOC(value(lbound(value,1),lbound(value,2),lbound(value,3),lbound(value,4)), addr)
#endif
        call c_view_set_external_data_ptr_only(view%voidptr, addr)
    !    call c_view_apply_type_shape(rv%voidptr, type, 4, extents)
    end subroutine view_set_array_data_ptr_double_4d
    ! splicer end class.View.additional_functions

    function name_is_valid(name) &
            result(SHT_rv)
        use iso_c_binding, only : C_BOOL
        character(*), intent(IN) :: name
        logical :: SHT_rv
        ! splicer begin function.name_is_valid
        SHT_rv = name .ne. " "
        ! splicer end function.name_is_valid
    end function name_is_valid

    ! splicer begin additional_functions
    ! splicer end additional_functions

    function datastore_eq(a,b) result (rv)
        use iso_c_binding, only: c_associated
        type(SidreDataStore), intent(IN) ::a,b
        logical :: rv
        if (c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function datastore_eq

    function datastore_ne(a,b) result (rv)
        use iso_c_binding, only: c_associated
        type(SidreDataStore), intent(IN) ::a,b
        logical :: rv
        if (.not. c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function datastore_ne

    function group_eq(a,b) result (rv)
        use iso_c_binding, only: c_associated
        type(SidreGroup), intent(IN) ::a,b
        logical :: rv
        if (c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function group_eq

    function group_ne(a,b) result (rv)
        use iso_c_binding, only: c_associated
        type(SidreGroup), intent(IN) ::a,b
        logical :: rv
        if (.not. c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function group_ne

    function buffer_eq(a,b) result (rv)
        use iso_c_binding, only: c_associated
        type(SidreBuffer), intent(IN) ::a,b
        logical :: rv
        if (c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function buffer_eq

    function buffer_ne(a,b) result (rv)
        use iso_c_binding, only: c_associated
        type(SidreBuffer), intent(IN) ::a,b
        logical :: rv
        if (.not. c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function buffer_ne

    function view_eq(a,b) result (rv)
        use iso_c_binding, only: c_associated
        type(SidreView), intent(IN) ::a,b
        logical :: rv
        if (c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function view_eq

    function view_ne(a,b) result (rv)
        use iso_c_binding, only: c_associated
        type(SidreView), intent(IN) ::a,b
        logical :: rv
        if (.not. c_associated(a%voidptr, b%voidptr)) then
            rv = .true.
        else
            rv = .false.
        endif
    end function view_ne

end module axom_sidre

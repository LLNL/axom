// wrapGroup.h
// This is generated code, do not edit
//
// Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
//
// Produced at the Lawrence Livermore National Laboratory
//
// LLNL-CODE-741217
//
// All rights reserved.
//
// This file is part of Axom.
//
// For details about use and distribution, please read axom/LICENSE.
//
/**
 * \file wrapGroup.h
 * \brief Shroud generated wrapper for Group class
 */
// For C users and C++ implementation

#ifndef WRAPGROUP_H
#define WRAPGROUP_H

#include <stddef.h>
#include "axom/sidre/interface/SidreTypes.h"
#include "typesSidre.h"

// splicer begin class.Group.CXX_declarations
// splicer end class.Group.CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// splicer begin class.Group.C_declarations
// splicer end class.Group.C_declarations

SIDRE_IndexType SIDRE_group_get_index(SIDRE_group * self);

const char * SIDRE_group_get_name(const SIDRE_group * self);

void SIDRE_group_get_name_bufferify(const SIDRE_group * self, char * SHF_rv, int NSHF_rv);

void SIDRE_group_get_path_bufferify(const SIDRE_group * self, char * SHF_rv, int NSHF_rv);

void SIDRE_group_get_path_name_bufferify(const SIDRE_group * self, char * SHF_rv, int NSHF_rv);

SIDRE_group * SIDRE_group_get_parent(const SIDRE_group * self, SIDRE_group * SHC_rv);

size_t SIDRE_group_get_num_groups(const SIDRE_group * self);

size_t SIDRE_group_get_num_views(const SIDRE_group * self);

SIDRE_datastore * SIDRE_group_get_data_store(const SIDRE_group * self, SIDRE_datastore * SHC_rv);

bool SIDRE_group_has_view(const SIDRE_group * self, const char * path);

bool SIDRE_group_has_view_bufferify(const SIDRE_group * self, const char * path, int Lpath);

bool SIDRE_group_has_child_view(const SIDRE_group * self, const char * name);

bool SIDRE_group_has_child_view_bufferify(const SIDRE_group * self, const char * name, int Lname);

SIDRE_IndexType SIDRE_group_get_view_index(const SIDRE_group * self, const char * name);

SIDRE_IndexType SIDRE_group_get_view_index_bufferify(const SIDRE_group * self, const char * name, int Lname);

const char * SIDRE_group_get_view_name(const SIDRE_group * self, SIDRE_IndexType idx);

void SIDRE_group_get_view_name_bufferify(const SIDRE_group * self, SIDRE_IndexType idx, char * SHF_rv, int NSHF_rv);

SIDRE_view * SIDRE_group_get_view_from_name(SIDRE_group * self, const char * path, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_get_view_from_name_bufferify(SIDRE_group * self, const char * path, int Lpath, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_get_view_from_index(SIDRE_group * self, const SIDRE_IndexType idx, SIDRE_view * SHC_rv);

SIDRE_IndexType SIDRE_group_get_first_valid_view_index(const SIDRE_group * self);

SIDRE_IndexType SIDRE_group_get_next_valid_view_index(const SIDRE_group * self, SIDRE_IndexType idx);

SIDRE_view * SIDRE_group_create_view_empty(SIDRE_group * self, const char * path, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_empty_bufferify(SIDRE_group * self, const char * path, int Lpath, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_from_type(SIDRE_group * self, const char * path, int type, SIDRE_IndexType num_elems, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_from_type_bufferify(SIDRE_group * self, const char * path, int Lpath, int type, SIDRE_IndexType num_elems, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_from_shape(SIDRE_group * self, const char * path, int type, int ndims, SIDRE_IndexType * shape, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_from_shape_bufferify(SIDRE_group * self, const char * path, int Lpath, int type, int ndims, SIDRE_IndexType * shape, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_into_buffer(SIDRE_group * self, const char * path, SIDRE_buffer * buff, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_into_buffer_bufferify(SIDRE_group * self, const char * path, int Lpath, SIDRE_buffer * buff, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_from_type_and_buffer(SIDRE_group * self, const char * path, int type, SIDRE_IndexType num_elems, SIDRE_buffer * buff, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_from_type_and_buffer_bufferify(SIDRE_group * self, const char * path, int Lpath, int type, SIDRE_IndexType num_elems, SIDRE_buffer * buff, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_from_shape_and_buffer(SIDRE_group * self, const char * path, int type, int ndims, SIDRE_IndexType * shape, SIDRE_buffer * buff, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_from_shape_and_buffer_bufferify(SIDRE_group * self, const char * path, int Lpath, int type, int ndims, SIDRE_IndexType * shape, SIDRE_buffer * buff, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_external(SIDRE_group * self, const char * path, void * external_ptr, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_external_bufferify(SIDRE_group * self, const char * path, int Lpath, void * external_ptr, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_from_type_external(SIDRE_group * self, const char * path, int type, SIDRE_IndexType num_elems, void * external_ptr, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_from_type_external_bufferify(SIDRE_group * self, const char * path, int Lpath, int type, SIDRE_IndexType num_elems, void * external_ptr, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_from_shape_external(SIDRE_group * self, const char * path, int type, int ndims, SIDRE_IndexType * shape, void * external_ptr, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_from_shape_external_bufferify(SIDRE_group * self, const char * path, int Lpath, int type, int ndims, SIDRE_IndexType * shape, void * external_ptr, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_and_allocate_nelems(SIDRE_group * self, const char * path, int type, SIDRE_IndexType num_elems, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_and_allocate_nelems_bufferify(SIDRE_group * self, const char * path, int Lpath, int type, SIDRE_IndexType num_elems, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_and_allocate_shape(SIDRE_group * self, const char * path, int type, int ndims, SIDRE_IndexType * shape, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_and_allocate_shape_bufferify(SIDRE_group * self, const char * path, int Lpath, int type, int ndims, SIDRE_IndexType * shape, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_scalar_int(SIDRE_group * self, const char * path, int value, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_scalar_int_bufferify(SIDRE_group * self, const char * path, int Lpath, int value, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_scalar_long(SIDRE_group * self, const char * path, long value, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_scalar_long_bufferify(SIDRE_group * self, const char * path, int Lpath, long value, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_scalar_float(SIDRE_group * self, const char * path, float value, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_scalar_float_bufferify(SIDRE_group * self, const char * path, int Lpath, float value, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_scalar_double(SIDRE_group * self, const char * path, double value, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_scalar_double_bufferify(SIDRE_group * self, const char * path, int Lpath, double value, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_string(SIDRE_group * self, const char * path, const char * value, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_create_view_string_bufferify(SIDRE_group * self, const char * path, int Lpath, const char * value, int Lvalue, SIDRE_view * SHC_rv);

void SIDRE_group_destroy_view(SIDRE_group * self, const char * path);

void SIDRE_group_destroy_view_bufferify(SIDRE_group * self, const char * path, int Lpath);

void SIDRE_group_destroy_view_and_data_name(SIDRE_group * self, const char * path);

void SIDRE_group_destroy_view_and_data_name_bufferify(SIDRE_group * self, const char * path, int Lpath);

void SIDRE_group_destroy_view_and_data_index(SIDRE_group * self, SIDRE_IndexType idx);

SIDRE_view * SIDRE_group_move_view(SIDRE_group * self, SIDRE_view * view, SIDRE_view * SHC_rv);

SIDRE_view * SIDRE_group_copy_view(SIDRE_group * self, SIDRE_view * view, SIDRE_view * SHC_rv);

bool SIDRE_group_has_group(SIDRE_group * self, const char * path);

bool SIDRE_group_has_group_bufferify(SIDRE_group * self, const char * path, int Lpath);

bool SIDRE_group_has_child_group(SIDRE_group * self, const char * name);

bool SIDRE_group_has_child_group_bufferify(SIDRE_group * self, const char * name, int Lname);

SIDRE_IndexType SIDRE_group_get_group_index(const SIDRE_group * self, const char * name);

SIDRE_IndexType SIDRE_group_get_group_index_bufferify(const SIDRE_group * self, const char * name, int Lname);

const char * SIDRE_group_get_group_name(const SIDRE_group * self, SIDRE_IndexType idx);

void SIDRE_group_get_group_name_bufferify(const SIDRE_group * self, SIDRE_IndexType idx, char * SHF_rv, int NSHF_rv);

SIDRE_group * SIDRE_group_get_group_from_name(SIDRE_group * self, const char * path, SIDRE_group * SHC_rv);

SIDRE_group * SIDRE_group_get_group_from_name_bufferify(SIDRE_group * self, const char * path, int Lpath, SIDRE_group * SHC_rv);

SIDRE_group * SIDRE_group_get_group_from_index(SIDRE_group * self, SIDRE_IndexType idx, SIDRE_group * SHC_rv);

SIDRE_IndexType SIDRE_group_get_first_valid_group_index(const SIDRE_group * self);

SIDRE_IndexType SIDRE_group_get_next_valid_group_index(const SIDRE_group * self, SIDRE_IndexType idx);

SIDRE_group * SIDRE_group_create_group(SIDRE_group * self, const char * path, SIDRE_group * SHC_rv);

SIDRE_group * SIDRE_group_create_group_bufferify(SIDRE_group * self, const char * path, int Lpath, SIDRE_group * SHC_rv);

void SIDRE_group_destroy_group_name(SIDRE_group * self, const char * path);

void SIDRE_group_destroy_group_name_bufferify(SIDRE_group * self, const char * path, int Lpath);

void SIDRE_group_destroy_group_index(SIDRE_group * self, SIDRE_IndexType idx);

SIDRE_group * SIDRE_group_move_group(SIDRE_group * self, SIDRE_group * grp, SIDRE_group * SHC_rv);

void SIDRE_group_print(const SIDRE_group * self);

bool SIDRE_group_is_equivalent_to(const SIDRE_group * self, const SIDRE_group * other);

void SIDRE_group_save(const SIDRE_group * self, const char * file_path, const char * protocol);

void SIDRE_group_save_bufferify(const SIDRE_group * self, const char * file_path, int Lfile_path, const char * protocol, int Lprotocol);

void SIDRE_group_load_0(SIDRE_group * self, const char * file_path, const char * protocol);

void SIDRE_group_load_0_bufferify(SIDRE_group * self, const char * file_path, int Lfile_path, const char * protocol, int Lprotocol);

void SIDRE_group_load_1(SIDRE_group * self, const char * file_path, const char * protocol, bool preserve_contents);

void SIDRE_group_load_1_bufferify(SIDRE_group * self, const char * file_path, int Lfile_path, const char * protocol, int Lprotocol, bool preserve_contents);

void SIDRE_group_load_external_data(SIDRE_group * self, const char * file_path);

void SIDRE_group_load_external_data_bufferify(SIDRE_group * self, const char * file_path, int Lfile_path);

bool SIDRE_group_rename(SIDRE_group * self, const char * new_name);

bool SIDRE_group_rename_bufferify(SIDRE_group * self, const char * new_name, int Lnew_name);

#ifdef __cplusplus
}
#endif

#endif  // WRAPGROUP_H

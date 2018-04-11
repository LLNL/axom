// wrapView.h
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
 * \file wrapView.h
 * \brief Shroud generated wrapper for View class
 */
// For C users and C++ implementation

#ifndef WRAPVIEW_H
#define WRAPVIEW_H

#include <stddef.h>
#include "sidre/SidreTypes.h"

// splicer begin class.View.CXX_declarations
// splicer end class.View.CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// declaration of shadow types
struct s_SIDRE_buffer;
typedef struct s_SIDRE_buffer SIDRE_buffer;
struct s_SIDRE_group;
typedef struct s_SIDRE_group SIDRE_group;
struct s_SIDRE_view;
typedef struct s_SIDRE_view SIDRE_view;

// splicer begin class.View.C_declarations
// splicer end class.View.C_declarations

SIDRE_IndexType SIDRE_view_get_index(SIDRE_view* self);

const char* SIDRE_view_get_name(const SIDRE_view* self);

void SIDRE_view_get_name_bufferify(const SIDRE_view* self, char* SHF_rv,
                                   int NSHF_rv);

void SIDRE_view_get_path_bufferify(const SIDRE_view* self, char* SHF_rv,
                                   int NSHF_rv);

void SIDRE_view_get_path_name_bufferify(const SIDRE_view* self, char* SHF_rv,
                                        int NSHF_rv);

SIDRE_group* SIDRE_view_get_owning_group(SIDRE_view* self);

bool SIDRE_view_has_buffer(const SIDRE_view* self);

SIDRE_buffer* SIDRE_view_get_buffer(SIDRE_view* self);

bool SIDRE_view_is_external(const SIDRE_view* self);

bool SIDRE_view_is_allocated(SIDRE_view* self);

bool SIDRE_view_is_applied(const SIDRE_view* self);

bool SIDRE_view_is_described(const SIDRE_view* self);

bool SIDRE_view_is_empty(const SIDRE_view* self);

bool SIDRE_view_is_opaque(const SIDRE_view* self);

bool SIDRE_view_is_scalar(const SIDRE_view* self);

bool SIDRE_view_is_string(const SIDRE_view* self);

int SIDRE_view_get_type_id(const SIDRE_view* self);

size_t SIDRE_view_get_total_bytes(const SIDRE_view* self);

size_t SIDRE_view_get_num_elements(const SIDRE_view* self);

size_t SIDRE_view_get_bytes_per_element(const SIDRE_view* self);

size_t SIDRE_view_get_offset(const SIDRE_view* self);

size_t SIDRE_view_get_stride(const SIDRE_view* self);

int SIDRE_view_get_num_dimensions(const SIDRE_view* self);

int SIDRE_view_get_shape(const SIDRE_view* self, int ndims,
                         SIDRE_SidreLength* shape);

void SIDRE_view_allocate_simple(SIDRE_view* self);

void SIDRE_view_allocate_from_type(SIDRE_view* self, int type,
                                   SIDRE_SidreLength num_elems);

void SIDRE_view_reallocate(SIDRE_view* self, SIDRE_SidreLength num_elems);

void SIDRE_view_attach_buffer_only(SIDRE_view* self, SIDRE_buffer* buff);

void SIDRE_view_attach_buffer_type(SIDRE_view* self, int type,
                                   SIDRE_SidreLength num_elems,
                                   SIDRE_buffer* buff);

void SIDRE_view_attach_buffer_shape(SIDRE_view* self, int type, int ndims,
                                    SIDRE_SidreLength* shape,
                                    SIDRE_buffer* buff);

void SIDRE_view_apply_0(SIDRE_view* self);

void SIDRE_view_apply_nelems(SIDRE_view* self, SIDRE_SidreLength num_elems);

void SIDRE_view_apply_nelems_offset(SIDRE_view* self,
                                    SIDRE_SidreLength num_elems,
                                    SIDRE_SidreLength offset);

void SIDRE_view_apply_nelems_offset_stride(SIDRE_view* self,
                                           SIDRE_SidreLength num_elems,
                                           SIDRE_SidreLength offset,
                                           SIDRE_SidreLength stride);

void SIDRE_view_apply_type_nelems(SIDRE_view* self, int type,
                                  SIDRE_SidreLength num_elems);

void SIDRE_view_apply_type_nelems_offset(SIDRE_view* self, int type,
                                         SIDRE_SidreLength num_elems,
                                         SIDRE_SidreLength offset);

void SIDRE_view_apply_type_nelems_offset_stride(SIDRE_view* self, int type,
                                                SIDRE_SidreLength num_elems,
                                                SIDRE_SidreLength offset,
                                                SIDRE_SidreLength stride);

void SIDRE_view_apply_type_shape(SIDRE_view* self, int type, int ndims,
                                 SIDRE_SidreLength* shape);

void SIDRE_view_set_scalar_int(SIDRE_view* self, int value);

void SIDRE_view_set_scalar_long(SIDRE_view* self, long value);

void SIDRE_view_set_scalar_float(SIDRE_view* self, float value);

void SIDRE_view_set_scalar_double(SIDRE_view* self, double value);

void SIDRE_view_set_string(SIDRE_view* self, const char* value);

void SIDRE_view_set_string_bufferify(SIDRE_view* self, const char* value,
                                     int Lvalue);

void SIDRE_view_set_external_data_ptr_only(SIDRE_view* self,
                                           void* external_ptr);

void SIDRE_view_set_external_data_ptr_type(SIDRE_view* self, int type,
                                           SIDRE_SidreLength num_elems,
                                           void* external_ptr);

void SIDRE_view_set_external_data_ptr_shape(SIDRE_view* self, int type,
                                            int ndims, SIDRE_SidreLength* shape,
                                            void* external_ptr);

const char* SIDRE_view_get_string(SIDRE_view* self);

void SIDRE_view_get_string_bufferify(SIDRE_view* self, char* name, int Nname);

int SIDRE_view_get_data_int(SIDRE_view* self);

long SIDRE_view_get_data_long(SIDRE_view* self);

float SIDRE_view_get_data_float(SIDRE_view* self);

double SIDRE_view_get_data_double(SIDRE_view* self);

void* SIDRE_view_get_void_ptr(const SIDRE_view* self);

void SIDRE_view_print(const SIDRE_view* self);

bool SIDRE_view_rename(SIDRE_view* self, const char* new_name);

bool SIDRE_view_rename_bufferify(SIDRE_view* self, const char* new_name,
                                 int Lnew_name);

#ifdef __cplusplus
}
#endif

#endif  // WRAPVIEW_H

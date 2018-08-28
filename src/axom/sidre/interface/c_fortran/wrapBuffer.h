// wrapBuffer.h
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
 * \file wrapBuffer.h
 * \brief Shroud generated wrapper for Buffer class
 */
// For C users and C++ implementation

#ifndef WRAPBUFFER_H
#define WRAPBUFFER_H

#include <stddef.h>
#include "axom/sidre/interface/SidreTypes.h"
#include "typesSidre.h"

// splicer begin class.Buffer.CXX_declarations
// splicer end class.Buffer.CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// splicer begin class.Buffer.C_declarations
// splicer end class.Buffer.C_declarations

SIDRE_IndexType SIDRE_buffer_get_index(const SIDRE_buffer* self);

size_t SIDRE_buffer_get_num_views(const SIDRE_buffer* self);

void* SIDRE_buffer_get_void_ptr(SIDRE_buffer* self);

int SIDRE_buffer_get_type_id(const SIDRE_buffer* self);

size_t SIDRE_buffer_get_num_elements(const SIDRE_buffer* self);

size_t SIDRE_buffer_get_total_bytes(const SIDRE_buffer* self);

size_t SIDRE_buffer_get_bytes_per_element(const SIDRE_buffer* self);

void SIDRE_buffer_describe(SIDRE_buffer* self, int type,
                           SIDRE_IndexType num_elems);

void SIDRE_buffer_allocate_existing(SIDRE_buffer* self);

void SIDRE_buffer_allocate_from_type(SIDRE_buffer* self, int type,
                                     SIDRE_IndexType num_elems);

void SIDRE_buffer_reallocate(SIDRE_buffer* self, SIDRE_IndexType num_elems);

void SIDRE_buffer_print(const SIDRE_buffer* self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPBUFFER_H

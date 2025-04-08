// wrapBuffer.h
// This file is generated by Shroud 0.13.0. Do not edit.
//
// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
/**
 * \file wrapBuffer.h
 * \brief Shroud generated wrapper for Buffer class
 */
// For C users and C++ implementation

#ifndef WRAPBUFFER_H
#define WRAPBUFFER_H

#include "wrapSidre.h"
#include "axom/sidre/interface/SidreTypes.h"
#ifdef __cplusplus
  #include <cstddef>
  #include "axom/sidre/core/SidreTypes.hpp"
#else
  #include <stddef.h>
#endif
#include "typesSidre.h"

// splicer begin class.Buffer.CXX_declarations
// splicer end class.Buffer.CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// splicer begin class.Buffer.C_declarations
// splicer end class.Buffer.C_declarations

SIDRE_IndexType SIDRE_Buffer_get_index(const SIDRE_Buffer* self);

size_t SIDRE_Buffer_get_num_views(const SIDRE_Buffer* self);

void* SIDRE_Buffer_get_void_ptr(SIDRE_Buffer* self);

SIDRE_TypeIDint SIDRE_Buffer_get_type_id(const SIDRE_Buffer* self);

size_t SIDRE_Buffer_get_num_elements(const SIDRE_Buffer* self);

size_t SIDRE_Buffer_get_total_bytes(const SIDRE_Buffer* self);

size_t SIDRE_Buffer_get_bytes_per_element(const SIDRE_Buffer* self);

void SIDRE_Buffer_describe(SIDRE_Buffer* self, SIDRE_TypeID type, SIDRE_IndexType num_elems);

void SIDRE_Buffer_allocate_existing(SIDRE_Buffer* self);

void SIDRE_Buffer_allocate_from_type(SIDRE_Buffer* self, SIDRE_TypeID type, SIDRE_IndexType num_elems);

void SIDRE_Buffer_reallocate(SIDRE_Buffer* self, SIDRE_IndexType num_elems);

void SIDRE_Buffer_print(const SIDRE_Buffer* self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPBUFFER_H

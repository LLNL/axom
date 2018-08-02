// wrapDataStore.h
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
 * \file wrapDataStore.h
 * \brief Shroud generated wrapper for DataStore class
 */
// For C users and C++ implementation

#ifndef WRAPDATASTORE_H
#define WRAPDATASTORE_H

#include <stddef.h>
#include "axom/sidre/interface/SidreTypes.h"

// splicer begin class.DataStore.CXX_declarations
// splicer end class.DataStore.CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// declaration of shadow types
struct s_SIDRE_buffer;
typedef struct s_SIDRE_buffer SIDRE_buffer;
struct s_SIDRE_datastore;
typedef struct s_SIDRE_datastore SIDRE_datastore;
struct s_SIDRE_group;
typedef struct s_SIDRE_group SIDRE_group;

// splicer begin class.DataStore.C_declarations
// splicer end class.DataStore.C_declarations

SIDRE_datastore* SIDRE_datastore_new();

void SIDRE_datastore_delete(SIDRE_datastore* self);

SIDRE_group* SIDRE_datastore_get_root(SIDRE_datastore* self);

size_t SIDRE_datastore_get_num_buffers(const SIDRE_datastore* self);

SIDRE_buffer* SIDRE_datastore_get_buffer(SIDRE_datastore* self,
                                         SIDRE_IndexType idx);

SIDRE_buffer* SIDRE_datastore_create_buffer_empty(SIDRE_datastore* self);

SIDRE_buffer* SIDRE_datastore_create_buffer_from_type(SIDRE_datastore* self,
                                                      int type,
                                                      SIDRE_SidreLength num_elems);

void SIDRE_datastore_destroy_buffer(SIDRE_datastore* self, SIDRE_IndexType id);

void SIDRE_datastore_print(const SIDRE_datastore* self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATASTORE_H

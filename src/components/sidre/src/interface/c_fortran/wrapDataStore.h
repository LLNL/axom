// wrapDataStore.h
// This is generated code, do not edit
//
// Copyright (c) 2015, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
//
// All rights reserved.
//
// This source code cannot be distributed without permission and
// further review from Lawrence Livermore National Laboratory.
//
/**
 * \file wrapDataStore.h
 * \brief Shroud generated wrapper for DataStore class
 */
// For C users and C++ implementation

#ifndef WRAPDATASTORE_H
#define WRAPDATASTORE_H

#include "sidre/SidreTypes.h"
#include "stdlib.h"

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types
struct s_SIDRE_buffer;
typedef struct s_SIDRE_buffer SIDRE_buffer;
struct s_SIDRE_datastore;
typedef struct s_SIDRE_datastore SIDRE_datastore;
struct s_SIDRE_group;
typedef struct s_SIDRE_group SIDRE_group;

// splicer begin class.DataStore.C_definition
// splicer end class.DataStore.C_definition

SIDRE_datastore * SIDRE_datastore_new();

void SIDRE_datastore_delete(SIDRE_datastore * self);

SIDRE_group * SIDRE_datastore_get_root(SIDRE_datastore * self);

SIDRE_buffer * SIDRE_datastore_get_buffer(SIDRE_datastore * self,
                                          SIDRE_IndexType idx);

SIDRE_buffer * SIDRE_datastore_create_buffer_empty(SIDRE_datastore * self);

SIDRE_buffer * SIDRE_datastore_create_buffer_from_type(SIDRE_datastore * self,
                                                       int type,
                                                       SIDRE_SidreLength num_elems);

void SIDRE_datastore_destroy_buffer(SIDRE_datastore * self, SIDRE_IndexType id);

size_t SIDRE_datastore_get_num_buffers(const SIDRE_datastore * self);

void SIDRE_datastore_print(const SIDRE_datastore * self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATASTORE_H

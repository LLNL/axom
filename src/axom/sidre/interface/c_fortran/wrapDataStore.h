// wrapDataStore.h
// This is generated code, do not edit
//
// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
/**
 * \file wrapDataStore.h
 * \brief Shroud generated wrapper for DataStore class
 */
// For C users and C++ implementation

#ifndef WRAPDATASTORE_H
#define WRAPDATASTORE_H

#include <stddef.h>
#include "axom/sidre/interface/SidreTypes.h"
#include "typesSidre.h"

// splicer begin class.DataStore.CXX_declarations
// splicer end class.DataStore.CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// splicer begin class.DataStore.C_declarations
// splicer end class.DataStore.C_declarations

SIDRE_datastore* SIDRE_datastore_new(SIDRE_datastore* SHC_rv);

void SIDRE_datastore_delete(SIDRE_datastore* self);

SIDRE_group* SIDRE_datastore_get_root(SIDRE_datastore* self,
                                      SIDRE_group* SHC_rv);

size_t SIDRE_datastore_get_num_buffers(const SIDRE_datastore* self);

SIDRE_buffer* SIDRE_datastore_get_buffer(SIDRE_datastore* self,
                                         SIDRE_IndexType idx,
                                         SIDRE_buffer* SHC_rv);

SIDRE_buffer* SIDRE_datastore_create_buffer_empty(SIDRE_datastore* self,
                                                  SIDRE_buffer* SHC_rv);

SIDRE_buffer* SIDRE_datastore_create_buffer_from_type(SIDRE_datastore* self,
                                                      int type,
                                                      SIDRE_IndexType num_elems,
                                                      SIDRE_buffer* SHC_rv);

void SIDRE_datastore_destroy_buffer(SIDRE_datastore* self, SIDRE_IndexType id);

bool SIDRE_datastore_generate_blueprint_index(SIDRE_datastore* self,
                                              const char* domain_path,
                                              const char* mesh_name,
                                              const char* index_path,
                                              int num_domains);

bool SIDRE_datastore_generate_blueprint_index_bufferify(SIDRE_datastore* self,
                                                        const char* domain_path,
                                                        int Ldomain_path,
                                                        const char* mesh_name,
                                                        int Lmesh_name,
                                                        const char* index_path,
                                                        int Lindex_path,
                                                        int num_domains);

void SIDRE_datastore_print(const SIDRE_datastore* self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATASTORE_H

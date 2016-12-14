// wrapDataBuffer.h
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
 * \file wrapDataBuffer.h
 * \brief Shroud generated wrapper for DataBuffer class
 */
// For C users and C++ implementation

#ifndef WRAPDATABUFFER_H
#define WRAPDATABUFFER_H

#include "sidre/SidreTypes.h"
#include "stdlib.h"

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types
struct s_SIDRE_databuffer;
typedef struct s_SIDRE_databuffer SIDRE_databuffer;

// splicer begin class.DataBuffer.C_definition
// splicer end class.DataBuffer.C_definition

SIDRE_IndexType SIDRE_databuffer_get_index(const SIDRE_databuffer * self);

size_t SIDRE_databuffer_get_num_views(const SIDRE_databuffer * self);

void SIDRE_databuffer_describe(SIDRE_databuffer * self, int type,
                               SIDRE_SidreLength num_elems);

void SIDRE_databuffer_allocate_existing(SIDRE_databuffer * self);

void SIDRE_databuffer_allocate_from_type(SIDRE_databuffer * self, int type,
                                         SIDRE_SidreLength num_elems);

void SIDRE_databuffer_reallocate(SIDRE_databuffer * self,
                                 SIDRE_SidreLength num_elems);

void * SIDRE_databuffer_get_void_ptr(SIDRE_databuffer * self);

int SIDRE_databuffer_get_type_id(const SIDRE_databuffer * self);

size_t SIDRE_databuffer_get_num_elements(const SIDRE_databuffer * self);

size_t SIDRE_databuffer_get_total_bytes(const SIDRE_databuffer * self);

size_t SIDRE_databuffer_get_bytes_per_element(const SIDRE_databuffer * self);

void SIDRE_databuffer_print(const SIDRE_databuffer * self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATABUFFER_H

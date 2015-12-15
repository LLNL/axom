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

SIDRE_IndexType SIDRE_databuffer_get_index(SIDRE_databuffer * self);

size_t SIDRE_databuffer_get_num_views(SIDRE_databuffer * self);

void SIDRE_databuffer_declare(SIDRE_databuffer * self, int type, SIDRE_SidreLength num_elems);

void SIDRE_databuffer_allocate_existing(SIDRE_databuffer * self);

void SIDRE_databuffer_allocate_from_type(SIDRE_databuffer * self, int type, SIDRE_SidreLength num_elems);

void SIDRE_databuffer_reallocate(SIDRE_databuffer * self, SIDRE_SidreLength num_elems);

void SIDRE_databuffer_set_external_data(SIDRE_databuffer * self, void * external_data);

bool SIDRE_databuffer_is_external(SIDRE_databuffer * self);

void * SIDRE_databuffer_get_void_ptr(SIDRE_databuffer * self);

int SIDRE_databuffer_get_type_id(SIDRE_databuffer * self);

size_t SIDRE_databuffer_get_num_elements(SIDRE_databuffer * self);

size_t SIDRE_databuffer_get_total_bytes(SIDRE_databuffer * self);

void SIDRE_databuffer_print(SIDRE_databuffer * self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATABUFFER_H

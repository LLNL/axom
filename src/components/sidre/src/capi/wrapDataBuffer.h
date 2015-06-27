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
// wrapDataBuffer.h
// For C users and C++ implementation

#ifndef WRAPDATABUFFER_H
#define WRAPDATABUFFER_H

#include "sidre/SidreTypes.h"
#include "stdlib.h"

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types
#ifdef EXAMPLE_WRAPPER_IMPL
typedef void ATK_databuffer;
#else
struct s_ATK_databuffer;
typedef struct s_ATK_databuffer ATK_databuffer;
#endif

// splicer begin class.DataBuffer.C_definition
// splicer end class.DataBuffer.C_definition

ATK_IndexType ATK_databuffer_get_index(ATK_databuffer * self);

size_t ATK_databuffer_get_num_views(ATK_databuffer * self);

void ATK_databuffer_declare(ATK_databuffer * self, int type, ATK_SidreLength len);

void ATK_databuffer_declare_external(ATK_databuffer * self, int type, ATK_SidreLength len);

void ATK_databuffer_set_external_data(ATK_databuffer * self, void * external_data);

void ATK_databuffer_allocate_existing(ATK_databuffer * self);

void ATK_databuffer_allocate_from_type(ATK_databuffer * self, int type, ATK_SidreLength len);

void ATK_databuffer_reallocate(ATK_databuffer * self, int type, ATK_SidreLength len);

bool ATK_databuffer_is_external(ATK_databuffer * self);

void * ATK_databuffer_get_data(ATK_databuffer * self);

size_t ATK_databuffer_get_total_bytes(ATK_databuffer * self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATABUFFER_H

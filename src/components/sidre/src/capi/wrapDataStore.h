//
// Copyright (c) 2015, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
//
// All rights reserved.
//
// This source code cannot be distributed without permission and
// further review from Lawrence Livermore National Laboratory.
//
// wrapDataStore.h
// For C users and C++ implementation

#ifndef WRAPDATASTORE_H
#define WRAPDATASTORE_H

#include "wrapDataBuffer.h"
#include "wrapDataGroup.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef EXAMPLE_WRAPPER_IMPL
typedef void ATK_datastore;
#else
struct s_ATK_datastore;
typedef struct s_ATK_datastore ATK_datastore;
#endif

ATK_datastore * ATK_datastore_new();

void ATK_datastore_delete(ATK_datastore * self);

ATK_databuffer * ATK_datastore_create_buffer(ATK_datastore * self);

ATK_datagroup * ATK_datastore_get_root(ATK_datastore * self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATASTORE_H

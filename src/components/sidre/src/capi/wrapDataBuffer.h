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

#include "sidre/DataTypes.h"

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

ATK_IDType ATK_databuffer_get_uid(ATK_databuffer * self);

ATK_databuffer * ATK_databuffer_declare(ATK_databuffer * self, ATK_TypeEnum type, long len);

ATK_databuffer * ATK_databuffer_allocate(ATK_databuffer * self);

void * ATK_databuffer_get_data(ATK_databuffer * self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATABUFFER_H

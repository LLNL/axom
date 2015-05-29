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

#ifdef __cplusplus
extern "C" {
#endif

#ifdef EXAMPLE_WRAPPER_IMPL
typedef void ATK_databuffer;
#else
struct s_ATK_databuffer;
typedef struct s_ATK_databuffer ATK_databuffer;
#endif

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATABUFFER_H

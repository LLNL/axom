// wrapSLIC.h
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
 * \file wrapSLIC.h
 * \brief Shroud generated wrapper for SLIC library
 */
// For C users and C++ implementation

#ifndef WRAPSLIC_H
#define WRAPSLIC_H

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types

// splicer begin C_definition
// splicer end C_definition

void SLIC_initialize();

bool SLIC_is_initialized();

void SLIC_finalize();

#ifdef __cplusplus
}
#endif

#endif  // WRAPSLIC_H

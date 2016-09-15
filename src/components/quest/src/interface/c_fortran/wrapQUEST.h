// wrapQUEST.h
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
 * \file wrapQUEST.h
 * \brief Shroud generated wrapper for QUEST library
 */
// For C users and C++ implementation

#ifndef WRAPQUEST_H
#define WRAPQUEST_H

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types

// splicer begin C_definition
// splicer end C_definition

void QUEST_initialize(MPI_Fint comm, const char * fileName, bool requiresDistance, int ndims, int maxElements, int maxLevels);

void QUEST_initialize_bufferify(MPI_Fint comm, const char * fileName, int LfileName, bool requiresDistance, int ndims, int maxElements, int maxLevels);

void QUEST_finalize();

double QUEST_distance(double x, double y, double z);

int QUEST_inside(double x, double y, double z);

void QUEST_mesh_min_bounds(double * coords);

void QUEST_mesh_max_bounds(double * coords);

void QUEST_mesh_center_of_mass(double * coords);

#ifdef __cplusplus
}
#endif

#endif  // WRAPQUEST_H

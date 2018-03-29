// wrapQUEST.h
// This is generated code, do not edit
//
// Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
//
// Produced at the Lawrence Livermore National Laboratory
//
// LLNL-CODE-741217
//
// All rights reserved.
//
// This file is part of Axom.
//
// For details about use and distribution, please read axom/LICENSE.
//
/**
 * \file wrapQUEST.h
 * \brief Shroud generated wrapper for QUEST library
 */
// For C users and C++ implementation

#ifndef WRAPQUEST_H
#define WRAPQUEST_H

#include "mpi.h"

// splicer begin CXX_declarations
// splicer end CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// declaration of wrapped types

// splicer begin C_declarations
// splicer end C_declarations

void QUEST_initialize(MPI_Fint comm, const char* fileName,
                      bool requiresDistance, int ndims, int maxElements,
                      int maxLevels);

void QUEST_initialize_bufferify(MPI_Fint comm, const char* fileName,
                                int LfileName, bool requiresDistance, int ndims,
                                int maxElements, int maxLevels);

double QUEST_distance(double x, double y, double z);

int QUEST_inside(double x, double y, double z);

void QUEST_mesh_min_bounds(double* coords);

void QUEST_mesh_max_bounds(double* coords);

void QUEST_mesh_center_of_mass(double* coords);

void QUEST_finalize();

#ifdef __cplusplus
}
#endif

#endif  // WRAPQUEST_H

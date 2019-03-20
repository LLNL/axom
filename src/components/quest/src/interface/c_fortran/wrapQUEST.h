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

#ifdef AXOM_USE_MPI
#include "mpi.h"
#endif

// splicer begin CXX_declarations
// splicer end CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// splicer begin C_declarations
// splicer end C_declarations

#ifdef AXOM_USE_MPI
void QUEST_initialize_mpi(MPI_Fint comm, const char* fileName,
                          bool requiresDistance, int ndims, int maxElements,
                          int maxLevels);
#endif

#ifdef AXOM_USE_MPI
void QUEST_initialize_mpi_bufferify(MPI_Fint comm, const char* fileName,
                                    int LfileName, bool requiresDistance,
                                    int ndims, int maxElements, int maxLevels);
#endif

#ifndef AXOM_USE_MPI
void QUEST_initialize_serial(const char* fileName, bool requiresDistance,
                             int ndims, int maxElements, int maxLevels);
#endif

#ifndef AXOM_USE_MPI
void QUEST_initialize_serial_bufferify(const char* fileName, int LfileName,
                                       bool requiresDistance, int ndims,
                                       int maxElements, int maxLevels);
#endif

double QUEST_distance_0(double x, double y);

double QUEST_distance_1(double x, double y, double z);

int QUEST_inside_0(double x, double y);

int QUEST_inside_1(double x, double y, double z);

void QUEST_mesh_min_bounds(double* coords);

void QUEST_mesh_max_bounds(double* coords);

void QUEST_mesh_center_of_mass(double* coords);

void QUEST_finalize();

#ifdef __cplusplus
}
#endif

#endif  // WRAPQUEST_H

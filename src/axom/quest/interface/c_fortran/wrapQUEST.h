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
#include "typesQUEST.h"

// splicer begin CXX_declarations
// splicer end CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// splicer begin C_declarations
// splicer end C_declarations

#ifdef AXOM_USE_MPI
void QUEST_initialize_mpi(MPI_Fint comm, const char * fileName, bool requiresDistance, int ndims, int maxElements, int maxLevels);
#endif

#ifdef AXOM_USE_MPI
void QUEST_initialize_mpi_bufferify(MPI_Fint comm, const char * fileName, int LfileName, bool requiresDistance, int ndims, int maxElements, int maxLevels);
#endif

#ifndef AXOM_USE_MPI
void QUEST_initialize_serial(const char * fileName, bool requiresDistance, int ndims, int maxElements, int maxLevels);
#endif

#ifndef AXOM_USE_MPI
void QUEST_initialize_serial_bufferify(const char * fileName, int LfileName, bool requiresDistance, int ndims, int maxElements, int maxLevels);
#endif

int QUEST_inside_0(double x, double y);

int QUEST_inside_1(double x, double y, double z);

void QUEST_mesh_min_bounds(double * coords);

void QUEST_mesh_max_bounds(double * coords);

void QUEST_mesh_center_of_mass(double * coords);

void QUEST_finalize();

#ifdef AXOM_USE_MPI
int QUEST_signed_distance_init_mpi(const char * file, MPI_Fint comm);
#endif

#ifdef AXOM_USE_MPI
int QUEST_signed_distance_init_mpi_bufferify(const char * file, int Lfile, MPI_Fint comm);
#endif

#ifndef AXOM_USE_MPI
int QUEST_signed_distance_init_serial(const char * file);
#endif

#ifndef AXOM_USE_MPI
int QUEST_signed_distance_init_serial_bufferify(const char * file, int Lfile);
#endif

bool QUEST_signed_distance_initialized();

void QUEST_signed_distance_get_mesh_bounds(double * lo, double * hi);

void QUEST_signed_distance_set_dimension(int dim);

void QUEST_signed_distance_set_closed_surface(bool status);

void QUEST_signed_distance_set_max_levels(int maxLevels);

void QUEST_signed_distance_set_max_occupancy(int maxOccupancy);

void QUEST_signed_distance_set_verbose(bool status);

void QUEST_signed_distance_use_shared_memory(bool status);

double QUEST_signed_distance_evaluate(double x, double y, double z);

void QUEST_signed_distance_finalize();

#ifdef __cplusplus
}
#endif

#endif  // WRAPQUEST_H

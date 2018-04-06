// wrapQUEST.cpp
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
#include "wrapQUEST.h"
#include <string>
#include "quest/quest.hpp"

// splicer begin CXX_definitions
// splicer end CXX_definitions

extern "C" {

// splicer begin C_definitions
// splicer end C_definitions

#ifdef AXOM_USE_MPI
void QUEST_initialize_mpi(MPI_Fint comm, const char* fileName,
                          bool requiresDistance, int ndims, int maxElements,
                          int maxLevels)
{
// splicer begin function.initialize_mpi
  MPI_Comm SHCXX_comm = MPI_Comm_f2c(comm);
  const std::string SH_fileName(fileName);
  axom::quest::initialize(SHCXX_comm, SH_fileName, requiresDistance, ndims,
                          maxElements, maxLevels);
  return;
// splicer end function.initialize_mpi
}
#endif  // ifdef AXOM_USE_MPI

#ifdef AXOM_USE_MPI
void QUEST_initialize_mpi_bufferify(MPI_Fint comm, const char* fileName,
                                    int LfileName, bool requiresDistance,
                                    int ndims, int maxElements, int maxLevels)
{
// splicer begin function.initialize_mpi_bufferify
  MPI_Comm SHCXX_comm = MPI_Comm_f2c(comm);
  const std::string SH_fileName(fileName, LfileName);
  axom::quest::initialize(SHCXX_comm, SH_fileName, requiresDistance, ndims,
                          maxElements, maxLevels);
  return;
// splicer end function.initialize_mpi_bufferify
}
#endif  // ifdef AXOM_USE_MPI

#ifndef AXOM_USE_MPI
void QUEST_initialize_serial(const char* fileName, bool requiresDistance,
                             int ndims, int maxElements, int maxLevels)
{
// splicer begin function.initialize_serial
  const std::string SH_fileName(fileName);
  axom::quest::initialize(SH_fileName, requiresDistance, ndims, maxElements,
                          maxLevels);
  return;
// splicer end function.initialize_serial
}
#endif  // ifndef AXOM_USE_MPI

#ifndef AXOM_USE_MPI
void QUEST_initialize_serial_bufferify(const char* fileName, int LfileName,
                                       bool requiresDistance, int ndims,
                                       int maxElements, int maxLevels)
{
// splicer begin function.initialize_serial_bufferify
  const std::string SH_fileName(fileName, LfileName);
  axom::quest::initialize(SH_fileName, requiresDistance, ndims, maxElements,
                          maxLevels);
  return;
// splicer end function.initialize_serial_bufferify
}
#endif  // ifndef AXOM_USE_MPI

double QUEST_distance_0(double x, double y)
{
// splicer begin function.distance_0
  double SHC_rv = axom::quest::distance(x, y);
  return SHC_rv;
// splicer end function.distance_0
}

double QUEST_distance_1(double x, double y, double z)
{
// splicer begin function.distance_1
  double SHC_rv = axom::quest::distance(x, y, z);
  return SHC_rv;
// splicer end function.distance_1
}

int QUEST_inside_0(double x, double y)
{
// splicer begin function.inside_0
  int SHC_rv = axom::quest::inside(x, y);
  return SHC_rv;
// splicer end function.inside_0
}

int QUEST_inside_1(double x, double y, double z)
{
// splicer begin function.inside_1
  int SHC_rv = axom::quest::inside(x, y, z);
  return SHC_rv;
// splicer end function.inside_1
}

void QUEST_mesh_min_bounds(double* coords)
{
// splicer begin function.mesh_min_bounds
  axom::quest::mesh_min_bounds(coords);
  return;
// splicer end function.mesh_min_bounds
}

void QUEST_mesh_max_bounds(double* coords)
{
// splicer begin function.mesh_max_bounds
  axom::quest::mesh_max_bounds(coords);
  return;
// splicer end function.mesh_max_bounds
}

void QUEST_mesh_center_of_mass(double* coords)
{
// splicer begin function.mesh_center_of_mass
  axom::quest::mesh_center_of_mass(coords);
  return;
// splicer end function.mesh_center_of_mass
}

void QUEST_finalize()
{
// splicer begin function.finalize
  axom::quest::finalize();
  return;
// splicer end function.finalize
}

}  // extern "C"

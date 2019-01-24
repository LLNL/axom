// wrapQUEST.cpp
// This is generated code, do not edit
//
// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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
#include <stdlib.h>
#include <string>
#include "axom/quest/interface/inout.hpp"
#include "axom/quest/interface/signed_distance.hpp"
#include "typesQUEST.h"

// splicer begin CXX_definitions
// splicer end CXX_definitions

extern "C" {

// splicer begin C_definitions
// splicer end C_definitions

#ifdef AXOM_USE_MPI
int QUEST_inout_init_mpi(const char* fileName, MPI_Fint comm)
{
// splicer begin function.inout_init_mpi
  const std::string SH_fileName(fileName);
  MPI_Comm SHCXX_comm = MPI_Comm_f2c(comm);
  int SHC_rv = axom::quest::inout_init(SH_fileName, SHCXX_comm);
  return SHC_rv;
// splicer end function.inout_init_mpi
}
#endif  // ifdef AXOM_USE_MPI

#ifdef AXOM_USE_MPI
int QUEST_inout_init_mpi_bufferify(const char* fileName, int LfileName,
                                   MPI_Fint comm)
{
// splicer begin function.inout_init_mpi_bufferify
  const std::string SH_fileName(fileName, LfileName);
  MPI_Comm SHCXX_comm = MPI_Comm_f2c(comm);
  int SHC_rv = axom::quest::inout_init(SH_fileName, SHCXX_comm);
  return SHC_rv;
// splicer end function.inout_init_mpi_bufferify
}
#endif  // ifdef AXOM_USE_MPI

#ifndef AXOM_USE_MPI
int QUEST_inout_init_serial(const char* fileName)
{
// splicer begin function.inout_init_serial
  const std::string SH_fileName(fileName);
  int SHC_rv = axom::quest::inout_init(SH_fileName);
  return SHC_rv;
// splicer end function.inout_init_serial
}
#endif  // ifndef AXOM_USE_MPI

#ifndef AXOM_USE_MPI
int QUEST_inout_init_serial_bufferify(const char* fileName, int LfileName)
{
// splicer begin function.inout_init_serial_bufferify
  const std::string SH_fileName(fileName, LfileName);
  int SHC_rv = axom::quest::inout_init(SH_fileName);
  return SHC_rv;
// splicer end function.inout_init_serial_bufferify
}
#endif  // ifndef AXOM_USE_MPI

bool QUEST_inout_initialized()
{
// splicer begin function.inout_initialized
  bool SHC_rv = axom::quest::inout_initialized();
  return SHC_rv;
// splicer end function.inout_initialized
}

int QUEST_inout_set_verbose(bool verbosity)
{
// splicer begin function.inout_set_verbose
  int SHC_rv = axom::quest::inout_set_verbose(verbosity);
  return SHC_rv;
// splicer end function.inout_set_verbose
}

bool QUEST_inout_evaluate_0(double x, double y)
{
// splicer begin function.inout_evaluate_0
  bool SHC_rv = axom::quest::inout_evaluate(x, y);
  return SHC_rv;
// splicer end function.inout_evaluate_0
}

bool QUEST_inout_evaluate_1(double x, double y, double z)
{
// splicer begin function.inout_evaluate_1
  bool SHC_rv = axom::quest::inout_evaluate(x, y, z);
  return SHC_rv;
// splicer end function.inout_evaluate_1
}

int QUEST_inout_mesh_min_bounds(double* coords)
{
// splicer begin function.inout_mesh_min_bounds
  int SHC_rv = axom::quest::inout_mesh_min_bounds(coords);
  return SHC_rv;
// splicer end function.inout_mesh_min_bounds
}

int QUEST_inout_mesh_max_bounds(double* coords)
{
// splicer begin function.inout_mesh_max_bounds
  int SHC_rv = axom::quest::inout_mesh_max_bounds(coords);
  return SHC_rv;
// splicer end function.inout_mesh_max_bounds
}

int QUEST_inout_mesh_center_of_mass(double* coords)
{
// splicer begin function.inout_mesh_center_of_mass
  int SHC_rv = axom::quest::inout_mesh_center_of_mass(coords);
  return SHC_rv;
// splicer end function.inout_mesh_center_of_mass
}

int QUEST_inout_get_dimension()
{
// splicer begin function.inout_get_dimension
  int SHC_rv = axom::quest::inout_get_dimension();
  return SHC_rv;
// splicer end function.inout_get_dimension
}

int QUEST_inout_finalize()
{
// splicer begin function.inout_finalize
  int SHC_rv = axom::quest::inout_finalize();
  return SHC_rv;
// splicer end function.inout_finalize
}

#ifdef AXOM_USE_MPI
int QUEST_signed_distance_init_mpi(const char* file, MPI_Fint comm)
{
// splicer begin function.signed_distance_init_mpi
  const std::string SH_file(file);
  MPI_Comm SHCXX_comm = MPI_Comm_f2c(comm);
  int SHC_rv = axom::quest::signed_distance_init(SH_file, SHCXX_comm);
  return SHC_rv;
// splicer end function.signed_distance_init_mpi
}
#endif  // ifdef AXOM_USE_MPI

#ifdef AXOM_USE_MPI
int QUEST_signed_distance_init_mpi_bufferify(const char* file, int Lfile,
                                             MPI_Fint comm)
{
// splicer begin function.signed_distance_init_mpi_bufferify
  const std::string SH_file(file, Lfile);
  MPI_Comm SHCXX_comm = MPI_Comm_f2c(comm);
  int SHC_rv = axom::quest::signed_distance_init(SH_file, SHCXX_comm);
  return SHC_rv;
// splicer end function.signed_distance_init_mpi_bufferify
}
#endif  // ifdef AXOM_USE_MPI

#ifndef AXOM_USE_MPI
int QUEST_signed_distance_init_serial(const char* file)
{
// splicer begin function.signed_distance_init_serial
  const std::string SH_file(file);
  int SHC_rv = axom::quest::signed_distance_init(SH_file);
  return SHC_rv;
// splicer end function.signed_distance_init_serial
}
#endif  // ifndef AXOM_USE_MPI

#ifndef AXOM_USE_MPI
int QUEST_signed_distance_init_serial_bufferify(const char* file, int Lfile)
{
// splicer begin function.signed_distance_init_serial_bufferify
  const std::string SH_file(file, Lfile);
  int SHC_rv = axom::quest::signed_distance_init(SH_file);
  return SHC_rv;
// splicer end function.signed_distance_init_serial_bufferify
}
#endif  // ifndef AXOM_USE_MPI

bool QUEST_signed_distance_initialized()
{
// splicer begin function.signed_distance_initialized
  bool SHC_rv = axom::quest::signed_distance_initialized();
  return SHC_rv;
// splicer end function.signed_distance_initialized
}

void QUEST_signed_distance_get_mesh_bounds(double* lo, double* hi)
{
// splicer begin function.signed_distance_get_mesh_bounds
  axom::quest::signed_distance_get_mesh_bounds(lo, hi);
  return;
// splicer end function.signed_distance_get_mesh_bounds
}

void QUEST_signed_distance_set_dimension(int dim)
{
// splicer begin function.signed_distance_set_dimension
  axom::quest::signed_distance_set_dimension(dim);
  return;
// splicer end function.signed_distance_set_dimension
}

void QUEST_signed_distance_set_closed_surface(bool status)
{
// splicer begin function.signed_distance_set_closed_surface
  axom::quest::signed_distance_set_closed_surface(status);
  return;
// splicer end function.signed_distance_set_closed_surface
}

void QUEST_signed_distance_set_max_levels(int maxLevels)
{
// splicer begin function.signed_distance_set_max_levels
  axom::quest::signed_distance_set_max_levels(maxLevels);
  return;
// splicer end function.signed_distance_set_max_levels
}

void QUEST_signed_distance_set_max_occupancy(int maxOccupancy)
{
// splicer begin function.signed_distance_set_max_occupancy
  axom::quest::signed_distance_set_max_occupancy(maxOccupancy);
  return;
// splicer end function.signed_distance_set_max_occupancy
}

void QUEST_signed_distance_set_verbose(bool status)
{
// splicer begin function.signed_distance_set_verbose
  axom::quest::signed_distance_set_verbose(status);
  return;
// splicer end function.signed_distance_set_verbose
}

void QUEST_signed_distance_use_shared_memory(bool status)
{
// splicer begin function.signed_distance_use_shared_memory
  axom::quest::signed_distance_use_shared_memory(status);
  return;
// splicer end function.signed_distance_use_shared_memory
}

double QUEST_signed_distance_evaluate(double x, double y, double z)
{
// splicer begin function.signed_distance_evaluate
  double SHC_rv = axom::quest::signed_distance_evaluate(x, y, z);
  return SHC_rv;
// splicer end function.signed_distance_evaluate
}

void QUEST_signed_distance_finalize()
{
// splicer begin function.signed_distance_finalize
  axom::quest::signed_distance_finalize();
  return;
// splicer end function.signed_distance_finalize
}

// Release C++ allocated memory.
void QUEST_SHROUD_memory_destructor(QUE_SHROUD_capsule_data* cap)
{
  cap->addr = NULL;
  cap->idtor = 0;    // avoid deleting again
}

}  // extern "C"

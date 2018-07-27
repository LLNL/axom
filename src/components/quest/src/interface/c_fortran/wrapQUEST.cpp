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
#include "quest/signed_distance.hpp"

// splicer begin CXX_definitions
// splicer end CXX_definitions

extern "C" {

// splicer begin C_definitions
// splicer end C_definitions

#ifdef AXOM_USE_MPI
void QUEST_initialize_mpi(MPI_Fint comm, const char * fileName, bool requiresDistance, int ndims, int maxElements, int maxLevels)
{
// splicer begin function.initialize_mpi
    MPI_Comm SHCXX_comm = MPI_Comm_f2c(comm);
    const std::string SH_fileName(fileName);
    axom::quest::initialize(SHCXX_comm, SH_fileName, requiresDistance, ndims, maxElements, maxLevels);
    return;
// splicer end function.initialize_mpi
}
#endif  // ifdef AXOM_USE_MPI

#ifdef AXOM_USE_MPI
void QUEST_initialize_mpi_bufferify(MPI_Fint comm, const char * fileName, int LfileName, bool requiresDistance, int ndims, int maxElements, int maxLevels)
{
// splicer begin function.initialize_mpi_bufferify
    MPI_Comm SHCXX_comm = MPI_Comm_f2c(comm);
    const std::string SH_fileName(fileName, LfileName);
    axom::quest::initialize(SHCXX_comm, SH_fileName, requiresDistance, ndims, maxElements, maxLevels);
    return;
// splicer end function.initialize_mpi_bufferify
}
#endif  // ifdef AXOM_USE_MPI

#ifndef AXOM_USE_MPI
void QUEST_initialize_serial(const char * fileName, bool requiresDistance, int ndims, int maxElements, int maxLevels)
{
// splicer begin function.initialize_serial
    const std::string SH_fileName(fileName);
    axom::quest::initialize(SH_fileName, requiresDistance, ndims, maxElements, maxLevels);
    return;
// splicer end function.initialize_serial
}
#endif  // ifndef AXOM_USE_MPI

#ifndef AXOM_USE_MPI
void QUEST_initialize_serial_bufferify(const char * fileName, int LfileName, bool requiresDistance, int ndims, int maxElements, int maxLevels)
{
// splicer begin function.initialize_serial_bufferify
    const std::string SH_fileName(fileName, LfileName);
    axom::quest::initialize(SH_fileName, requiresDistance, ndims, maxElements, maxLevels);
    return;
// splicer end function.initialize_serial_bufferify
}
#endif  // ifndef AXOM_USE_MPI

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

void QUEST_mesh_min_bounds(double * coords)
{
// splicer begin function.mesh_min_bounds
    axom::quest::mesh_min_bounds(coords);
    return;
// splicer end function.mesh_min_bounds
}

void QUEST_mesh_max_bounds(double * coords)
{
// splicer begin function.mesh_max_bounds
    axom::quest::mesh_max_bounds(coords);
    return;
// splicer end function.mesh_max_bounds
}

void QUEST_mesh_center_of_mass(double * coords)
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

#ifdef AXOM_USE_MPI
void QUEST_signed_distance_init_mpi(const char * file, MPI_Fint comm)
{
// splicer begin function.signed_distance_init_mpi
    const std::string SH_file(file);
    MPI_Comm SHCXX_comm = MPI_Comm_f2c(comm);
    axom::quest::signed_distance_init(SH_file, SHCXX_comm);
    return;
// splicer end function.signed_distance_init_mpi
}
#endif  // ifdef AXOM_USE_MPI

#ifdef AXOM_USE_MPI
void QUEST_signed_distance_init_mpi_bufferify(const char * file, int Lfile, MPI_Fint comm)
{
// splicer begin function.signed_distance_init_mpi_bufferify
    const std::string SH_file(file, Lfile);
    MPI_Comm SHCXX_comm = MPI_Comm_f2c(comm);
    axom::quest::signed_distance_init(SH_file, SHCXX_comm);
    return;
// splicer end function.signed_distance_init_mpi_bufferify
}
#endif  // ifdef AXOM_USE_MPI

#ifndef AXOM_USE_MPI
void QUEST_signed_distance_init_serial(const char * file)
{
// splicer begin function.signed_distance_init_serial
    const std::string SH_file(file);
    axom::quest::signed_distance_init(SH_file);
    return;
// splicer end function.signed_distance_init_serial
}
#endif  // ifndef AXOM_USE_MPI

#ifndef AXOM_USE_MPI
void QUEST_signed_distance_init_serial_bufferify(const char * file, int Lfile)
{
// splicer begin function.signed_distance_init_serial_bufferify
    const std::string SH_file(file, Lfile);
    axom::quest::signed_distance_init(SH_file);
    return;
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

void QUEST_signed_distance_get_mesh_bounds(double * lo, double * hi)
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

}  // extern "C"

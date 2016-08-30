// wrapQUEST.cpp
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
// wrapQUEST.cpp
#include "wrapQUEST.h"
#include <string>
#include "quest/quest.hpp"

extern "C" {
namespace quest {

void QUEST_initialize(MPI_Fint comm, const char * fileName, bool requiresDistance, int ndims, int maxElements, int maxLevels)
{
// splicer begin function.initialize
const std::string SH_fileName(fileName);
initialize(MPI_Comm_f2c(comm), SH_fileName, requiresDistance, ndims, maxElements, maxLevels);
return;
// splicer end function.initialize
}

void QUEST_initialize_bufferify(MPI_Fint comm, const char * fileName, int LfileName, bool requiresDistance, int ndims, int maxElements, int maxLevels)
{
// splicer begin function.initialize_bufferify
const std::string SH_fileName(fileName, LfileName);
initialize(MPI_Comm_f2c(comm), SH_fileName, requiresDistance, ndims, maxElements, maxLevels);
return;
// splicer end function.initialize_bufferify
}

void QUEST_finalize()
{
// splicer begin function.finalize
finalize();
return;
// splicer end function.finalize
}

double QUEST_distance(double x, double y, double z)
{
// splicer begin function.distance
double rv = distance(x, y, z);
return rv;
// splicer end function.distance
}

int QUEST_inside(double x, double y, double z)
{
// splicer begin function.inside
int rv = inside(x, y, z);
return rv;
// splicer end function.inside
}

double QUEST_mesh_min_x()
{
// splicer begin function.mesh_min_x
double rv = mesh_min_x();
return rv;
// splicer end function.mesh_min_x
}

double QUEST_mesh_min_y()
{
// splicer begin function.mesh_min_y
double rv = mesh_min_y();
return rv;
// splicer end function.mesh_min_y
}

double QUEST_mesh_min_z()
{
// splicer begin function.mesh_min_z
double rv = mesh_min_z();
return rv;
// splicer end function.mesh_min_z
}

double QUEST_mesh_max_x()
{
// splicer begin function.mesh_max_x
double rv = mesh_max_x();
return rv;
// splicer end function.mesh_max_x
}

double QUEST_mesh_max_y()
{
// splicer begin function.mesh_max_y
double rv = mesh_max_y();
return rv;
// splicer end function.mesh_max_y
}

double QUEST_mesh_max_z()
{
// splicer begin function.mesh_max_z
double rv = mesh_max_z();
return rv;
// splicer end function.mesh_max_z
}

double QUEST_mesh_center_of_mass_x()
{
// splicer begin function.mesh_center_of_mass_x
double rv = mesh_center_of_mass_x();
return rv;
// splicer end function.mesh_center_of_mass_x
}

double QUEST_mesh_center_of_mass_y()
{
// splicer begin function.mesh_center_of_mass_y
double rv = mesh_center_of_mass_y();
return rv;
// splicer end function.mesh_center_of_mass_y
}

double QUEST_mesh_center_of_mass_z()
{
// splicer begin function.mesh_center_of_mass_z
double rv = mesh_center_of_mass_z();
return rv;
// splicer end function.mesh_center_of_mass_z
}

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace quest
}  // extern "C"

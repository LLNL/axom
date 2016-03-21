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

void QUEST_initialize(int mpicomm, const char * fileName, int ndims, int maxElements, int maxLevels)
{
// splicer begin function.initialize
std::string SH_fileName(fileName);
initialize(mpicomm, SH_fileName, ndims, maxElements, maxLevels);
return;
// splicer end function.initialize
}

void QUEST_initialize_bufferify(int mpicomm, const char * fileName, int LfileName, int ndims, int maxElements, int maxLevels)
{
// splicer begin function.initialize_bufferify
std::string SH_fileName(fileName, LfileName);
initialize(mpicomm, SH_fileName, ndims, maxElements, maxLevels);
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

// splicer begin additional_functions
// splicer end additional_functions

}  // namespace quest
}  // extern "C"

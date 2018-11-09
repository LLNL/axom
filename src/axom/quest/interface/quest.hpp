/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef QUEST_HPP_
#define QUEST_HPP_

// C/C++ includes
#include <string>

#include "axom/config.hpp"
#include "axom/mint.hpp"

#ifdef AXOM_USE_MPI
#include "mpi.h"
#endif


namespace axom
{
namespace quest
{

/*!
 * \brief Initializes quest.
 * \param [in] comm communicator to use (when running in parallel)
 * \param [in] fileName the name of the file to read in the surface.
 * \param [in] input_mesh pointer to already-read-in mint::Mesh
 * \param [in] requiresDistance flag to determine which structure to build.
 * \param [in] ndims the surface's spatial dimension
 * \param [in] maxElements max elements per bucket.
 * \param [in] maxLevel max levels of subdivision.
 * \param [in] deleteMesh if true, finalize() will delete the mesh;
 *       if false, the user is responsible for deleting the mesh
 * \note If requiresDistance is true, we will build a structure that supports
 *       signed distance queries in addition to containment queries.
 *       Otherwise, we build a structure that only supports containment queries.
 */
#ifdef AXOM_USE_MPI
void initialize( MPI_Comm comm, const std::string& fileName,
                 bool requiresDistance, int ndims, int maxElements,
                 int maxLevels );
void initialize( MPI_Comm comm, mint::Mesh*& input_mesh,
                 bool requiresDistance, int ndims, int maxElements,
                 int maxLevels, bool deleteMesh = false );
#else
void initialize( const std::string& fileName,
                 bool requiresDistance, int ndims, int maxElements,
                 int maxLevels );
void initialize( mint::Mesh*& input_mesh,
                 bool requiresDistance, int ndims, int maxElements,
                 int maxLevels, bool deleteMesh = false );
#endif

/*!
 * \brief Checks if the given point is inside or outside.
 * \param [in] x the x-coordinate of the point in query.
 * \param [in] y the y-coordinate of the point in query.
 * \param [in] z the z-coordinate of the point in query.
 * \return return_value a return value of 1 indicates inside, 0 otherwise.
 */
int inside(double x, double y, double z=0.0);

/*!
 * \brief Checks if the given set of points are inside or outside.
 * \param [in]  xyz user-supplied array of coordinates stored
 * \param [out] in user-supplied array where to store result for each point.
 * \param [in]  npoints total number of points.
 * \pre xyz != nullptr
 * \pre  in != nullptr
 */
void inside( const double* xyz, int* in, int npoints);



/*!
 * \brief Gets coordinates of the minimum corner of the mesh's bounding box
 * \param [out] coords user-supplied array to store the coordinates
 * \pre coords must be preallocated with room for at least three doubles.
 */
void mesh_min_bounds(double* coords);

/*!
 * \brief Gets coordinates of the maximum corner of the mesh's bounding box
 * \param [out] coords user-supplied array to store the coordinates
 * \pre coords must be preallocated with room for at least three doubles.
 */
void mesh_max_bounds(double* coords);

/*!
 * \brief Gets coordinates of the mesh's center of mass
 * \param [out] coords user-supplied array to store the coordinates
 * \pre coords must be preallocated with room for at least three doubles.
 * \note The center of mass is computed as the average vertex position
 */
void mesh_center_of_mass(double* coords);



/*!
 * \brief Finalize quest.
 */
void finalize();

} // end namespace quest
} // end namespace axom

#endif /* QUEST_HPP_ */

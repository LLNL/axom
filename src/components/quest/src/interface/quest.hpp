/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#ifndef QUEST_HPP_
#define QUEST_HPP_

// C/C++ includes
#include <string>

#ifdef USE_MPI
#include "mpi.h"
#endif


namespace quest
{

/*!
 *******************************************************************************
 * \brief Initializes quest.
 * \param [in] comm communicator to use (if running in parallel)
 * \param [in] fileName the name of the file to read in the surface.
 * \param [in] maxElements max elements per bucket.
 * \param [in] maxLevel max levels of subdivision.
 *******************************************************************************
 */
#ifdef USE_MPI
void initialize( MPI_Comm comm, const std::string& fileName,
                 int ndims, int maxElements, int maxLevels );
#else
void initialize( const std::string& fileName, int ndims,
                 int maxElements, int maxLevels );
#endif

/*!
 *******************************************************************************
 * \brief Computes the signed distance of the given point to the surface.
 * \param [in] x the x-coordinate of the point in query.
 * \param [in] y the y-coordinate of the point in query.
 * \param [in] z the z-coordinate of the point in query.
 * \return dist the signed distance function.
 *******************************************************************************
 */
double distance(double x, double y, double z=0.0);

/*!
 *******************************************************************************
 * \brief Computes the signed distance for a set of points to the surface.
 * \param [in]  xyz user-supplied array of coordinates.
 * \param [out] dist user-supplied array where to store the signed distance.
 * \param [in]  npoints total number of points.
 * \pre xyz  != ATK_NULLPTR
 * \pre dist != ATK_NULLPTR
 * \pre npoints >= 0
 *******************************************************************************
 */
void distance( const double* xyz, double* dist, int npoints);

/*!
 *******************************************************************************
 * \brief Checks if the given point is inside or outside.
 * \param [in] x the x-coordinate of the point in query.
 * \param [in] y the y-coordinate of the point in query.
 * \param [in] z the z-coordinate of the point in query.
 * \return return_value a return value of -1 indicates inside, +1 otherwise.
 *******************************************************************************
 */
int inside(double x, double y, double z=0.0);

/*!
 *******************************************************************************
 * \brief Checks if the given set of points are inside or outside.
 * \param [in]  xyz user-supplied array of coordinates stored
 * \param [out] in user-supplied array where to store result for each point.
 * \param [in]  npoints total number of points.
 * \pre xyz != ATK_NULLPTR
 * \pre  in != ATK_NULLPTR
 *******************************************************************************
 */
void inside( const double* xyz, int* in, int npoints);

/*!
 *******************************************************************************
 * \brief Finalize quest.
 *******************************************************************************
 */
void finalize();

} /* end namespace quest */

#endif /* QUEST_HPP_ */

/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

#ifndef QUEST_INOUT_INTERFACE_HPP_
#define QUEST_INOUT_INTERFACE_HPP_

// Axom includes
#include "axom/config.hpp"  // for compile-time configuration options

// Quest includes
#include "axom/quest/interface/internal/mpicomm_wrapper.hpp" // MPI_COMM_SELF

// C/C++ includes
#include <string>   // for std::string

namespace axom
{

// Forward Mint declarations
namespace mint
{
class Mesh;
}

namespace quest
{
constexpr int QUEST_INOUT_SUCCESS = 0;
constexpr int QUEST_INOUT_FAILED = -1;

/// \name Functions related to initializing and finalizing the InOut query
/// @{

int inout_init(const std::string& file, MPI_Comm comm = MPI_COMM_SELF);
int inout_init(mint::Mesh*& mesh, MPI_Comm comm = MPI_COMM_SELF);

int inout_finalize();

/// @}


/*!
 * \brief Predicate to test whether the inout query has been initialized
 *
 * \return True if the inout query has been initialized, false otherwise.
 */
bool inout_initialized();



/*!
 * \brief Tests if the point (\a x, \a y, \a z) is inside the contained volume
 *
 * \param x The x-coordinate of the query point
 * \param y The y-coordinate of the query point
 * \param z The z-coordinate of the query point
 * \return True if the point is within the contained volume, false otherwise.
 */
bool inout_inside(double x, double y, double z=0.);
//void inout_evaluate(
//    const double* x,
//    const double* y,
//    const double* z,
//    int npoints,
//    bool* res);

/*!
 * \brief Returns the lower coordinates of a bounding box containing the mesh
 *
 * \param coords A buffer for the coordinates
 * \pre \a coords != nullptr and has sufficient storage for the coordinates
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 */
int inout_mesh_min_bounds(double* coords);

/*!
 * \brief Returns the upper coordinates of a bounding box containing the mesh
 *
 * \param coords A buffer for the coordinates
 * \pre \a coords != nullptr and has sufficient storage for the coordinates
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 */
int inout_mesh_max_bounds(double* coords);

/*!
 * \brief Returns the center of mass of the mesh
 *
 * The function computes a discrete center of mass defined by the average of the
 * mesh coordinates rather than a continuous center of mass defined by the mesh
 * faces.
 *
 * \param coords A buffer for the coordinates
 * \pre \a coords != nullptr and has sufficient storage for the coordinates
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 */
int inout_mesh_center_of_mass(double* coords);

/// \name Options for InOut query
/// \note These must be called before initializing the query
/// @{

/*!
 * \brief Sets the logging verbosity
 *
 * \param verbosity True for more verbose, false
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 */
int inout_set_verbose(bool verbosity);

/*!
 * \brief Sets the spatial dimension of the mesh
 *
 * \param dimension The spatial dimension
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 */
int inout_set_dimension(int dimension);

/// @}



} // end namespace quest
} // end namespace axom

#endif /* QUEST_INOUT_INTERFACE_HPP_ */

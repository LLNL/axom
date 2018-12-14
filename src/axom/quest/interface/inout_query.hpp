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

/*!
 * \brief Initialize the inout query from a mesh file
 *
 * \param [in] file Path to an STL file containing the surface mesh
 * \param comm The MPI communicator
 * \return rc Return code indicating success of the operation
 *
 * \pre inout_initialized() == false
 * \post inout_initialized() == true when rc is QUEST_INOUT_SUCCESS
 */
int inout_init(const std::string& file, MPI_Comm comm = MPI_COMM_SELF);

/*!
 * \brief Initialize the inout query using a pre-loaded mesh
 *
 * \param [inout] mesh The input mesh. Note: This pointer will be updated
 * during this invocation
 * \param comm The MPI communicator
 * \return rc Return code indicating success of the operation
 *
 * \pre inout_initialized() == false
 * \post inout_initialized() == true when rc is QUEST_INOUT_SUCCESS
 * \note The underlying data structure modifies the input mesh (e.g.
 * by welding vertices) and updates the \a mesh pointer. It is the user's
 * responsibility to update any other pointers to this same mesh.
 */
int inout_init(mint::Mesh*& mesh, MPI_Comm comm = MPI_COMM_SELF);

/*!
 * \brief Finalizes the inout query
 */
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

/*!
 * \brief Tests an array of points for containment
 *
 * Upon successful completion, entries in array \a res
 * will be \a true for points that are inside the enclosed volume
 * and \a false otherwise.
 *
 * \param [in] x Array of x-coordinates for the query points
 * \param [in] y Array of y-coordinates for the query points
 * \param [in] z Array of z-coordinates for the query points
 * \param [in] npoints The number of points to test
 * \param [out] res An array of results. Each entry has value \a 1
 * if the corresponding point is inside and \a 0 otherwise.
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 *
 *  \pre When \a npoints is greater than zero, arrays \a x, \a y, \a z
 *  and \a res are not nullptr and contain sufficient data/space for
 *  \a npoints points.
 */
int inout_inside(const double* x,const double* y,const double* z,
                 int npoints, int* res);

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
 * \warning The quest inout_query is only currently defined
 * for 3D meshes. This function is in anticipation of
 * support for 2D meshes.
 */
int inout_set_dimension(int dimension);

/// @}



} // end namespace quest
} // end namespace axom

#endif /* QUEST_INOUT_INTERFACE_HPP_ */

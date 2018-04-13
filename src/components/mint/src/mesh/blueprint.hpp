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

#ifndef MINT_MESH_BLUEPRINT_HPP_
#define MINT_MESH_BLUEPRINT_HPP_

// C/C++ includes
#include <string> // for std::string

/*!
 * \file
 *
 * \brief Provides convenience methods for navigating and querying a
 *  Sidre group hierarchy according to the conventions outlined in the
 *  <a href="http://llnl-conduit.readthedocs.io/en/latest/"> mesh blueprint</a>.
 *
 */

namespace axom
{

// Sidre Forward Declarations
namespace sidre
{
class Group;
}

namespace mint
{

namespace blueprint
{

/*!
 * \brief Checks if the given root group conforms to the mesh blueprint.
 *
 * \param [in] group pointer to the blueprint root group.
 * \return status true if the group conforms to the blueprint, else, false.
 *
 * \note A root group conforms to the computational mesh blueprint if it
 *  satisfies the following conditions:
 *  (a) it has a "coordsets" child group
 *  (b) it has a "topologies" child group
 *  (c) it has a "fields" child group
 *
 * \note This method emits warnings for each missing group it encounters.
 *
 * \pre group != AXOM_NULLPTR
 */
bool validRootGroup( const sidre::Group* group );

/*!
 * \brief Checks if the given topology group conforms to the blueprint.
 *
 * \param [in] topo pointer to a topology group.
 * \return status true if the group conforms to the blueprint, else, false.
 *
 * \note This routine performs cursory mesh-agnostic checks, i.e., checks that
 *  the necessary child views exists etc., but, it does not validate the content
 *  of those view at this level.
 *
 * \note A valid topology group conforms to the computational mesh blueprint if
 *  it satisfies the following conditions:
 *  (a) it has a "type" view
 *  (b) it has a "coordset" view
 *  (c) the "type" and "coordset" views are both string views.
 *
 * \pre topo != AXOM_NULLPTR
 */
bool validTopologyGroup( const sidre::Group* topo );

/*!
 * \brief Checks if the given coordset group conforms to the blueprint.
 *
 * \param [in] coordset pointer to a coordset group.
 * \return status true if the group conforms to the blueprint, else, false.
 *
 * \note This routine performs cursory mesh-agnostic checks, i.e., checks that
 *  the necessary child views exists etc., but, it does not validate the content
 *  of those view at this level.
 *
 * \note A valid coordset group conforms to the computational mesh blueprint
 *  if  it satisfies the following conditions:
 *  (a) it has a "type" view
 *  (b) the "type" view is a string view
 *
 * \pre coordset != AXOM_NULLPTR
 */
bool validCoordsetGroup( const sidre::Group* coordset );

/*!
 * \brief Returns the requested mesh topology group.
 *
 * \param [in] group the blueprint mesh root group.
 * \param [in] topo optional argument to specify the topology name.
 *
 * \note If a topology is not specified with the second argument, the
 *  code returns the first group under the topologies group.
 *
 * \return topo pointer to the requested topology group or AXOM_NULLPTR, if
 *  the requested topology group does not exists or is invalid.
 *
 * \pre group != AXOM_NULLPTR
 * \post topo != AXOM_NULLPTR
 * \post blueprint::validTopologyGrou( topo ) == true
 *
 */
const sidre::Group* getTopologyGroup( const sidre::Group* group,
                                      const std::string& topo="" );

/*!
 * \brief Returns the coordset group associated with the given topology group.
 *
 * \param [in] group the mesh blueprint root group.
 * \param [in] topology the topology group associated with the requested
 *  coordset.
 *
 * \return coordset pointer to the coordset associated with given topology,
 *  else false.
 *
 * \pre blueprint::validRootGroup( group )
 * \pre blueprint::validTopologyGroup( topology )
 * \post blueprint::validCoordsetGroup( coordset )
 */
const sidre::Group* getCoordsetGroup( const sidre::Group* group,
                                      const sidre::Group* topology );

/*!
 * \brief Returns the mesh type and dimension given a root group that conforms
 *  to the computational mesh blueprint.
 *
 * \param [out] mesh_type the mint mesh type
 * \param [out] dimension the mesh dimension
 * \param [in] group pointer to the root group
 * \param [in] topo optional argument corresponding to the mesh topology.
 *
 * \note If a topology is not specified, the code assumes that the first group
 *  under the 'topologies' group corresponds to the mesh topology.
 *
 * \pre group != AXOM_NULLPTR
 *
 * \see MeshTypes
 */
/// @{

void getMeshTypeAndDimension( int& mesh_type, int& dimension,
                              const sidre::Group* group,
                              const std::string& topology );

void getMeshTypeAndDimension( int& mesh_type, int& dimension,
                              const sidre::Group* group );

/// @}

} /* namespace blueprint */

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_MESH_BLUEPRINT_HPP_ */

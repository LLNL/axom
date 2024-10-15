// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_MESH_HELPERS__HPP_
#define AXOM_QUEST_MESH_HELPERS__HPP_

#include "axom/config.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"

#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#endif

/*!
 * \brief Some helper functions to create and process mesh data
 */

namespace axom
{
namespace quest
{
namespace util
{
#ifdef AXOM_USE_MFEM
/*!
 * /brief Creates a 2D Cartesian mfem mesh lying within a given bounding box
 *
 * \param bbox The bounding box for the mesh
 * \param res The resolution of the mesh
 * \param polynomial_order The polynomial order for the mesh curvature
 * \param reorder_space_filling Reorder the mesh elements according to a space filling curve
 *
 * \return A pointer to the generated mesh
 * \note The user will be responsible for the mesh's memory
 */
mfem::Mesh* make_cartesian_mfem_mesh_2D(const primal::BoundingBox<double, 2>& bbox,
                                        const primal::NumericArray<int, 2>& res,
                                        int polynomial_order,
                                        bool reorder_space_filling = true);

/*!
 * /brief Creates a 3D Cartesian mfem mesh lying within a given bounding box
 *
 * \param bbox The bounding box for the mesh
 * \param res The resolution of the mesh
 * \param polynomial_order The polynomial order for the mesh curvature
 * \param reorder_space_filling Reorder the mesh elements according to a space filling curve
 *
 * \return A pointer to the generated mesh
 * \note The user will be responsible for the mesh's memory
 */
mfem::Mesh* make_cartesian_mfem_mesh_3D(const primal::BoundingBox<double, 3>& bbox,
                                        const primal::NumericArray<int, 3>& res,
                                        int polynomial_order,
                                        bool reorder_space_filling = true);
#endif

#if defined(AXOM_USE_SIDRE)
/*!
  /brief Creates a 3D Cartesian mfem mesh lying within a given bounding box

  \param meshGrp Put the mesh in this Group
  \param bbox The bounding box for the mesh
  \param res The resolution of the mesh
  \param topologyName Name of the blueprint topoloyy to use.
  \param coordsetName Name of the blueprint coordset to use.

  \return The meshGrp pointer
*/
axom::sidre::Group* make_structured_blueprint_box_mesh(
  axom::sidre::Group* meshGrp,
  const primal::BoundingBox<double, 3>& bbox,
  const primal::NumericArray<int, 3>& res,
  const std::string& topologyName = "mesh",
  const std::string& coordsetName = "coords");

axom::sidre::Group* make_unstructured_blueprint_box_mesh(
  axom::sidre::Group* meshGrp,
  const primal::BoundingBox<double, 3>& bbox,
  const primal::NumericArray<int, 3>& res,
  const std::string& topologyName = "mesh",
  const std::string& coordsetName = "coords");

void convert_blueprint_structured_explicit_to_unstructured(
  axom::sidre::Group* meshGrp,
  const std::string& topoName);

  #if defined(AXOM_USE_CONDUIT)
/*!
  \brief Check if blueprint mesh is valid.
*/
bool verifyBlueprintMesh(const axom::sidre::Group* meshGrp, conduit::Node info);
  #endif

#endif

}  // namespace util
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_MESH_HELPERS__HPP_

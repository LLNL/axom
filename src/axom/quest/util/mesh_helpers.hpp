// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_QUEST_MESH_HELPERS__HPP_
#define AXOM_QUEST_MESH_HELPERS__HPP_

#include "axom/config.hpp"
#include "axom/primal.hpp"

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

}  // namespace util
}  // namespace quest
}  // namespace axom

#endif  // AXOM_QUEST_MESH_HELPERS__HPP_

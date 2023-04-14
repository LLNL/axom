// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file MarchingCubes.hpp
 *
 * \brief Consists of class implementing marching cubes algorithm to
 * compute iso-surface from a scalar field in a blueprint mesh.
 */

#ifndef AXOM_PRIMAL_MARCHINGCUBES_H_
#define AXOM_PRIMAL_MARCHINGCUBES_H_

// Axom includes
#include "axom/config.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"

// Conduit includes
#include "conduit_node.hpp"

// C++ includes
#include <string>

namespace axom
{
namespace quest
{
/*!
  @brief Enum for runtime execution policy

  The policy implicitly selects the execution space and allocator id.
*/
enum class MarchingCubesRuntimePolicy
{
  seq = 0,
  omp = 1,
  cuda = 2
};

template <int DIM, typename ExecSpace = axom::SEQ_EXEC>
struct MarchingCubesImpl;

/*!
 * \@brief Class implementing marching cubes algorithm on a single
 * structured mesh domain within a multi-domain mesh.
 *
 * Implementation is for 2D (marching squares) and 3D (marching
 * cubes).
 *
 * The input mesh is a Conduit::Node following the Mesh Blueprint
 * convention.
 *
 * Usage example:
 * @beginverbatim
 *   MarchingCubes computationalMesh(meshNode, "coords");
 *   axom::mint::UnstructuredMesh<axom::mint::_SINGLe_SHAPE>
 *     surfaceMesh(3, min::CellType::Triangle);
 *   double contourValue = 0.0;
 *   set_function_field("my_function");
 *   mc.compute_iso_surface(contourValue);
 *   mc.populate_surface_mesh( surfaceMesh, "cellIdField");
 * @endverbatim
 *
 * To avoid confusion between the two meshes, we refer to the mesh
 * with the computational data as "parent" and the generated mesh the
 * "surface".
 *
 * \sa MarchingCubes
 */
class MarchingCubesSingleDomain
{
public:
  /*!
   * \brief Constructor for applying algorithm in multi-level mesh context.
   *
   * \param [in] dom Blueprint single-domain mesh containing scalar field.
   * \param [in] coordsetName Name of blueprint point coordinate set.
   * \param [in] maskField Cell-based std::int32_t mask field.  If provided,
   *             cells where this field evaluates to false are skipped.
   *
   * Some data from \a dom may be cached by the constructor.
   * Any change to it after the constructor leads to undefined behavior.
   * See set_domain(const conduit::Node &) for requirements on \a dom.
   *
   * The mesh coordinates should be contiguous.  See
   * conduit::blueprint::is_contiguous().  In the future, this
   * requirement may be relaxed, possibly at the cost of a
   * transformation and storage of the temporary contiguous layout.
   */
  MarchingCubesSingleDomain(const conduit::Node &dom,
                            const std::string &coordsetName,
                            const std::string &maskfield);

  //! @brief Spatial dimension of domain.
  int dimension() const { return m_ndim; }

  /*!
    @brief Set the field containing the nodal function.
    \param [in] fcnField Name of node-based scalar function values.
  */
  void set_function_field(const std::string &fcnField);

  /*!
   * \brief Computes the iso-surface.
   *
   * \param [in] contourVal iso-contour value
   *
   * Compute iso-surface using the marching cubes algorithm.
   */
  void compute_iso_surface(double contourVal = 0.0);

  //!@brief Get number of cells in the generated contour mesh.
  axom::IndexType get_surface_cell_count() const;
  //!@brief Get number of nodes in the generated contour mesh.
  axom::IndexType get_surface_node_count() const
  {
    return m_ndim * get_surface_cell_count();
  }

  /*!
    @brief Put generated surface in a mint::UnstructuredMesh.
    @param mesh Output mesh
    @param cellIdField Name of field to store the prent cell ids.
      If omitted, the data is not copied.
  */
  void populate_surface_mesh(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> &mesh,
    const std::string &cellIdField = {});

private:
  // MarchingCubesRuntimePolicy m_runtimePolicy;
  /*!
    \brief Computational mesh as a conduit::Node.

    Plan for supporting blueprint mesh in sidre is to shallow-copy the
    mesh to a conduit::Node on the heap, and point m_dom to it.
  */
  const conduit::Node *m_dom;
  int m_ndim;

  const std::string m_coordsetPath;
  std::string m_fcnPath;
  std::string m_maskPath;

  std::shared_ptr<MarchingCubesImpl<2, axom::execution_space<axom::SEQ_EXEC>>> m_impl2d;
  std::shared_ptr<MarchingCubesImpl<3, axom::execution_space<axom::SEQ_EXEC>>> m_impl3d;

  /*!
   * \brief Set the blueprint single-domain mesh.
   *
   * Some data from \a dom may be cached by the constructor.
   */
  void set_domain(const conduit::Node &dom);

};  // class MarchingCubesSingleDomain

/*!
 * \@brief Class implementing marching cubes algorithm.
 *
 * \sa MarchingCubesSingleDomain
 */
class MarchingCubes
{
public:
  /*!
   * \brief Constructor sets up computational mesh and data for running the
   * marching cubes algorithm.
   *
   * \param [in] bpMesh Blueprint multi-domain mesh containing scalar field.
   * \param [in] coordsetName Name of blueprint point coordinate set.
   * \param [in] maskField Cell-based std::int32_t mask field.  If provided,
   *             cells where this field evaluates to false are skipped.
   *
   * Some data from \a bpMesh may be cached by the constructor.
   * Any change to it after the constructor leads to undefined behavior.
   *
   * The mesh coordinates should be contiguous.  See
   * conduit::blueprint::is_contiguous().  In the future, this
   * requirement may be relaxed, possibly at the cost of a
   * transformation and storage of the temporary contiguous layout.
   */
  MarchingCubes(const conduit::Node &bpMesh,
                const std::string &coordsetName,
                const std::string &maskField = {});

  /*!
    @brief Set the field containing the nodal function.
    \param [in] fcnField Name of node-based scalar function values.
  */
  void set_function_field(const std::string &fcnField);

  /*!
   \brief Computes the iso-surface.
   \param [in] contourVal iso-contour value
   */
  void compute_iso_surface(double contourVal = 0.0);

  /*!
    @brief Put generated surface in a mint::UnstructuredMesh.
    @param mesh Output mesh
    @param cellIdField Name of field to store the (axom::IndexType)
      parent cell ids. If omitted, the data is not provided.
    @param domainIdField Name of field to store the (axom::IndexType)
      parent domain ids. If omitted, the data is not provided.

    If the fields aren't in the mesh, they will be created.
  */
  void populate_surface_mesh(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> &mesh,
    const std::string &cellIdField = {},
    const std::string &domainIdField = {});

private:
  MarchingCubesRuntimePolicy m_runtimePolicy;

  //! @brief Single-domain implementations.
  axom::Array<std::shared_ptr<MarchingCubesSingleDomain>> m_singles;
  int m_ndim;
  const std::string m_coordsetPath;
  std::string m_fcnPath;
  std::string m_maskPath;

  /*
    Non-state variables we don't want to have to propagate down the
    stack repeatedly during computation.
  */
  mutable axom::mint::Mesh *m_surfaceMesh;

  // Use simple pointers for now.  Later, maybe sidre.

  void set_mesh(const conduit::Node &bpMesh);
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_PRIMAL_ISOSURFACE_H_

// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file MarchingCubesAlgo.hpp
 *
 * \brief Consists of class implementing marching cubes algorithm to
 * compute iso-surface from a scalar field in a blueprint mesh.
 */

#ifndef AXOM_PRIMAL_MARCHINGCUBES_H_
#define AXOM_PRIMAL_MARCHINGCUBES_H_

// Axom includes
#include "axom/config.hpp"
#include "axom/mint/mesh/Mesh.hpp"

// Conduit includes
#include "conduit_node.hpp"

// C++ includes
#include <string>

namespace axom
{
namespace quest
{
/*!
 * \@brief Class implementing marching cubes algorithm on a single
 * structured mesh domain.
 *
 * Implementation is for 2D (marching squares) and 3D.
 *
 * \sa MarchingCubesAlgo
 */
class MarchingCubesAlgo1
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
  MarchingCubesAlgo1(const conduit::Node &dom,
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
    @brief Set the output surface mesh object.
  */
  void set_output_mesh(axom::mint::Mesh *surfaceMesh);

  /*!
    @brief Save the originating cell index in output mesh.

    For each contour mesh cell generated, the 1D index of the
    generating computational mesh cell is set in the \a cellIdField of
    the contour mesh.  If the field doesn't exist, it will be created.
    To disable, set \a outputCellIds to empty.
  */
  void set_cell_id_field(const std::string &cellIdField)
  {
    m_cellIdField = cellIdField;
  }

  /*!
   * \brief Computes the iso-surface.
   *
   * \param [in] contourVal iso-contour value
   *
   * Compute iso-surface using the marching cubes algorithm.
   */
  void compute_iso_surface(double contourVal = 0.0);

private:
  /*!
    \brief Computational mesh as a conduit::Node.

    Plan for supporting blueprint mesh in sidre is to shallow-copy the
    mesh to a conduit::Node on the heap, and point m_dom to it.
  */
  const conduit::Node *m_dom;
  int m_ndim;
  axom::Array<axom::IndexType> m_cShape;  //!< @brief Shape of cell-centered data arrays in m_dom.
  axom::Array<axom::IndexType> m_logicalOrigin;  //!< @brief First domain cell in each direction.

  const std::string m_coordsetPath;
  std::string m_fcnPath;
  std::string m_maskPath;

  /*
    Non-state variables we don't want to have to propagate down the
    stack repeatedly during computation.
  */
  mutable axom::mint::Mesh *m_surfaceMesh;

  /* @brief Field name for recording computational-mesh cell id.

     This is for looking up the computational-mesh cell id that
     contributed a surface-mesh cell.
  */
  std::string m_cellIdField;

  double _contourVal;

  /*!
   * \brief Set the blueprint single-domain mesh.
   *
   * Some data from \a dom may be cached by the constructor.
   */
  void set_domain(const conduit::Node &dom);

  void contourCell2D(double xx[4], double yy[4], double cellValues[4]);

  void contourCell3D(double xx[8], double yy[8], double zz[8], double f[8]);

  void linear_interp(int edgeIdx,
                     const double *xx,
                     const double *yy,
                     const double *zz,
                     const double *nodeValues,
                     double *xyz);

  int computeIndex(const double *f);

};  // class MarchingCubesAlgo1

/*!
 * \@brief Class implementing marching cubes algorithm.
 *
 * \sa MarchingCubesAlgo1
 */
class MarchingCubesAlgo
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
  MarchingCubesAlgo(const conduit::Node &bpMesh,
                    const std::string &coordsetName,
                    const std::string &maskField = {});

  /*!
    @brief Set the field containing the nodal function.
    \param [in] fcnField Name of node-based scalar function values.
  */
  void set_function_field(const std::string &fcnField);

  /*!
    @brief Set the output surface mesh object.
  */
  void set_output_mesh(axom::mint::Mesh *surfaceMesh);

  /*!
    @brief Save the originating 1D cell index as an int in the specified
    cell-centered field of the mesh.

    If the field doesn't exist, it will be created.  To disable, set
    \a domainIdField to empty.
  */
  void set_cell_id_field(const std::string &cellIdField)
  {
    m_cellIdField = cellIdField;
    for(auto &s : m_sd)
    {
      s->set_cell_id_field(cellIdField);
    }
  }

  /*!
    @brief Save the originating domain index as an int in the specified
    cell-centered field of the mesh.

    If the field doesn't exist, it will be created.  To disable, set
    \a domainIdField to empty.
  */
  void set_domain_id_field(const std::string &domainIdField)
  {
    m_domainIdField = domainIdField;
  }

  /*!
 * \brief Computes the iso-surface.
 *
 * \param [in] contourVal iso-contour value
 */
  void compute_iso_surface(double contourVal = 0.0);

private:
  //! @brief Single-domain implementations.
  axom::Array<std::shared_ptr<MarchingCubesAlgo1>> m_sd;
  int m_ndim;
  const std::string m_coordsetPath;
  std::string m_fcnPath;
  std::string m_maskPath;

  /*
    Non-state variables we don't want to have to propagate down the
    stack repeatedly during computation.
  */
  mutable axom::mint::Mesh *m_surfaceMesh;

  //! @brief See MarchingCubesAlgo1::m_cellIdField.
  std::string m_cellIdField;
  /* @brief Field name for recording computational-mesh domain id.

     This is for looking up the computational-mesh domain id that
     contributed a surface-mesh cell.
  */
  std::string m_domainIdField;

  // Use simple pointers for now.  Later, maybe sidre.

  void set_mesh(const conduit::Node &bpMesh);
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_PRIMAL_ISOSURFACE_H_

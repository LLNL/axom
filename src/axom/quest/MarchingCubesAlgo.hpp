// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file MarchingCubesAlgo.hpp
 *
 * \brief Consists of class implementing marching cubes algorithm to
 * compute iso-surface from a scalar field in a blueprint mesh.
 */

#ifndef AXOM_PRIMAL_ISOSURFACE_H_
#define AXOM_PRIMAL_ISOSURFACE_H_

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
   * \param [in] fieldName Name of scalar field data.
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
                     const std::string &fcnField,
                     const std::string &maskfield);

  //! @brief Spatial dimension of domain.
  int dimension() const
  {
    return _ndim;
  }

  /*!
    @brief Set the output sufrace mesh object.
  */
  void set_output_mesh(axom::mint::Mesh *surfaceMesh);

  /*!
    @brief Save the originating cell index in output mesh.

    For each contour mesh cell generated, the 1D index of the
    generating computational mesh cell is set in the \a cellIdField of
    the contour mesh.  If the field doesn't exist, it will be created.
    To disable, set \a outputCellIds to empty.
  */
  void set_cell_id_field(const std::string& cellIdField)
  {
    _cellIdField = cellIdField;
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
    mesh to a conduit::Node on the heap, and point _dom to it.
  */
  const conduit::Node *_dom;
  int _ndim;
  axom::Array<axom::IndexType> _logicalSize; //! @brief Number of cells in each direction
  axom::Array<axom::IndexType> _logicalOrigin; //! @brief First domain cell in each direction.

  const std::string _coordsetPath;
  const std::string _fcnPath;
  const std::string _maskPath;

  /*
    Non-state variables we don't want to have to propagate down the
    stack repeatedly during computation.
  */
  mutable axom::mint::Mesh* _surfaceMesh;

  std::string _cellIdField;

  double _contourVal;

  /*!
   * \brief Set the blueprint single-domain mesh.
   *
   * Some data from \a dom may be cached by the constructor.
   */
  void set_domain(const conduit::Node &dom);

  void contourCell2D( double xx[4], double yy[4], double cellValues[4]);

  void contourCell3D( double xx[8], double yy[8], double zz[8], double f[8]);

  void linear_interp( int edgeIdx,
                      const double* xx, const double* yy, const double* zz,
                      const double* f, double* xyz );

  int computeIndex( const double* f );

}; // class MarchingCubesAlgo1

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
   * \param [in] fieldName Name of node-based scalar function values.
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
                    const std::string &fcnField,
                    const std::string &maskField={});

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
  void set_cell_id_field(const std::string& cellIdField)
  {
    _cellIdField = cellIdField;
    for(auto &s : _sd)
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
  void set_domain_id_field(const std::string& domainIdField)
  {
    _domainIdField = domainIdField;
  }


/*!
 * \brief Computes the iso-surface.
 *
 * \param [in] contourVal iso-contour value
 */
  void compute_iso_surface(double contourVal = 0.0);

private:
  //! @brief Single-domain implementations.
  axom::Array<std::shared_ptr<MarchingCubesAlgo1>> _sd;
  int _ndim;
  const std::string _coordsetPath;
  const std::string _fcnPath;
  const std::string _maskPath;

  /*
    Non-state variables we don't want to have to propagate down the
    stack repeatedly during computation.
  */
  mutable axom::mint::Mesh *_surfaceMesh;

  std::string _cellIdField;
  std::string _domainIdField;

  // Use simple pointers for now.  Later, maybe sidre.

  void set_mesh(const conduit::Node &bpMesh);
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_PRIMAL_ISOSURFACE_H_

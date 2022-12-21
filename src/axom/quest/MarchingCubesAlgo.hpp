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
   * \param [in] fieldName Name of scalar field data.
   * \param [in] maskField Cell-based mask field.  If provided, only cells
   *             where this field evaluates to true are evaluated.
   *
   * See set_domain(const conduit::Node &) for requirements on \a dom.
   */
  MarchingCubesAlgo1(const conduit::Node &dom,
                     const std::string &valueField,
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
    @brief Set a container for saving the output cell ids.

    For each surface mesh cell generated, the generating computational mesh
    cell id is appended to \a outputCellIds.  To disable, set \a outputCellIds to nullptr.
  */
  void set_output_cell_ids(axom::Array<int>* outputCellIds)
  {
    _outputCellIds = outputCellIds;
  }

  /*!
   * \brief Computes the iso-surface.
   *
   * \param [in] isoValue iso-contour value
   *
   * Compute iso-surface using the marching cubes algorithm.
   */
  void compute_iso_surface(double isoValue = 0.0);

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

  const std::string _valueField;
  const std::string _maskField;

  /*
    Non-state variables we don't want to have to propagate down the
    stack repeatedly during computation.
  */
  mutable axom::mint::Mesh* _surfaceMesh;
  mutable axom::Array<int>* _outputCellIds;
  double _isoValue{0.0};

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
 * \param [in] fieldName Name of node-based scalar function values.
 * \param [in] maskField Cell-based mask field.  If provided, cells
 *             where this field evaluates to false are skipped.
 *
 * Some data from \a bpMesh may be cached by the constructor.
 */
  MarchingCubesAlgo(const conduit::Node &bpMesh,
                    const std::string &valueField,
                    const std::string &maskField={});

  /*!
    @brief Set the output surface mesh object.
  */
  void set_output_mesh(axom::mint::Mesh *surfaceMesh);


  /*!
    @brief Set a container for saving the output domain ids.

    For each surface mesh cell generated, the generating computational mesh
    domain id is appended to \a outputDomainIds.
    To disable, set \a outputDomainIds to nullptr.
  */
  void set_output_domain_ids(axom::Array<int>* outputDomainIds)
  {
    _outputDomainIds = outputDomainIds;
  }


  /*!
    @brief Set the container for saving the output cell ids.

    For each surface mesh cell generated, the generating computational mesh
    cell id is appended to \a outputCellIds.  To disable, set \a outputCellIds to nullptr.
  */
  void set_output_cell_ids(axom::Array<int>* outputCellIds)
  {
    _outputCellIds = outputCellIds;
    for(auto &s : _sd)
    {
      s->set_output_cell_ids(outputCellIds);
    }
  }


/*!
 * \brief Computes the iso-surface.
 *
 * \param [in] isoValue iso-contour value
 *
 * Compute iso-surface using the marching cubes algorithm.
 */
  void compute_iso_surface(double isoValue = 0.0);

private:
  //! @brief Single-domain implementations.
  axom::Array<std::shared_ptr<MarchingCubesAlgo1>> _sd;
  int _ndim;
  const std::string _valueField;
  const std::string _maskField;

  /*
    Non-state variables we don't want to have to propagate down the
    stack repeatedly during computation.
  */
  mutable axom::mint::Mesh *_surfaceMesh;
  mutable axom::Array<int>* _outputCellIds;
  mutable axom::Array<int>* _outputDomainIds;

  // Use simple pointers for now.  Later, maybe sidre.

  void set_mesh(const conduit::Node &bpMesh);
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_PRIMAL_ISOSURFACE_H_

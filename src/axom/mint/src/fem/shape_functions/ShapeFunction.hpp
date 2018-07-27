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

#ifndef MINT_SHAPEFUNCTION_HPP_
#define MINT_SHAPEFUNCTION_HPP_

namespace axom
{
namespace mint
{

/*!
 * \brief The ShapeFunction class defines the shape functions, \f$ N(\xi)^e \f$,
 *  for a given reference element \f$ \overline{\Omega^e} \f$
 *
 *  ShapeFunction is a templated class that employs the Curiously Recurring
 *  Template Pattern (CRTP) to provide a unified API for all underlying
 *  ShapeFunction implementations.
 *
 * \see FiniteElement
 * \see FEBasisTypes.hpp
 */
template < typename ShapeType >
class ShapeFunction
{
public:

/// \name Generic ShapeFunction API
///@{

  /*!
   * \brief Returns the underlying cell type, e.g., MINT_QUAD, etc.
   * \return cellType the cell type
   * \see CellType.hpp
   */
  static CellType cellType( ) { return ShapeType::getCellType(); };

  /*!
   * \brief Returns the Finite Element basis family type
   * \return type the basis function type, e.g., MINT_LAGRANGE_BASIS
   * \see FEBasisTypes.hpp
   */
  static int type( ) { return ShapeType::getType(); };

  /*!
   * \brief Returns the number of degrees of freedom
   * \return ndofs the number of degrees of freedom
   * \post ndofs >= 1
   */
  static int numDofs( ) { return ShapeType::getNumDofs(); };

  /*!
   * \brief Returns the maximum number of iterations for the Newton-Raphson
   * \return N the maximum number of Newton-Raphson iterations
   */
  static int maxNewtonIters( ) { return ShapeType::getMaxNewtonIters( ); };

  /*!
   * \brief Returns the dimension of the reference element
   * \return ndims the dimension of the reference element
   * \post ndims >= 1
   */
  static int dimension( ) { return ShapeType::getDimension( ); };

  /*!
   * \brief Returns the min coordinate of the reference element
   * \return min the min coordinate of the reference element
   */
  static double min( ) { return ShapeType::getMin( ); };

  /*!
   * \brief Returns the max coordinate of the reference element
   * \return max the max coordinate of the reference element
   */
  static double max( ) { return ShapeType::getMax( ); };

  /*!
   * \brief Returns the center of the reference element.
   *
   * \param [out] center buffer (ndims long) to store the centroid
   * \pre center != AXOM_NULLPTR
   */
  static void center( double* center ) { ShapeType::getCenter( center ); };

  /*!
   * \brief Returns the coordinates of the reference element.
   *
   * \param [out] coords buffer (ndims*ndofs long) to store the coordinates
   * \pre coords != AXOM_NULLPTR
   *
   * \note THe coordinates are arranged in column-major flat array layout
   */
  static void coords( double* coords ) { ShapeType::getCoords( coords ); };

  /*!
   * \brief Evaluates the ShapeFunction at the given natural coordinates.
   *
   * \param [in] nc natural coordinates at which to compute the shape functions
   * \param [out] phi buffer (ndofs long) to store the shape function.
   *
   * \pre nc != AXOM_NULLPTR
   * \pre phi != AXOM_NULLPTR
   */
  static void evaluate( const double* nc, double* phi )
  {
    ShapeType::computeShape( nc, phi );
  }

  /*!
   * \brief Evaluates the first derivatives of the shape function at the given
   *  natural coordinates.
   *
   * \param [in] nc natural coordinates at which to compute the derivatives
   * \param [out] phidot buffer (ndofs*ndims long) to store the derivatives
   *
   * \pre nc != AXOM_NULLPTR
   * \pre phidot != AXOM_NULLPTR
   */
  static void derivatives( const double* nc, double* phidot )
  {
    ShapeType::computeDerivatives( nc, phidot );
  }

///@}
};

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_SHAPEFUNCTION_HPP_ */

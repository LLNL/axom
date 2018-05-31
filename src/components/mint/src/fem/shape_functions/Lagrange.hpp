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

#ifndef MINT_LAGRANGE_SHAPEFUNCTION_HPP_
#define MINT_LAGRANGE_SHAPEFUNCTION_HPP_

// Axom includes
#include "axom/Macros.hpp" // For AXOM_STATIC_ASSERT(), AXOM_NOT_USED()

// Mint includes
#include "mint/CellTypes.hpp"
#include "mint/FEBasisTypes.hpp"
#include "mint/ShapeFunction.hpp"

// Slic includes
#include "slic/slic.hpp"

namespace axom
{
namespace mint
{

/*!
 * \brief Defines the Lagrange family of Finite Elements
 *
 * \tparam CellType the cell type of the element, e.g., mint::QUAD, etc.
 *
 * \note This is the default implementation. Only stubs are defined at this
 *  level.This class is specialized according to cell type.
 */
template < int CellType >
class Lagrange : public ShapeFunction< Lagrange< CellType > >
{

public:

  /*!
   * \brief Returns the cell type of this instance.
   * \return c the cell type, e.g., MINT_QUAD etc.
   * \see CellTypes.hpp
   *
   * \note This method is implemented in specialized instances.
   */
  static int getCellType()
  {
    AXOM_STATIC_ASSERT( (CellType >= 0) && (CellType < mint::NUM_CELL_TYPES) );
    Lagrange< CellType >::checkCellType();
    return CellType;
  }

  /*!
   * \brief Returns the Basis function type.
   * \return b the basis thpe, e.g., MINT_LAGRANGE_BASIS
   * \see FEBasisTypes.hpp
   *
   * \note This method is implemented in specialized instances.
   */
  static int getType()
  {
    AXOM_STATIC_ASSERT( (CellType >= 0) && (CellType < mint::NUM_CELL_TYPES) );
    Lagrange< CellType >::checkCellType();
    return MINT_UNDEFINED_BASIS;
  }

  /*!
   * \brief Returns the number of degrees of freedom on this Finite Element.
   * \return nodfs the number of degrees of freedom.
   *
   * \note This method is implemented in specialized instances.
   */
  static int getNumDofs()
  {
    AXOM_STATIC_ASSERT( (CellType >= 0) && (CellType < mint::NUM_CELL_TYPES) );
    Lagrange< CellType >::checkCellType();
    return 0;
  }

  /*!
   * \brief Returns the max number of newton iterations for this Finite Element
   * \return N the max number of newton iterations
   * \post N >= 1
   *
   * \note This method is implemented in specialized instances.
   */
  static int getMaxNewtonIters()
  {
    AXOM_STATIC_ASSERT( (CellType >= 0) && (CellType < mint::NUM_CELL_TYPES) );
    Lagrange< CellType >::checkCellType();
    return 0;
  }

  /*!
   * \brief Returns the dimension of the reference element.
   * \return dim the dimension of the reference element
   * \post dim >= 1 && dim <= 3
   *
   * \note This method is implemented in specialized instances.
   */
  static int getDimension()
  {
    AXOM_STATIC_ASSERT( (CellType >= 0) && (CellType < mint::NUM_CELL_TYPES) );
    Lagrange< CellType >::checkCellType();
    return 0;
  }

  /*!
   * \brief Returns the min coordinate of the element's reference space
   * \return min the min coordinate of the element's reference space, e.g., 0
   *
   * \note This method is implemented in specialized instances.
   */
  static int getMin()
  {
    AXOM_STATIC_ASSERT( (CellType >= 0) && (CellType < mint::NUM_CELL_TYPES) );
    Lagrange< CellType >::checkCellType();
    return 0;
  }

  /*!
   * \brief Returns the max coordinate of the element's reference space
   * \return max the max coordinate of the element's reference space, e.g., 1
   *
   * \note This method is implemented in specialized instances.
   */
  static int getMax()
  {
    AXOM_STATIC_ASSERT( (CellType >= 0) && (CellType < mint::NUM_CELL_TYPES) );
    Lagrange< CellType >::checkCellType();
    return 0;
  }

  /*!
   * \brief Returns the centroid of the reference element.
   *
   * \param [out] center ndims long buffer to store the centroid
   * \pre center != AXOM_NULLPTR
   *
   * \note This method is implemented in specialized instances.
   */
  static void getCenter( double* AXOM_NOT_USED(center) )
  {
    AXOM_STATIC_ASSERT( (CellType >= 0) && (CellType < mint::NUM_CELL_TYPES) );
    Lagrange< CellType >::checkCellType();
  }

  /*!
   * \brief Returns the coordinates of the reference element.
   *
   * \param [out] coords ndims*ndofs long buffer to store the coordinates
   * \pre coords != AXOM_NULLPTR
   *
   * \note The coordinates are arranged in column-major flat array layout.
   * \note This method is implemented in specialized instances.
   */
  static void getCoords( double* AXOM_NOT_USED(coords) )
  {
    AXOM_STATIC_ASSERT( (CellType >= 0) && (CellType < mint::NUM_CELL_TYPES) );
    Lagrange< CellType >::checkCellType();
  }

  /*!
   * \brief Computes the shape functions of the Finite Element at the given
   *  natural coordinates, \f$ \xi \f$
   *
   * \param [in] nc natural coordinates at which to compute the shape functions
   * \param [out] phi buffer (ndofs long) to store the shape functions
   *
   * \pre nc != AXOM_NULLPTR
   * \pre phi != AXOM_NULLPTR
   *
   * \note This method is implemented in specialized instances.
   */
  static void computeShape( const double* AXOM_NOT_USED(nc),
                            double* AXOM_NOT_USED(phi) )
  {
    AXOM_STATIC_ASSERT( (CellType >= 0) && (CellType < mint::NUM_CELL_TYPES) );
    Lagrange< CellType >::checkCellType();
  }

  /*!
   * \brief Computes the shape function first derivatives for the Finite
   *  Element at the given natural coordinates, \f$ \xi \f$
   *
   * \param [in] nc natural coordinates at which to compute the derivatives
   * \param [out] phidot buffer (ndofs*ndims long) for the derivatives
   *
   * \pre nc != AXOM_NULLPTR
   * \pre phidot != AXOM_NULLPTR
   *
   * \note This method is implemented in specialized instances.
   */
  static void computeDerivatives( const double* AXOM_NOT_USED(nc),
                                  double* AXOM_NOT_USED(phidot) )
  {
    AXOM_STATIC_ASSERT( (CellType >= 0) && (CellType < mint::NUM_CELL_TYPES) );
    Lagrange< CellType >::checkCellType();
  }

private:

  /*!
   * \brief Checks if the CellType is valid and supported in the Lagrange basis.
   */
  static void checkCellType( )
  {
    if ( (CellType >= 0) && ( CellType < mint::NUM_CELL_TYPES ) )
    {
      SLIC_ERROR( "Lagrange shape functions for [" <<
                  cell_info[ CellType ].name  << "] are not defined!" );
    }
    else
    {
      SLIC_ERROR( "Invalid CellType: " << CellType );
    }

  }

};

} /* namespace mint */
} /* namespace axom */

// Template Specializations for Lagrange
#include "mint/lagrange_hexa_27.hpp"
#include "mint/lagrange_hexa_8.hpp"
#include "mint/lagrange_prism_6.hpp"
#include "mint/lagrange_pyra_5.hpp"
#include "mint/lagrange_quad_4.hpp"
#include "mint/lagrange_quad_9.hpp"
#include "mint/lagrange_tetra_4.hpp"
#include "mint/lagrange_tri_3.hpp"

#endif /* MINT_LAGRANGE_SHAPEFUNCTION_HPP_ */

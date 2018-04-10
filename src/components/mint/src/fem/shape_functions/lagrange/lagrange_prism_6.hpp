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

#ifndef MINT_PRISM_6_HPP_
#define MINT_PRISM_6_HPP_

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

static const double PRISM_ONE_THIRD =  1.0/3.0;

/*!
 * \brief Lagrange Finite Element definition for the Linear Prism
 *
 * \verbatim
 *
 * prism_6:
 *
 *       1_ _ _ _ _ _ _4
 *      / \           / \
 *     /   \         /   \
 *    /_ _ _\_ _ _ _/_ _ _\
 *   0       2      3      5
 *
 * \endverbatim
 *
 * \see ShapeFunction
 */
template < >
class Lagrange< mint::PRISM > :
  public ShapeFunction< Lagrange< mint::PRISM > >
{
public:

  static int getCellType() { return mint::PRISM; }

  static int getType() { return MINT_LAGRANGE_BASIS; }

  static int getNumDofs() { return 6; }

  static int getMaxNewtonIters() { return 16; }

  static int getDimension() { return 3; }

  static double getMin() { return 0; }

  static double getMax() { return 1; }

  static void getCenter( double* center )
  {
    SLIC_ASSERT( center != AXOM_NULLPTR );
    center[ 0 ] = center[ 1 ] = PRISM_ONE_THIRD;
    center[ 2 ] = 0.5;
  }

  static void getCoords( double* coords )
  {
    SLIC_ASSERT( coords != AXOM_NULLPTR );

    // node 0
    coords[ 0 ] = 0.0;
    coords[ 1 ] = 0.0;
    coords[ 2 ] = 0.0;

    // node 1
    coords[ 3 ] = 1.0;
    coords[ 4 ] = 0.0;
    coords[ 5 ] = 0.0;

    // node 2
    coords[ 6 ] = 0.0;
    coords[ 7 ] = 1.0;
    coords[ 8 ] = 0.0;

    // node 3
    coords[ 9 ]  = 0.0;
    coords[ 10 ] = 0.0;
    coords[ 11 ] = 1.0;

    // node 4
    coords[ 12 ] = 1.0;
    coords[ 13 ] = 0.0;
    coords[ 14 ] = 1.0;

    // node 5
    coords[ 15 ] = 0.0;
    coords[ 16 ] = 1.0;
    coords[ 17 ] = 1.0;
  }

  static void computeShape( const double* xr, double* phi )
  {
    SLIC_ASSERT(  xr != AXOM_NULLPTR );
    SLIC_ASSERT(  phi != AXOM_NULLPTR );

    const double r  = xr[0];
    const double s  = xr[1];
    const double t  = xr[2];
    const double rm = 1.-r;
    const double tm = 1.-t;
    const double rm_minus_s = rm-s;

    phi[ 0 ] = rm_minus_s * tm;
    phi[ 1 ] = r * tm;
    phi[ 2 ] = s * tm;

    phi[ 3 ] = rm_minus_s * t;
    phi[ 4 ] = r * t;
    phi[ 5 ] = s * t;
  }

  static void computeDerivatives( const double* xr, double* phidot )
  {
    SLIC_ASSERT(  xr != AXOM_NULLPTR );
    SLIC_ASSERT(  phidot != AXOM_NULLPTR );

    const double r  = xr[0];
    const double s  = xr[1];
    const double t  = xr[2];

    // r-derivatives
    phidot[ 0 ] = -1.0 + t;
    phidot[ 1 ] =  1.0 - t;
    phidot[ 2 ] =  0.0;
    phidot[ 3 ] = -t;
    phidot[ 4 ] =  t;
    phidot[ 5 ] =  0.0;

    // s-derivatives
    phidot[ 6  ] = -1.0 + t;
    phidot[ 7  ] =  0.0;
    phidot[ 8  ] =  1.0 - t;
    phidot[ 9  ] = -t;
    phidot[ 10 ] = 0.0;
    phidot[ 11 ] = t;

    // t-derivatives
    phidot[ 12 ] = -1.0 + r + s;
    phidot[ 13 ] = -r;
    phidot[ 14 ] = -s;
    phidot[ 15 ] =  1.0 - r - s;
    phidot[ 16 ] =  r;
    phidot[ 17 ] =  s;
  }
};

} /* namespace mint */
} /* namespace axom */
#endif /* MINT_PRISM_6_HPP_ */

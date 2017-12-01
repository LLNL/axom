/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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

#ifndef AXOM_NUMERICS_POLY_SOLVE_HPP_
#define AXOM_NUMERICS_POLY_SOLVE_HPP_

#include "axom/Types.hpp" // for AXOM_NULLPTR

// C/C++ includes
#include <cassert> // for assert()

namespace axom
{
namespace numerics
{

/*!
 * \brief Solves a linear equation of the form \f$ ax + b = 0 \f$.
 *
 * \param [in] coeff Equation coefficients: coeff[i] multiplies \f$ x^i \f$.
 * \param [out] roots The roots of the equation.
 * \param [out] numRoots The number of roots found.
 * \return 0 for success
 *
 * \pre coeff and roots point to arrays of at least 2.
 */
int solve_linear( const double* coeff, double* roots, int& numRoots );

/*!
 * \brief Solves a quadratic equation of the form \f$ ax^2 + bx + c = 0 \f$.
 *
 * \param [in] coeff Equation coefficients: coeff[i] multiplies \f$ x^i \f$.
 * \param [out] roots The roots of the equation.
 * \param [out] numRoots The number of roots found.
 * \return 0 for success
 *
 * \pre coeff and roots point to arrays of at least 3.
 */
int solve_quadratic( const double* coeff, double* roots, int& numRoots );

/*!
 * \brief Solves a cubic equation of the form \f$ ax^3 + bx^2 + cx + d = 0 \f$.
 *
 * \param [in] coeff Equation coefficients: coeff[i] multiplies \f$ x^i \f$.
 * \param [out] roots The roots of the equation.
 * \param [out] numRoots The number of roots found.
 * \return 0 for success
 *
 * \pre coeff and roots point to arrays of at least 4.
 */
int solve_cubic( const double* coeff,  double* roots, int &numRoots );

} /* end namespace numerics */
} /* end namespace axom */

//------------------------------------------------------------------------------
// Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace numerics
{

int solve_linear( const double* coeff, double* roots, int& numRoots )
{
  int status = 0;

  return status;
}

int solve_quadratic( const double* coeff, double* roots, int& numRoots )
{
  int status = 0;

  return status;
}

int solve_cubic( const double* coeff, double* roots, int& numRoots )
{
  int status = 0;

  return status;
}

} /* end namespace numerics */
} /* end namespace axom */

#endif // AXOM_NUMERICS_POLY_SOLVE_HPP_

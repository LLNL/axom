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
#include "axom_utils/Utilities.hpp" // for isNearlyEqual()

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
  int status = -1;

  if (utilities::isNearlyEqual(coeff[1], 0.)) {
    if (utilities::isNearlyEqual(coeff[0], 0.)) {
      // Infinite solutions: a horizontal line on the X-axis.
      status = 0;
      numRoots = -1;
    } else {
      // No solutions: a horizontal line not on the X-axis.
      numRoots = 0;
    }
  } else {
    // One solution, where the line crosses the X-axis.
    status = 0;
    numRoots = 1;
    roots[0] = -coeff[0] / coeff[1];
  }

  return status;
}

int solve_quadratic( const double* coeff, double* roots, int& numRoots )
{
  int status = -1;

  if (utilities::isNearlyEqual(coeff[2], 0.)) {
    // If this system is nearly linear, solve it as such.
    return solve_linear(coeff, roots, numRoots);
  }

  double discriminant = coeff[1]*coeff[1] - 4*coeff[2]*coeff[0];

  if (utilities::isNearlyEqual(discriminant, 0.)) {
    // One unique real root
    status = 0;
    numRoots = 1;
    roots[0] = roots[1] = -coeff[1] / (2*coeff[2]);
  } else if (discriminant < 0) {
    // No unique real roots
    numRoots = 0;
  } else {
    // Two real roots
    status = 0;
    numRoots = 2;
    roots[0] = (-coeff[1] + std::sqrt(discriminant)) / (2*coeff[2]);
    roots[1] = (-coeff[1] - std::sqrt(discriminant)) / (2*coeff[2]);
  }

  return status;
}

int solve_cubic( const double* coeff, double* roots, int& numRoots )
{
  int status = -1;

  return status;
}

} /* end namespace numerics */
} /* end namespace axom */

#endif // AXOM_NUMERICS_POLY_SOLVE_HPP_

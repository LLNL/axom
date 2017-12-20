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
 * \brief Solves a quadratic equation of the form \f$ ax^2 + bx + c = 0 \f$,
 * using the quadratic formula.
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
 * A closed-form solution for cubic equations was published in Cardano's
 * *Ars Magna* of 1545.  This can be summarized as follows:
 *
 * Start with the cubic equation
 * 
 *   \f[ x^3 + bx^2 + cx + d = 0. \f]
 *
 * Define the following:
 *
 *   \f[ p = b^2 - 3c \f]
 *   \f[ q = -\frac{27}{2}d - b^3 + \frac{9}{2}cb \f]
 *   \f[ t = 2p^{-3/2}q \f]
 *   \f[ y = \frac{3}{\sqrt{p}}(x + \frac{1}{3}b) \f]
 *
 * Then the original cubic equation can be written as \f$ y^3 - 3y = t \f$
 * with a solution \f$ y = \frac{1}{u} + u \f$, with
 *
 *   \f[ u = \sqrt[3]{\frac{t}{2} \pm \sqrt{\frac{t^2}{4} - 1}}. \f]
 *
 * Because the cubic is an odd function, there can be either one, two, or three
 * real roots.  In the case of one real root, the root can either be tripled or
 * single with two complex roots.  In the case of two real roots, one will be
 * doubled.  If the discriminant \f$ d = -27t^2 - 4(-3^3) \f$ is zero, the
 * equation has doubled or tripled roots.  If \f$ d > 0, \f$ there are three
 * distinct real roots, and if \f$ d > 0, \f$ there is one real root.
 *
 * See J. Kopp, Efficient numerical diagonalization of hermitian 3x3 matrices,
 * Int.J.Mod.Phys. C19:523-548, 2008 (https://arxiv.org/abs/physics/0610206)
 * and G. A. Korn and T. M. Korn, "Mathematical Handbook for Scientists and
 * Engineers," QA37 K84 1968 in the library; section 1.8-3 "Cubic Equations",
 * p. 23.
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
    // No real roots
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

inline double cuberoot(double x)
{
  // pow(x, y) returns NaN for negative finite x and noninteger y.
  if (x < 0) {
    return -pow(-x, 1. / 3.);
  } else {
    return  pow( x, 1. / 3.);
  }
}

int solve_cubic( const double* coeff, double* roots, int& numRoots )
{
  int status = -1;

  // Here I use variable names as presented in Korn & Korn:
  // x^3 + ax^2 + bx + c = 0
  double cubecoeff = coeff[3];
  double a = coeff[2];
  double b = coeff[1];
  double c = coeff[0];

  if (utilities::isNearlyEqual(cubecoeff, 0.)) {
    // If this system is nearly quadratic, solve it as such.
    return solve_quadratic(coeff, roots, numRoots);
  }

  a = a / cubecoeff;
  b = b / cubecoeff;
  c = c / cubecoeff;

  // Note that p and q differ by a multiplicative constant from Korn,
  // because they're always used with that multiplication.
  double p = (-a*a + 3*b) / 9;                 //  1/3 Korn's p
  double q = (a*(-2*a*a + 9*b) - 27*c) / 54;   // -1/2 Korn's q

  double Q = p*p*p + q*q;           // actual discriminant == -108Q
  // the term chVar occurs because we've changed variables
  // (x = y - a/3) and we need to change back to x.
  double chVar = -a / 3.;

  if (utilities::isNearlyEqual(Q, 0.)) {
    // We have three real roots, and at least two are equal.
    if (utilities::isNearlyEqual(q, 0.)) {
      numRoots = 1;
    } else {
      numRoots = 2;
    }
    status = 0;

    double cuberootq = cuberoot(q);
    roots[0] = chVar + 2*cuberootq;
    roots[1] = chVar -   cuberootq;
    roots[2] = roots[1];
  } else if (Q > 0) {
    // We have one real root, and two complex roots.
    // Right now we're calculating the real root only, but the complex roots
    // can be easily added by un-commenting the calculation.
    numRoots = 1;
    status = 0;

    double sqrtQ = sqrt(Q);
    double A = cuberoot(q + sqrtQ);
    double B = cuberoot(q - sqrtQ);

    roots[0] = chVar + A + B;
    roots[1] = 0;
    roots[2] = 0;
    // double imagpart = sqrt(3) * (A - B)/2;
    // double stavg = (A + B)/2;
    // Here note use of imaginary i:
    // roots[1] = chVar - stavg - i*imagpart;
    // roots[2] = chVar - stavg + i*imagpart;
  } else {
    // Q < 0, and we have three distinct real roots.
    // The case for Q == 0 is a special case of the case for Q > 0, and
    // both correspond to Cardano's method proper, reported in section
    // 1.8-3 of Korn, where the quadratic term ax^2 is eliminated by a
    // change of variable x = y - a/3.  For "irreducible" cubics (this
    // case, where the variable substitution didn't eliminate the quadratic),
    // we use the trigonometric solution to the cubic equation, reported
    // in section 1.8-4a of Korn.
    numRoots = 3;
    status = 0;

    double alpha = acos(q / sqrt(-p*p*p));
    double m = 2 * sqrt(-p);

    roots[0] = chVar + m*cos(alpha / 3.);
    roots[1] = chVar - m*cos((alpha + M_PI) / 3.);
    roots[2] = chVar - m*cos((alpha - M_PI) / 3.);
  }

  return status;
}

} /* end namespace numerics */
} /* end namespace axom */

#endif // AXOM_NUMERICS_POLY_SOLVE_HPP_

// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_EVAL_INTEGRAL_IMPL_HPP_
#define PRIMAL_EVAL_INTEGRAL_IMPL_HPP_

// Axom includes
#include "axom/config.hpp"  // for compile-time configuration options
#include "axom/primal.hpp"

// MFEM includes
#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#else
  #error "Primal's integral evaluation functions require mfem library."
#endif

// C++ includes
#include <cmath>

namespace axom
{
namespace primal
{
namespace detail
{
/*!
 * \brief Computes the tangent vector to a Bezier curve.
 * 
 * Uses the de Casteljau algorithm. Functionality will eventually
 * be supplanted by Bezier curve class method.
 * 
 * \param [in] c the Bezier curve
 * \param [in] t on [0, 1] the parameter value on which to evaluate the derivative
 * \return the tangent vector
 */
template <typename T>
primal::Vector<double, 2> temp_dt(const primal::BezierCurve<double, 2>& c, T t)
{
  Vector<double, 2> val;

  const int ord = c.getOrder();
  std::vector<T> dCarray(ord + 1);

  // Run de Casteljau algorithm on each dimension
  for(int i = 0; i < 2; ++i)
  {
    for(int p = 0; p <= ord; ++p)
    {
      dCarray[p] = c[p][i];
    }

    // stop one step early and take difference of last two values
    for(int p = 1; p <= ord - 1; ++p)
    {
      const int end = ord - p;
      for(int k = 0; k <= end; ++k)
      {
        dCarray[k] = (1 - t) * dCarray[k] + t * dCarray[k + 1];
      }
    }
    val[i] = ord * (dCarray[1] - dCarray[0]);
  }

  return val;
}

/*!
 * \brief Evaluate a scalar field line integral on a single Bezier curve.
 *
 * Evaluate the scalar field line integral with MFEM integration rule
 *
 * \param [in] c the Bezier curve object
 * \param [in] integrand the lambda function representing the integrand. 
 * Must accept a 2D point as input and return a double
 * \param [in] quad the mfem integration rule containing nodes and weights
 * \return the value of the integral
 */
inline double evaluate_line_integral_component(
  const primal::BezierCurve<double, 2>& c,
  std::function<double(Point2D)> integrand,
  const mfem::IntegrationRule* quad)
{
  // Store intermediate values
  Point2D x_q;
  Vector2D dx_q;

  // Store/compute quadrature result
  double full_quadrature = 0.0;
  for(int q = 0; q < quad->GetNPoints(); q++)
  {
    // Get intermediate quadrature point
    //  at which to evaluate tangent vector
    x_q = c.evaluate(quad->IntPoint(q).x);
    dx_q = detail::temp_dt(c, quad->IntPoint(q).x);

    full_quadrature += quad->IntPoint(q).weight * integrand(x_q) * dx_q.norm();
  }

  return full_quadrature;
}

/*!
 * \brief Evaluate a vector field line integral on a single Bezier curve.
 *
 * Evaluate the vector field line integral with MFEM integration rule
 *
 * \param [in] c the Bezier curve object
 * \param [in] integrand the lambda function representing the integrand. 
 * Must accept a 2D point as input and return a 2D axom vector object
 * \param [in] quad the mfem integration rule containing nodes and weights
 * \return the value of the integral
 */
inline double evaluate_line_integral_component(
  const primal::BezierCurve<double, 2>& c,
  std::function<Vector2D(Point2D)> vec_field,
  const mfem::IntegrationRule* quad)
{
  // Store intermediate values.
  Point2D x_q;
  Vector2D dx_q, func_val;

  // Store/compute quadrature result
  double full_quadrature = 0.0;
  for(int q = 0; q < quad->GetNPoints(); q++)
  {
    // Get intermediate quadrature point
    //  on which to evaluate dot product
    x_q = c.evaluate(quad->IntPoint(q).x);
    dx_q = detail::temp_dt(c, quad->IntPoint(q).x);
    func_val = vec_field(x_q);

    full_quadrature +=
      quad->IntPoint(q).weight * Vector2D::dot_product(func_val, dx_q);
  }

  return full_quadrature;
}

/*!
 * \brief Evaluate the integral across a single component of the curved polygon.
 *
 * Uses a Spectral Mesh-Free Quadrature derived from Green's theorem, evaluating
 * the area integral as a line integral of the antiderivative over the curve.
 *
 * \param [in] cs the array of Bezier curve objects that bound the region
 * \param [in] integrand the lambda function representing the integrand. 
 * Must accept a 2D point as input and return a double
 * \param [in] The lower bound of integration for the antiderivatives
 * \param [in] quad_Q the quadrature rule for the line integral
 * \param [in] quad_P the quadrature rule for the antiderivative
 * \return the value of the integral, which is mathematically meaningless.
 */
template <class Lambda>
double evaluate_area_integral_component(const primal::BezierCurve<double, 2>& c,
                                        Lambda&& integrand,
                                        double int_lb,
                                        const mfem::IntegrationRule* quad_Q,
                                        const mfem::IntegrationRule* quad_P)
{
  // Store some intermediate values
  Point2D x_q, x_qxi;
  double antiderivative = 0.0;

  // Store/compute quadrature result
  double full_quadrature = 0.0;
  for(int q = 0; q < quad_Q->GetNPoints(); q++)
  {
    // Get intermediate quadrature point
    //  on which to evaluate antiderivative
    x_q = c.evaluate(quad_Q->IntPoint(q).x);

    // Evaluate the antiderivative at x_q, add it to full quadrature
    for(int xi = 0; xi < quad_P->GetNPoints(); xi++)
    {
      // Define interior quadrature points
      x_qxi[0] = x_q[0];
      x_qxi[1] = (x_q[1] - int_lb) * quad_P->IntPoint(xi).x + int_lb;

      antiderivative =
        quad_P->IntPoint(xi).weight * (x_q[1] - int_lb) * integrand(x_qxi);

      full_quadrature += quad_Q->IntPoint(q).weight *
        detail::temp_dt(c, quad_Q->IntPoint(q).x)[0] * -antiderivative;
    }
  }

  return full_quadrature;
}

}  // end namespace detail
}  // end namespace primal
}  // end namespace axom

#endif
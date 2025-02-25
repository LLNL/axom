// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internals.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_ELVIRA_ALGORITHM_HPP_
#define AXOM_MIR_ELVIRA_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

// Include these directly for now.
#include "axom/mir/MIRAlgorithm.hpp"
#include "axom/mir/utilities/ExtractZones.hpp"
#include "axom/mir/utilities/MakeZoneCenters.hpp"
#include "axom/mir/utilities/MergeMeshes.hpp"
#include "axom/mir/utilities/NodeToZoneRelationBuilder.hpp"
#include "axom/mir/utilities/RecenterField.hpp"
#include "axom/mir/utilities/ZoneListBuilder.hpp"
#include "axom/mir/views/dispatch_coordset.hpp"
#include "axom/mir/views/MaterialView.hpp"

#include <conduit/conduit.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

#include <algorithm>
#include <string>

// Uncomment to save inputs and outputs.
// #define AXOM_ELVIRA_DEBUG

#if defined(AXOM_ELVIRA_DEBUG)
  #include <conduit/conduit_relay_io_blueprint.hpp>
#endif

namespace axom
{
namespace mir
{
namespace elvira
{

enum class Direction
{
  VERTICAL = 0,
  HORIZONTAL = 1
};

struct Result2D
{
  int plane {2};
  int difference_used {0};
  double columns[3] {0., 0., 0.};
  double normal[3][2] {{0., 0.}, {0., 0.}, {0., 0.}};
};

/*!
 * \brief 2nd order 2-d mir function: returns the normal to the slope in 
 *        the direction away from the material.
 *
 * \param[out] result The result object.
 * \param vf   An array containing volume fractions for the current material
 *             arranged logically like this, according to stencil indices:
 *                       0 - lower left neighbor
 *                       1 - lower neighbor
 *                       2 - lower right neighbor
 *                       3 - left neighbor
 *                       4 - zone in question
 *                       5 - right neighbor
 *                       6 - upper left neighbor
 *                       7 - upper neighbor
 *                       8 - upper right neighbor
 *             Note that the actual vf data can be 9 or 27 elements depending
 *             on the dimensionality.
 *
 * \param ivf The indices to use to pull data out of vf.
 * \param direction hv direction to use
 *
 * \return The result object will contain these values:
 *         normal  : normal (x,y) pointing away from the material
 *         err     : the L2 error of the answer returned
 *
 * \note Adapted from J. Grandy's Overlink code
 */
AXOM_HOST_DEVICE
void elvira2d(Result2D &result, const double *vf, const int *ivf, Direction direction)
{
  const double jb = vf[ivf[0]] + vf[ivf[1]] + vf[ivf[2]];  // bottom row
  const double jm = vf[ivf[3]] + vf[ivf[4]] + vf[ivf[5]];  // middle row
  const double jt = vf[ivf[6]] + vf[ivf[7]] + vf[ivf[8]];  // top row

  const double il = vf[ivf[0]] + vf[ivf[3]] + vf[ivf[6]];  // left column
  const double im = vf[ivf[1]] + vf[ivf[4]] + vf[ivf[7]];  // middle column
  const double ir = vf[ivf[2]] + vf[ivf[5]] + vf[ivf[8]];  // right column

  double slope[3] {0., 0., 0.};  // backward, cent, and forward differencing
  double slopev = 0.;            // Sets overall sign of normal.
  int ihv;
  if(direction == Direction::VERTICAL)
  {
    slope[0] = (im - il);        // backward differencing
    slope[1] = (ir - il) / 2.0;  // central differencing
    slope[2] = (ir - im);        // forward differencing
    slopev = (jb - jt);

    result.columns[0] = il;
    result.columns[1] = im;
    result.columns[2] = ir;

    ihv = 0;
  }
  else
  {
    slope[0] = (jm - jb);        // backward differencing
    slope[1] = (jt - jb) / 2.0;  // central differencing
    slope[2] = (jt - jm);        // forward differencing
    slopev = (il - ir);

    result.columns[0] = jb;
    result.columns[1] = jm;
    result.columns[2] = jt;

    ihv = 1;
  }

  const double sign = (slopev > 0.0) ? (1.) : (-1.);

  // return the normal in the direction away from the material
  // i = loop over three difference schemes.
  for(int i = 0; i < 3; i++)
  {
    result.normal[i][1 - ihv] = sign;   // steep direction
    result.normal[i][ihv] = -slope[i];  // shallow
  }
}

/*!
 * \brief Get chi squared for actual and fitted vol fracs.
 *
 * \param vf actual volume fractions.
 * \param vfs trial volume fractions.
 * \param ivf index list of volume fractions to be compared.
 * \param k The number of volume fractions to compare.
 *
 * \return chi squared (measure of error)
 *
 * \note Adapted from J. Grandy's Overlink code
 */
AXOM_HOST_DEVICE
inline double elvira_chisq(const double *vf,
                           const double *vfs,
                           const int *ivf,
                           int k)
{
  double chisq = 0.0;
  for(; k--;)
  {
    const double c = vf[ivf[k]] - vfs[ivf[k]];
    chisq += c * c;
  }
  return chisq;
}

/*!
 * \brief Normalize 2d vector, append zero Z component.
 *
 * Notes: Assumes n2 is not (0,0) .
 *
 * \note Adapted from J. Grandy's Overlink code
 */
AXOM_HOST_DEVICE
inline void norm2d(double n2[2], double n3[3])
{
  constexpr double PTINY = 1.e-15;
  // Compute 2D vector magnitude.
  double magnitude = sqrt(n2[0] * n2[0] + n2[1] * n2[1]);
  double magnitude2 = magnitude * magnitude;
  double magnitudeReciprocal = magnitude / (magnitude2 + PTINY);
  // normalize vector, make it 3D
  n3[0] = n2[0] * magnitudeReciprocal;
  n3[1] = n2[1] * magnitudeReciprocal;
  n3[2] = (magnitude2 <= 0.) ? 1.0 : 0.0;
}

/*!
 * \brief Return sorted, positive components of vector length 3.
 */
AXOM_HOST_DEVICE
inline void n_sort(const double n[3],   // input vector
                   double na[3])  // output vector
{
  na[0] = axom::utilities::abs(n[0]);
  na[1] = axom::utilities::abs(n[1]);
  na[2] = axom::utilities::abs(n[2]);

  /* Sort the components. */
  if(na[0] > na[2])
  {
    axom::utilities::swap(na[0], na[2]);
  }
  if(na[1] > na[2])
  {
    axom::utilities::swap(na[1], na[2]);
  }
  if(na[0] > na[1])
  {
    axom::utilities::swap(na[0], na[1]);
  }
}

/*!
 *
 * \brief Compute volume fraction of unit cube at origin below plane.
 *
 * \param d  Plane displacement
 * \param n1 Sorted positive normal components of plane 1
 * \param n2 Sorted positive normal components of plane 2
 * \param n3 Sorted positive normal components of plane 3 (highest)
 *
 * \return Volume fraction
 *
 * \note Plane equation is (n,x) + d = 0.
 *       Cross section formulas from Youngs (1987).
 */
AXOM_HOST_DEVICE
double vf_1cube(double d,   // plane displacement
                double n1,  // sorted positive normal components of plane
                double n2,  // sorted positive normal components of plane
                double n3)  // sorted positive normal components of plane
{
  double vf = 0.0, dd1, dd2, dd3;
  enum
  {
    cut_triangle,
    cut_quada,
    cut_pent,
    cut_hex,
    cut_quadb
  } xsec;
  constexpr double one6 = 1. / 6.;
  double nbd = n1 + n2 + n3;

  int ifdhi = 0;
  if(d > 0.5 * nbd)
  {
    d = nbd - d;
    ifdhi = 1; /* symmetric at vf = 0.5 */
  }
  if(d <= 0.0)
  {
    vf = 0.0;
  } /* For numerical safety.   */
  else
  {
    /* Do the vf calculation per Youngs.   */
    /* We have already guaranteed d>0. */

    /* Select cross section type.   */
    if(d <= n1)
    {
      xsec = cut_triangle;
    }
    else if(d <= n2)
    {
      xsec = cut_quada;
    }
    else if(d <= n3)
    {
      if(d <= (n1 + n2))
      {
        xsec = cut_pent;
      }
      else
      {
        xsec = cut_quadb;
      }
    }
    else
    {
      xsec = cut_hex;
    }

    /* Compute volume fraction based on type.   */
    switch(xsec)
    {
    case cut_triangle:

      vf = one6 * d * d * d / (n1 * n2 * n3);
      break;

    case cut_quada:

      vf = one6 * (3.0 * d * (d - n1) + n1 * n1) / (n2 * n3);
      break;

    case cut_pent:

      dd1 = d - n1;
      dd2 = d - n2;
      vf =
        one6 * (d * d * d - dd1 * dd1 * dd1 - dd2 * dd2 * dd2) / (n1 * n2 * n3);
      break;

    case cut_hex:

      dd1 = d - n1;
      dd2 = d - n2;
      dd3 = d - n3;
      vf = one6 *
        (d * d * d - dd1 * dd1 * dd1 - dd2 * dd2 * dd2 - dd3 * dd3 * dd3) /
        (n1 * n2 * n3);
      break;

    case cut_quadb:

      vf = (d - 0.5 * (n1 + n2)) / n3;
      break;
    default:
      break;
    }
  }
  if(ifdhi)
  {
    vf = 1.0 - vf;
  } /* vf=0.5 symmetry.   */
  constexpr double eps = 1.0e-15;
  constexpr double lower = eps;
  constexpr double upper = 1. - eps;
  if(vf < lower) vf = 0.0; /* Truncate. */
  if(vf > upper) vf = 1.0; /* Truncate. */
  return (vf);
}

/*!
 * \brief Solve a CUBic from 4 equally spaced Points.
 *
 * \param x An array containing 4 x values.[0] and [3] are interval endpoints, [1],[2] unused.
 * \param y An array containing 4 y values.Dependent variable at 4 equally spaced points
 *
 * \return A root of the cubic. 
 *
 * \note  
 *   x[0], x[3] are interval endpoints, with function values y[0], y[3].
 *   y[1], y[2] are function at 2/3x[0] + 1/3x[3] and 1/3x[0] + 2/3x[3].
 *
 */
AXOM_HOST_DEVICE
double cub4p(const double *x, const double *y)
{
  double dstar;
  double e0 = 0.0, e1 = 0.0, e2 = 0.0, ep, em, ea, eb;
  double de0, dde0, de0sq, dde0sq;
  double de1, dde1, de1sq, dde1sq, demin, ddemax;
  double y0, y1, y2, y3;
  double w0 = 0.0, w1 = 0.0, w2;
  double b, y30, y12, a, tol2;
  double tolsec = 1.0e-28;
  double one6 = 0.16666666666666667;
  int ifsec, ifbis;

  y0 = y[0];
  y1 = y[1];
  y2 = y[2];
  y3 = y[3];

  tol2 = (y0 * y0 + y1 * y1 + y2 * y2 + y3 * y3) * tolsec;

  /* Initialize secant search.   */
  /* The default covers end cases where vks never 
         changes sign, as can occur for vol fraction 
         very near zero or one.   */
  ifbis = 1;
  if(y2 * y3 < 0.0)
  {
    e0 = one6;
    w0 = y2;
    e1 = 0.5;
    w1 = y3;
  }
  else if(y1 * y2 < 0.0)
  {
    e0 = -one6;
    w0 = y1;
    e1 = one6;
    w1 = y2;
  }
  else if(y0 * y1 < 0.0)
  {
    e0 = -0.5;
    w0 = y0;
    e1 = -one6;
    w1 = y1;
  }

  // The quantities y1, y2, y3 have already been squared.  These
  // test for not equal to 0 (clears compiler warning).
  else if(y3 * y3 <= 0.0)
  {
    e2 = 0.5;
    ifbis = 0;
  }
  else if(y2 * y2 <= 0.0)
  {
    e2 = one6;
    ifbis = 0;
  }
  else if(y1 * y1 <= 0.0)
  {
    e2 = -one6;
    ifbis = 0;
  }
  else
  {
    e2 = -0.5;
    ifbis = 0;
  }

  ifsec = 0;

  /* Initialize with one Newton iteration.   */
  if(ifbis)
  {
    /* Precompute these. */
    y30 = y3 - y0;
    y12 = y1 - y2;
    b = 4.5 * y30 + 13.5 * y12;

    ep = e0 + 0.5;
    em = e0 - 0.5;
    ea = e0 + one6;
    eb = e0 - one6;
    a = 4.5 * (ep * y3 - em * y0) + 13.5 * (eb * y1 - ea * y2);

    de0 = y30 + (ep + em) * a + ep * em * b;
    dde0 = 2.0 * (a + (ep + em) * b);
    de0sq = de0 * de0;
    dde0sq = dde0 * dde0;

    ep = e1 + 0.5;
    em = e1 - 0.5;
    ea = e1 + one6;
    eb = e1 - one6;
    a = 4.5 * (ep * y3 - em * y0) + 13.5 * (eb * y1 - ea * y2);

    de1 = y30 + (ep + em) * a + ep * em * b;
    dde1 = 2.0 * (a + (ep + em) * b);
    de1sq = de1 * de1;
    dde1sq = dde1 * dde1;

    /* This division is safe. */
    e2 = e0 + (e1 - e0) * w0 / (w0 - w1);
  }
  else
  {
    // This clears up compiler warnings for uninitialized vars.
    e0 = e1 = w0 = w1 = 0.0;
    de0 = de0sq = dde0sq = de1 = de1sq = dde1sq = b = y30 = 0.0;
  }

  /* Do the bisection until ok for secant.   */
  while(ifbis)
  {
    ep = e2 + 0.5;
    em = e2 - 0.5;
    ea = e2 + one6;
    eb = e2 - one6;

    a = 4.5 * (ep * y3 - em * y0) + 13.5 * (eb * y1 - ea * y2);
    w2 = ep * y3 - em * y0 + ep * em * a;

    /* If we have converged, quit. */
    if(w2 * w2 < tol2)
    {
      ifbis = ifsec = 0;
    }

    else
    {
      /* Else, do the bisection, rotate. */
      /* Set correct interval to (e0,e1). */
      if(w1 * w2 < 0.0)
      {
        e0 = e1;
        w0 = w1;
        de0 = de1;
        dde0 = dde1;
        de0sq = de1sq;
        dde0sq = dde1sq;
      }

      e1 = e2;
      w1 = w2;

      de1 = y30 + (ep + em) * a + ep * em * b;
      dde1 = 2.0 * (a + (ep + em) * b);
      de1sq = de1 * de1;
      dde1sq = dde1 * dde1;

      /* Test for secant applicability.   */
      if(de0 * de1 > 0.0)
      { /* dy/de same sign. */

        demin = (de0sq < de1sq) ? de0sq : de1sq;
        ddemax = (dde0sq > dde1sq) ? dde0sq : dde1sq;

        /* Test on second derivative.   */
        /* This works since equation is cubic.  */
        if(demin > (e1 - e0) * (e1 - e0) * ddemax)
        {
          ifbis = 0;
          ifsec = 1;
        }
      }

      /* Next estimate.    */
      if(ifsec)
        e2 = e0 + (e1 - e0) * w0 / (w0 - w1);
      else
        e2 = 0.5 * (e0 + e1);
    }
  }

  /* Refine with secant search.   */
  while(ifsec)
  {
    ep = e2 + 0.5;
    em = e2 - 0.5;
    ea = e2 + one6;
    eb = e2 - one6;

    a = 4.5 * (ep * y3 - em * y0) + 13.5 * (eb * y1 - ea * y2);

    /* Evaluate the volume at e2 as w2. */
    w2 = ep * y3 - em * y0 + ep * em * a;

    if(w2 * w2 < tol2)
    {
      ifsec = 0;
    }
    else
    {
      e0 = e1;
      e1 = e2;
      w0 = w1;
      w1 = w2;
      e2 = e0 + (e1 - e0) * w0 / (w0 - w1);
    }
  }

  dstar = x[0] + (e2 + 0.5) * (x[3] - x[0]);

  return (dstar);
}

/*!
 *
 * \brief Calculate displacement of plane through unit cube located 
 *        between (1,1,1) and (2,2,2) given volume fraction.  
 *
 * \param n    Plane normal
 * \param vf13 Volume fraction
 *
 * \return Plane displacement, d, such that (n,x)+d=0 defines plane.
 *
 * \note "Lowest d" node depends on signs of n components, hence additional
 *         offsets from the lowest d node to origin.  
 *
 * \note Adapted from J. Grandy's Overlink code
 */
AXOM_HOST_DEVICE
double d_3cube(const double n[3], double vf13)
{
  int i, ik, ifdhi;
  double d[5], vc[5], v[4], x[4];
  double na[3], n1, n2, n3, nbd, dstar;
  double one3 = 0.3333333333333333;
  double two3 = 2.0 * one3;

  /* Get sorted positive normal components.   */
  n_sort(n, na);
  n1 = na[0];
  n2 = na[1];
  n3 = na[2];
  nbd = n1 + n2 + n3;

  /* Remove cases of bad volume fractions. */
  if(vf13 <= 0.0)
  {
    dstar = 0.0;
  }
  else if(vf13 >= 1.0)
  {
    dstar = nbd;
  }
  else
  {
    /* Set up ascending sequence of d's, where plane passes through 
      cube nodes.    */
    d[0] = 0.0;
    d[1] = n1;
    d[2] = n2;
    if((n1 + n2) < n3)
    {
      d[3] = n1 + n2;
    }
    else
    {
      d[3] = n3;
    }
    d[4] = 0.5 * nbd;

    vc[0] = 0.0;
    vc[4] = 0.5;

    /* Symmetrize: For volume fractions > 0.5, use complements.   */
    if(vf13 > 0.5)
    {
      ifdhi = 1;
      vc[0] = 1.0;
      for(i = 0; i < 4; i++)
      {
        d[i] = nbd - d[i];
      }
    }
    else
    {
      ifdhi = 0;
    }

    /* Get vol fraction for each d, where plane passes through node. */
    for(i = 1; i <= 3; i++)
    {
      vc[i] = vf_1cube(d[i], n1, n2, n3);
    }

    /* Find volume fraction's interval.   */
    for(ik = 0; ik < 4; ik++)
    {
      if((vc[ik] - vf13) * (vc[ik + 1] - vf13) <= 0)
      {
        break;
      }
    }

    if(ik == 4) return (d[4]); /* Exact boundary. */

    /* Cross section is a simple tetrahedron.   */
    if(ik == 0)
    {
      if(ifdhi)
      {
        dstar = d[0] - pow(6.0 * (1.0 - vf13) * n1 * n2 * n3, one3);
      }
      else
      {
        dstar = pow(6.0 * vf13 * n1 * n2 * n3, one3);
      }
    }

    /* More complicated: set up, solve the cubic.   */
    else
    {
      x[0] = d[ik];
      x[3] = d[ik + 1];
      v[0] = vc[ik] - vf13;
      v[3] = vc[ik + 1] - vf13;
      x[1] = two3 * x[0] + one3 * x[3];
      v[1] = vf_1cube(x[1], n1, n2, n3) - vf13;
      x[2] = one3 * x[0] + two3 * x[3];
      v[2] = vf_1cube(x[2], n1, n2, n3) - vf13;

      dstar = cub4p(x, v);
    }

    if(dstar < 0.0) dstar = 0.0;
    if(dstar > nbd) dstar = nbd;
  }

  dstar += (n[0] + n[1] + n[2]); /* (1,1,1) offset: center zone. */

  /* Additional offset for lowest d node of cube. */
  for(i = 0; i < 3; i++)
    if(n[i] < 0.0)
    {
      dstar += n[i];
    }

  dstar = -dstar; /* conventional (n,x)+d=0.   */
  return (dstar);
}

/*!
 * \brief Given a plane through 3 cubed stencil, find volume fraction 
 *        inside half space away from which normal points ("below" normal)
 *        for the zones that are center, face, edges. 
 *
 * \param n The normal to use.
 * \param pd The displacement from the lower left corner of the zone.
 * \param[out] vfs An array of volume fractions for the ivf zone ids.
 * \param ivf A list of stencil indices.
 * \param k The number of stencil indices in \a ivf.
 *
 * \note  Plane equation is (n,x) + d = 0.
 *
 * \note Adapted from J. Grandy's Overlink code
 */
AXOM_HOST_DEVICE
void vf_3cube(double n[3],  // plane normal (pointing outward).
              double pd,    // -displacement.
              double *vfs,
              const int *ivf,
              int k)
{
  /* find node of lowest d */
  /* as 0 or 1 offset for each coord.  */
  const double delx0 = (n[0] < 0.0) ? 1.0 : 0.0;
  const double dely0 = (n[1] < 0.0) ? 1.0 : 0.0;
  const double delz0 = (n[2] < 0.0) ? 1.0 : 0.0;

  double na[3];
  n_sort(n, na);

  for(int i = 0; i < k; i++)
  {
    // Compute lower left coordinate of stencil box.
    // The values are 0,1,2
    const double sx = static_cast<double>(ivf[i] % 3);
    const double sy = static_cast<double>((ivf[i] % 9) / 3);
    const double sz = static_cast<double>(ivf[i] / 9);

    const double x0 = sx + delx0;
    const double y0 = sy + dely0;
    const double z0 = sz + delz0;
    const double d = -(pd + (n[0] * x0 + n[1] * y0 + n[2] * z0));

    const double vf = vf_1cube(d, na[0], na[1], na[2]);

    vfs[ivf[i]] = vf;
  }
}

/*!
 * \brief 2nd order 2-d mir function.
 *
 * \param vf The input volume fraction stencil. This is 9 elements.
 * \param[out] n The output normal (x,y,0) pointing away from the material.
 *
 * \note Adapted from J. Grandy's Overlink code
 */
AXOM_HOST_DEVICE
void elvira2xy(const double *vf, double n[3])
{
  // These are indices into the volume fractions that pull out values in the XY plane.
  const int ivf[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  // Pick some positions out of ivf.
  constexpr int center = 4;
  constexpr int bottom = 1;
  constexpr int top = 7;
  constexpr int left = 3;
  constexpr int right = 5;

  // Pick a direction.
  double variance[2];
  variance[0] = ((vf[center] - vf[bottom]) * (vf[center] - vf[bottom]) +
                 (vf[center] - vf[top]) * (vf[center] - vf[top]));
  variance[1] = ((vf[center] - vf[left]) * (vf[center] - vf[left]) +
                 (vf[center] - vf[right]) * (vf[center] - vf[right]));
  const auto direction =
    (variance[0] < variance[1]) ? Direction::HORIZONTAL : Direction::VERTICAL;

  // Compute normals
  Result2D result;
  elvira2d(result, vf, ivf, direction);
  result.plane = 2;

  double chmin = 0.0;
  double vfs[9];
  const double chfac[3] = {1.0, 0.25, 1.0};

  // Choose difference schemes.
  for(int diff = 0; diff < 3; diff++)
  {
    // Turn the normal into normalized 3D vector.
    double n3[3];
    norm2d(result.normal[diff], n3);

    // Compute displacement for center stencil for this normal.
    double d = d_3cube(n3, vf[center]);

    // Compute the vfs for the various boxes in the stencil using current n3.
    vf_3cube(n3, d, vfs, ivf, 9);

    // Compute diffs between computed, actual vf's.
    double chisq = elvira_chisq(vf, vfs, ivf, 9);
    // Prefer central difference (since a low value will lower the error)
    chisq *= chfac[diff];

    // Keep the normal with the least error
    if((diff == 0) || (chisq < chmin))
    {
      chmin = chisq;
      n[0] = n3[0];
      n[1] = n3[1];
      n[2] = n3[2];
    }
  }
}

/*!
 * \brief Multiply the normal times the Jacobian and normalize.
 *
 * \param normal The normal to transform, nx,ny,nz
 * \param jac    The jacobian to multiply by
 *
 * \note Adapted from J. Grandy's Overlink code
 */
AXOM_HOST_DEVICE
void transform(double *normal, const double jac[3][3])
{
  SLIC_ASSERT(normal != nullptr);

  double norm = 0.0;
  double dfvsum = 0.0;  // stores nx^2 + ny^2 + nz^2
  double dfv[3], delfv[3];

  dfv[0] = normal[2];  // BJW: Why switch components? Rotation?
  dfv[1] = normal[1];
  dfv[2] = normal[0];

  for(int f = 0; f < 3; f++)
  {
    delfv[f] = dfv[0] * jac[0][f] + dfv[1] * jac[1][f] + dfv[2] * jac[2][f];

    norm += delfv[f] * delfv[f];
    dfvsum += dfv[f] * dfv[f];
  }

  // dfvsum tests small variations in vol frac ( see mira1_youngs )
  if((dfvsum > 1.0e-24) && (norm > 0.0))
  {
    norm = 1.0 / sqrt(norm);
    for(int f = 0; f < 3; f++)
    {
      delfv[f] = delfv[f] * norm;
    }
  }
  else
  {
    delfv[0] = 1.0; /* An arbitrary slope. */
    delfv[1] = 0.0;
    delfv[2] = 0.0;
  }

  /* Store the slope.  */
  normal[0] = delfv[0];
  normal[1] = delfv[1];
  normal[2] = delfv[2];
}

/*!
 * \brief Compute the jacobian from an input stencil of x,y,z coordinates derived
 *        from zone centroids.
 *
 * \param xcst A stencil of X coordinate values. 9 values if ndims == 2, 27 values if ndims == 3.
 * \param ycst A stencil of Y coordinate values. 9 values if ndims == 2, 27 values if ndims == 3.
 * \param zcst A stencil of Z coordinate values. 9 values if ndims == 2, 27 values if ndims == 3.
 * \param ndims The number of dimensions (2 or 3).
 * \param[out] jac A 3x3 Jacobian matrix.
 *
 * \note  Adapted from J. Grandy's BasicStencil.cc:52 
 */
AXOM_HOST_DEVICE
void computeJacobian(const double *xcst,
                     const double *ycst,
                     const double *zcst,
                     int ndims,
                     double jac[3][3])
{
  double del[3][3], det;
  int f, f0, f1, f2, g0, g1, g2;
  int perm1[3] = {1, 2, 0};
  int perm2[3] = {2, 0, 1};

  int idx_2D[6] = {4, 4, 1, 7, 3, 5};
  int idx_3D[6] = {4, 22, 10, 16, 12, 14};

  int *idx = (ndims == 3) ? idx_3D : idx_2D;

  /*
  *  Note face convention : opposite face pairs (01), (23), (45) form
  *  a left handed orientation so that det is negative for uniform
  *  positively oriented meshes.
  */
  for(f = 0; f < 3; f++)
  {
    f1 = idx[f * 2];
    f2 = idx[f * 2 + 1];

    del[0][f] = 0.5 * (xcst[f2] - xcst[f1]);
    del[1][f] = 0.5 * (ycst[f2] - ycst[f1]);
    del[2][f] = 0.5 * (zcst[f2] - zcst[f1]);
  }  // END for all face pairs

  /* Construct the Left Jacobian matrix. */
  /* Neglect the determinant since we renormalize. */
  /* This is jac = inverse(del) = 1/det(del) Adj(del).  */
  /* (Note g0 <--> f0 in jac).    */
  for(f0 = 0; f0 < 3; f0++)
  {
    f1 = perm1[f0];
    f2 = perm2[f0];

    for(g0 = 0; g0 < 3; g0++)
    {
      g1 = perm1[g0];
      g2 = perm2[g0];
      jac[g0][f0] = del[f1][g1] * del[f2][g2] - del[f1][g2] * del[f2][g1];
    }
  }

  /* Check the sign of the determinant. */
  det = del[0][0] * jac[0][0] + del[0][1] * jac[1][0] + del[0][2] * jac[2][0];

  /* We still need the sign of 1/det in Cramer's rule so we
    * change the sign here.  We don't need a normalized jacobian since
    * the mix slope normal is normalized later.    */
  if(det < 0.0)
  {
    for(f0 = 0; f0 < 3; f0++)
    {
      for(int f3 = 0; f3 < 3; f3++)
      {
        jac[f0][f3] *= -1.0;
      }
    }

  }  // END if negative determinant

  /* Fix 2D stencil */
  if(ndims == 2)
  {
    jac[2][0] = jac[2][1] = jac[0][2] = jac[1][2] = 0.0;
    jac[2][2] = 1.0;
  }
}

AXOM_HOST_DEVICE
constexpr int getStencilSize(int ndims)
{
  return (ndims == 3) ? 27 : ((ndims == 2) ? 9 : 3);
}

/*!
 * \brief Create normals for the material interface fragments in the selected zone.
 *
 * \param szIndex The index of the zone in the selected, sorted zones list. We probably can remove this...
 * \param matCount The number of materials in the current zone.
 * \param offset An offset for the current zone, used to index into the material and vector fragment data.
 * \param zoneMatStencilView The material volume fraction stencil data for all fragments.
 * \param fragmentVectorsView The normals for all fragments.
 *
 * \note Calling this function will update some vectors in the \a fragmentVectorsView.
 */
AXOM_HOST_DEVICE
void elvira(axom::IndexType szIndex,
            int matCount,
            int offset,
            axom::ArrayView<double> zoneMatStencilView,
            axom::ArrayView<double> fragmentVectorsView,
            int iskip,
            int ndims)
{
  constexpr int StencilSize = getStencilSize(ndims);
  constexpr int numVectorComponents = 3;

  for(int m = 0; m < matCount; m++)
  {
    const auto fragmentIndex = offset + m;

    // Get the pointer for this interface's normal
    double *normal =
      fragmentVectorsView.data() + fragmentIndex * numVectorComponents;

    // Get the pointer for this material's volume fraction stencil.
    double *vol_fracs = zoneMatStencilView.data() + fragmentIndex * StencilSize;

    if(m != iskip)
    {
      if(ndims == 2)
      {
        elvira2xy(vol_fracs, normal);
      }
      else
      {
        // elvira3d(vol_fracs, normal);
      }
    }
    else
    {
      normal[0] = 1.;
      normal[1] = 0.;
      normal[2] = 0.;
    }
  }
}

}  // namespace elvira

/*!
 * \brief Implements Elvira algorithm for 2D structured meshes.
 */
template <typename ExecSpace, typename IndexPolicy, typename CoordsetView, typename MatsetView>
class ElviraAlgorithm2D : public axom::mir::MIRAlgorithm
{
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;

public:
  using TopologyView = axom::mir::views::StructuredTopologyView<IndexPolicy>;
  using ConnectivityType = typename TopologyView::ConnectivityType;

  // Make sure we only instantiate this algorithm for 2D meshes.
  static_assert(IndexPolicy::dimension() == 2, "IndexPolicy must be 2D");

  /*!
   * \brief Constructor
   *
   * \param topoView The topology view to use for the input data.
   * \param coordsetView The coordset view to use for the input data.
   * \param matsetView The matset view to use for the input data.
   */
  ElviraAlgorithm2D(const TopologyView &topoView,
                    const CoordsetView &coordsetView,
                    const MatsetView &matsetView)
    : axom::mir::MIRAlgorithm()
    , m_topologyView(topoView)
    , m_coordsetView(coordsetView)
    , m_matsetView(matsetView)
  { }

  /// Destructor
  virtual ~ElviraAlgorithm2D() = default;

  /*!
   * \brief Perform material interface reconstruction on a single domain.
   *
   * \param[in] n_topo The Conduit node containing the topology that will be used for MIR.
   * \param[in] n_coordset The Conduit node containing the coordset.
   * \param[in] n_fields The Conduit node containing the fields.
   * \param[in] n_matset The Conduit node containing the matset.
   * \param[in] n_options The Conduit node containing the options that help govern MIR execution.
   *
   * \param[out] n_newTopo A node that will contain the new clipped topology.
   * \param[out] n_newCoordset A node that will contain the new coordset for the clipped topology.
   * \param[out] n_newFields A node that will contain the new fields for the clipped topology.
   * \param[out] n_newMatset A Conduit node that will contain the new matset.
   * 
   */
  virtual void executeDomain(const conduit::Node &n_topo,
                             const conduit::Node &n_coordset,
                             const conduit::Node &n_fields,
                             const conduit::Node &n_matset,
                             const conduit::Node &n_options,
                             conduit::Node &n_newTopo,
                             conduit::Node &n_newCoordset,
                             conduit::Node &n_newFields,
                             conduit::Node &n_newMatset) override
  {
    namespace bputils = axom::mir::utilities::blueprint;
    using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
    //const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    AXOM_ANNOTATE_SCOPE("ElviraAlgorithm2D");

    // Copy the options to make sure they are in the right memory space.
    conduit::Node n_options_copy;
    bputils::copy<ExecSpace>(n_options_copy, n_options);
    n_options_copy["topology"] = n_topo.name();

    // Get selected zones from the options.
    bputils::SelectedZones<ExecSpace> selectedZones(
      m_topologyView.numberOfZones(),
      n_options_copy);
    const auto selectedZonesView = selectedZones.view();

    // Partition the selected zones into clean, mixed lists.
    axom::Array<axom::IndexType> cleanZones, mixedZones;
    bputils::ZoneListBuilder<ExecSpace, TopologyView, MatsetView> zlb(
      m_topologyView,
      m_coordsetView);
    zlb.execute(selectedZonesView, cleanZones, mixedZones);

    if(cleanZones.size() > 0 && mixedZones.size() > 0)
    {
      // Some clean, some mixed.
    }
    else if(mixedZones.size() > 0)
    {
      // All mixed.
    }
    else if(cleanZones.size() > 0)
    {
      // All clean. We could just pass the inputs through.
      n_newTopo.set_external(n_topo);
      n_newCoordset.set_external(n_coordset);
      n_newFields.set_external(n_fields);
      n_newMatset.set_external(n_matset);
    }
  }

  void processMixedZones(const axom::ArrayView<axom::IndexType> mixedZonesView,
                         const conduit::Node &n_topo,
                         const conduit::Node &n_coordset,
                         const conduit::Node &n_fields,
                         const conduit::Node &n_matset,
                         const conduit::Node &n_options,
                         conduit::Node &n_newTopo,
                         conduit::Node &n_newCoordset,
                         conduit::Node &n_newFields,
                         conduit::Node &n_newMatset)
  {
    AXOM_ANNOTATE_SCOPE("processMixedZones");
    namespace bputils = axom::mir::utilities::blueprint;
    using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
    using reduce_policy =
      typename axom::execution_space<ExecSpace>::reduce_policy;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    constexpr int ndims = TopologyView::dimension();

    // Count the number of fragments we'll make for the mixed zones.
    // We also
    AXOM_ANNOTATE_BEGIN("counting");
    RAJA::ReduceSum<reduce_policy, axom::IndexType> num_reduce(0);

    const auto nzones = mixedZonesView.size();
    axom::Array<axom::IndexType> matCount(nzones, nzones, allocatorID);
    axom::Array<axom::IndexType> matZone(nzones, nzones, allocatorID);
    auto matCountView = matCount.view();
    auto matZoneView = matZone.view();

    const auto matsetView = m_matsetView;
    axom::for_all<ExecSpace>(
      mixedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = mixedZonesView[szIndex];
        const auto nmats = matsetView.numberOfMaterials(zoneIndex);
        matCountView[szIndex] = nmats;
        matZoneView[szIndex] = zoneIndex;
        num_reduce += nmats;
      });
    const auto numFragments = num_reduce.get();
    AXOM_ANNOTATE_END("counting");

    //--------------------------------------------------------------------------
    // Sort the zones by the mat count. This should make adjacent zones in the
    // list more likely to have the same number of materials.
    AXOM_ANNOTATE_BEGIN("sorting");
    RAJA::sort_pairs<loop_policy>(RAJA::make_span(matCountView.data(), nzones),
                                  RAJA::make_span(matZoneView.data(), nzones));
    AXOM_ANNOTATE_END("sorting");

    //--------------------------------------------------------------------------
    AXOM_ANNOTATE_BEGIN("offsets");
    axom::Array<axom::IndexType> matOffset(nzones, nzones, allocatorID);
    auto matOffsetView = matOffset.view();
    axom::exclusive_scan<ExecSpace>(matCountView, matOffsetView);
    AXOM_ANNOTATE_END("offsets");

    //--------------------------------------------------------------------------
    // NOTE: this makes zone centers for all zones.
    AXOM_ANNOTATE_BEGIN("centroids");
    bputils::MakeZoneCenters<ExecSpace, TopologyView, CoordsetView> zc(
      m_topologyView,
      m_coordsetView);
    conduit::Node n_zcfield;
    zc.execute(n_topo, n_coordset, n_zcfield);
    using zc_value = typename CoordsetView::value_type;
    axom::ArrayView<zc_value> xview, yview, zview;
    xview = bputils::make_array_view<zc_value>(n_zcfield["values/x"]);
    yview = bputils::make_array_view<zc_value>(n_zcfield["values/y"]);
    if(n_zcfield.has_path("values/z"))
    {
      zview = bputils::make_array_view<zc_value>(n_zcfield["values/z"]);
    }
    AXOM_ANNOTATE_END("centroids");

    //--------------------------------------------------------------------------
    // Retrieve stencil data for each zone material.
    AXOM_ANNOTATE_BEGIN("stencil");
    constexpr int StencilSize = elvira::getStencilSize(ndims);
    const auto numFragmentsStencil = numFragments * StencilSize;

    axom::Array<double> zoneMatStencil(numFragmentsStencil,
                                       numFragmentsStencil,
                                       allocatorID);
    auto zoneMatStencilView = zoneMatStencil.view();

    // Sorted material ids for each zone.
    axom::Array<typename MatsetView::IndexType> sortedMaterialIds(numFragments,
                                                                  numFragments,
                                                                  allocatorID);
    auto sortedMaterialIdsView = sortedMaterialIds.view();

    // Coordinate stencil data for each zone.
    const auto nzonesStencil = nzones * StencilSize;
    axom::Array<double> xcStencil(nzonesStencil, nzonesStencil, allocatorID);
    axom::Array<double> ycStencil(nzonesStencil, nzonesStencil, allocatorID);
    axom::Array<double> zcStencil(nzonesStencil, nzonesStencil, allocatorID);
    auto xcStencilView = xcStencil.view();
    auto ycStencilView = ycStencil.view();
    auto zcStencilView = zcStencil.view();

    // Traverse the selected zones based on how many materials there are in a zone.
    TopologyView deviceTopologyView(m_topologyView);
    axom::for_all<ExecSpace>(
      matZoneView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        // The selected zone index in the whole mesh.
        const auto zoneIndex = matZoneView[szIndex];
        const auto matCount = matCountView[szIndex];
        // Where to begin writing this zone's fragment data.
        const auto offset = matOffsetView[szIndex];

        // Get materials for this zone from the matset.
        typename MatsetView::IDList ids;
        typename MatsetView::VFList vfs;
        matsetView.zoneMaterials(zoneIndex, ids, vfs);

        // Sort the materials by the volume fraction. Save sorted ids in sortedMaterialIdsView.
        axom::utilities::sort_multiple(vfs.data(), ids.data());
        for(axom::IndexType m = 0; m < matCount; m++)
        {
          const auto fragmentIndex = offset + m;
          sortedMaterialIdsView[fragmentIndex] = ids[m];
        }

        // Retrieve the stencil data from neighbor zones.
        auto logical = deviceTopologyView.indexing().IndexToLogicalIndex(zoneIndex);
        for(int si = 0; si < StencilSize; si++)
        {
          // Stencil neighbor logical index.
          typename TopologyView::LogicalIndex neighbor(logical);
          const int disp[3] = {(si % 3) - 1, ((si % 9) / 3) - 1, (si / 9) - 1};
          for(int d = 0; d < ndims; d++)
          {
            neighbor[d] += disp[d];
          }

          const auto coordIndex = szIndex * StencilSize + si;

          if(deviceTopologyView.indexing().contains(neighbor))
          {
            const auto neighborIndex =
              deviceTopologyView.indexing().LogicalIndexToIndex(neighbor);

            for(axom::IndexType m = 0; m < matCount; m++)
            {
              // Ask the neighbor zone for how much of the current material it contains.
              const auto fragmentIndex = offset + m;
              double vf = 0.;
              matsetView.zoneContainsMaterial(
                neighborIndex,
                sortedMaterialIdsView[fragmentIndex],
                vf);

              // Store the vf into the stencil for the current material.
              const auto destIndex = fragmentIndex * StencilSize + si;
              zoneMatStencilView[destIndex] = vf;
            }

            // coord stencil
            xcStencilView[coordIndex] = xview[neighborIndex];
            ycStencilView[coordIndex] = yview[neighborIndex];
            zcStencilView[coordIndex] = zview.empty() ? 0. : zview[neighborIndex];
          }
          else
          {
            // All of the material contributions for the neighbor are 0.
            for(axom::IndexType m = 0; m < matCount; m++)
            {
              // Store the vf into the stencil for the current material.
              const auto fragmentIndex = offset + m;
              const auto destIndex = fragmentIndex * StencilSize + si;
              zoneMatStencilView[destIndex] = 0.;
            }

            // coord stencil
            xcStencilView[coordIndex] = 0.;
            ycStencilView[coordIndex] = 0.;
            zcStencilView[coordIndex] = 0.;
          }
        }
      });
    AXOM_ANNOTATE_END("stencil");

    AXOM_ANNOTATE_BEGIN("vectors");
    constexpr int numVectorComponents = 3;
    const auto vecSize = numFragments * numVectorComponents;
    axom::Array<double> fragmentVectors(vecSize, vecSize, allocatorID);
    auto fragmentVectorsView = fragmentVectors.view();

    axom::for_all<ExecSpace>(
      matZoneView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        // The selected zone index in the whole mesh.
        const auto zoneIndex = matZoneView[szIndex];
        const auto matCount = matCountView[szIndex];
        // Where to begin writing this zone's fragment data.
        const auto offset = matOffsetView[szIndex];

        // The following lines are adapted from mira1.c:262

        // Compute Jacobian here since we have coordinate stencil data.
        double jac[3][3];
        const auto coordIndex = szIndex * StencilSize;
        const double *xcStencil = xcStencilView.data() + coordIndex;
        const double *ycStencil = ycStencilView.data() + coordIndex;
        const double *zcStencil = zcStencilView.data() + coordIndex;
        elvira::computeJacobian(xcStencil, ycStencil, zcStencil, ndims, jac);

        // Produce normal for each material in this zone.
        int iskip = 0;
        axom::mir::elvira::elvira(zoneIndex,
               matCount,
               offset,
               zoneMatStencilView,
               fragmentVectorsView,
               iskip,
               ndims);

        // Transform the normals.
        for(axom::IndexType i = 0; i < matCount; i++)
        {
          double *normal =
            fragmentVectorsView.data() + ((offset + i) * numVectorComponents);
          elvira::transform(normal, jac);
        }
      });
    AXOM_ANNOTATE_END("vectors");

    // intermediate...

    // Now that we have all of the normals for each zone, clip the zones using those normals, making a fragment
    // shape that has the right VF. then save the shape, take the remaining piece and apply the next normal to it.

    // Add in clean zones - do we end up with a vector of polyhedral zones here that might be easier to pass into
    // the mapping stage? Or, do we make a BG topology for the polyhedra and make sure the mapper can eat that?

    // Make matset for combined mesh
  }

private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  MatsetView m_matsetView;
};

}  // end namespace mir
}  // end namespace axom

#endif

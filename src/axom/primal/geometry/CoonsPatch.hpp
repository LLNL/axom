// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file CoonsPatch.hpp
 *
 * \brief A CoonsPatch primitive
 */

#ifndef AXOM_PRIMAL_COONSPATCH_HPP_
#define AXOM_PRIMAL_COONSPATCH_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include "axom/primal/operators/squared_distance.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T>
class CoonsPatch;

/*! \brief Overloaded output operator for Coons Patches*/
template <typename T>
std::ostream& operator<<(std::ostream& os, const CoonsPatch<T>& cPatch);

/*!
 * \class CoonsPatch
 *
 * \brief A 3D Coons patch defined by an array of Bezier curves
 * \tparam T the coordinate type, e.g., double, float, etc.
 *
 * Performs transfinite interpolation between a sequence of Bezier curves.
 * Generalizes the standard Triangle and Quadrilateral Coons patches to N curves.
 * The patch is parameterized by Cartesian coordinates (u, v) defined on a
 * regular N-gon of radius 1.0. 
 * 
 * Contains a closed CurvePolygon object.
 */
template <typename T>
class CoonsPatch
{
public:
	using CurvedPolygonType = primal::CurvedPolygon<T, 3>;
};

}  // namespace primal
}  // namespace axom

#endif
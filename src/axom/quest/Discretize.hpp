// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_DISCRETIZE_HPP_
#define QUEST_DISCRETIZE_HPP_

// Axom includes
#include "axom/core/Macros.hpp"  // for axom macros

// Geometry
#include "axom/primal/geometry/Sphere.hpp"  
#include "axom/primal/geometry/Octahedron.hpp"
#include "axom/primal/geometry/Ray.hpp"

// C/C++ includes
#include <vector>  // for std::vector

namespace axom
{
namespace quest
{

/// \name Discretize primitive shapes to linear shapes
/// @{

using SphereType = primal::Sphere<double, 3>;
using OctType = primal::Octahedron<double, 3>;
using RayType = primal::Ray<double, 3>;
using TwoDPointType = primal::Point<double, 2>;

/*!
 * \brief Given a primitive sphere and a refinement level, return a list
 *   of Octahedra approximating the shape.
 * \param [in] s The sphere to approximate
 * \param [in] levels The number of refinements to perform
 * \param [out] out The collection of octahedra representing s
 *
 * This routine generates O(4^level) octahedra.  That's exponential growth.
 * Use appropriate caution.
 */
void discretize(const SphereType & s, int levels, std::vector<OctType> & out);

/*!
 * \brief Given a 2D polyline revolved around the positive X-axis, return a list
 *   of Octahedra approximating the shape.
 * \param [in] polyline The polyline to revolve around the X-axis
 * \param [in] levels The number of refinements to perform
 * \param [out] out The collection of octahedra representing the revolved
 *   polyline
 *
 * This routine generates n*O(2^level) octahedra, where n is the number of
 * segments in \a polyline (one less than the length).
 * That's exponential growth.  Use appropriate caution.
 */
void discretize(std::vector<TwoDPointType> & polyline,
                int levels,
                std::vector<OctType> & out);

/// @}

}  // end namespace quest
}  // end namespace axom

#endif /* QUEST_DISCRETIZE_HPP_ */

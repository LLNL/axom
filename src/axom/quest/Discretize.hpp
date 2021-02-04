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
 * \brief Return a discretized unit sphere.
 *
 * The octahedra returned by this routine were hand-computed.  This routine
 * is useful for testing or for use until the general discretize routine is
 * ready.
 */
void discretized_sphere(std::vector<OctType> & out);
   
/// @}

}  // end namespace quest
}  // end namespace axom

#endif /* QUEST_DISCRETIZE_HPP_ */

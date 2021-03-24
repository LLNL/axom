// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/Discretize.hpp"
#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/operators/squared_distance.hpp"

#include <cmath>

namespace axom
{
namespace quest
{
constexpr double PTINY = 1e-80;

using PointType = primal::Point<double, 3>;
using NAType = primal::NumericArray<double, 3>;

/* Project a Point onto a sphere.
 * Do we want to modify the passed-in point, or return a new one?
 */
PointType project_to_shape(const PointType &p, const SphereType &sphere)
{
  const double *ctr = sphere.getCenter();
  double dist2 = primal::squared_distance(ctr, p.data(), 3);
  double dist = sqrt(dist2);
  double drat = sphere.getRadius() * dist / (dist2 + PTINY);
  double dratc = drat - 1.0;
  return PointType::make_point(drat * p[0] - dratc * ctr[0],
                               drat * p[1] - dratc * ctr[1],
                               drat * p[2] - dratc * ctr[2]);
}

/* Return an octahedron whose six points lie on the given sphere.
 */
OctType from_sphere(const SphereType &sphere)
{
  NAType center(sphere.getCenter());
  NAType ihat({1., 0., 0.});
  NAType jhat({0., 1., 0.});
  NAType khat({0., 0., 1.});

  NAType dp = center + ihat;
  NAType dq = center + jhat;
  NAType dr = center + khat;
  NAType ds = center - ihat;
  NAType dt = center - jhat;
  NAType du = center - khat;

  PointType P = project_to_shape(PointType(dp), sphere);
  PointType Q = project_to_shape(PointType(dq), sphere);
  PointType R = project_to_shape(PointType(dr), sphere);
  PointType S = project_to_shape(PointType(ds), sphere);
  PointType T = project_to_shape(PointType(dt), sphere);
  PointType U = project_to_shape(PointType(du), sphere);

  return OctType(P, Q, R, S, T, U);
}

/* Given a sphere, a parent octahedron with vertices lying on the
 * sphere, and vertex indices s, t, u defining a face on that
 * octahedron, return a new oct sharing the face (s,t,u) and all other
 * faces looking "outward" toward the sphere.
 */
OctType new_inscribed_oct(const SphereType &sphere, OctType &o, int s, int t, int u)
{
  PointType P = PointType::midpoint(o[t], o[u]);
  PointType Q = PointType::midpoint(o[s], o[u]);
  PointType R = PointType::midpoint(o[s], o[t]);

  P = project_to_shape(P, sphere);
  Q = project_to_shape(Q, sphere);
  R = project_to_shape(R, sphere);

  return OctType(P, Q, R, o[s], o[t], o[u]);
}

/* Given a primitive shape and level of refinement, return the list of
 * octahedra that approximate the shape at that LOR.
 *
 * As noted in the doxygen, this function chops a shape into O(4^levels)
 * octahedra.  That is exponential growth.  Use appropriate caution.
 */
void discretize(const SphereType &sphere, int levels, std::vector<OctType> &out)
{
  // How many octahedra will we generate?  One oct for the central octahedron
  // (level 0), eight more for its faces (level 1), and four more to refine
  // each exposed face of the last generation.
  //
  // Total = 1 + 8*sum[i=0 to levels-1](4^i)
  int octcount = 1;
  for(int level = levels; level > 0; --level)
  {
    octcount *= 4;
    if(level == 1)
    {
      octcount *= 2;
    }
    octcount += 1;
  }
  out.reserve(octcount);

  // Establish an octahedron with all of its points lying on the sphere (level 0)
  out.emplace_back(from_sphere(sphere));

  // last_gen indexes to an octahedron of the last generation.
  int last_gen = 0;

  // max_last_gen indexes to the last oct of the last generation.
  int max_last_gen = 0;

  enum
  {
    P = 0,
    Q,
    R,
    S,
    T,
    U
  };

  // Refine: add an octahedron to each exposed face.
  for(int level = 0; level < levels; ++level)
  {
    max_last_gen = out.size();
    while(last_gen <= max_last_gen)
    {
      // Octahedra are defined by points P, Q, R, S, T, U (indexes 0--5).
      // Point oct[i] is opposite point oct[(i+3)%6].
      // Convention for new octahedra: P, Q, R are new points from midpoints,
      // S, T, U are inherited from the last gen.
      /* newoct[0] uses (P,Q,R)-old. */
      /* newoct[1] uses (T,P,R)-old. */
      /* newoct[2] uses (P,U,Q)-old. */
      /* newoct[3] uses (R,Q,S)-old. */
      out.push_back(new_inscribed_oct(sphere, out[last_gen], P, Q, R));
      out.push_back(new_inscribed_oct(sphere, out[last_gen], T, P, R));
      out.push_back(new_inscribed_oct(sphere, out[last_gen], P, U, Q));
      out.push_back(new_inscribed_oct(sphere, out[last_gen], R, Q, S));
      if(last_gen == 0)
      {
        out.push_back(new_inscribed_oct(sphere, out[last_gen], P, T, U));
        out.push_back(new_inscribed_oct(sphere, out[last_gen], Q, U, S));
        out.push_back(new_inscribed_oct(sphere, out[last_gen], T, R, S));
        out.push_back(new_inscribed_oct(sphere, out[last_gen], U, T, S));
      }

      last_gen += 1;
    }
  }
}

}  // end namespace quest
}  // end namespace axom

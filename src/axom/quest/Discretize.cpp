// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/Discretize.hpp"
#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/operators/squared_distance.hpp"

#include <cmath>

namespace
{
enum
{
  P = 0,
  Q,
  R,
  S,
  T,
  U
};
}

namespace axom
{
namespace quest
{
constexpr double PTINY = 1e-80;

using PointType = primal::Point<double, 3>;
using NAType = primal::NumericArray<double, 3>;

/* ------------------------------------------------------------ */
/* Project a Point onto a sphere.
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

PointType rescale_YZ(const PointType &p, double new_dst)
{
  double cur_dst = sqrt(p[1] * p[1] + p[2] * p[2]);
  PointType retval;
  retval[0] = p[0];
  retval[1] = p[1] * new_dst / (cur_dst + PTINY);
  retval[2] = p[2] * new_dst / (cur_dst + PTINY);
  return retval;
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

/* Return an octahedron whose six points lie on the truncated cone
 * described by rotating the line segment ab around the positive X-axis
 * in a right-handed way (thumb points out the axis, fingers spin segment
 * toward wrist).
 */
OctType from_segment(const TwoDPointType &a, const TwoDPointType &b)
{
  const double SQ_3_2 = sqrt(3.) / 2.;

  PointType p {a[0], 0., a[1]};
  PointType q {b[0], -SQ_3_2 * b[1], -0.5 * b[1]};
  PointType r {a[0], -SQ_3_2 * a[1], -0.5 * a[1]};
  PointType s {b[0], SQ_3_2 * b[1], -0.5 * b[1]};
  PointType t {a[0], SQ_3_2 * a[1], -0.5 * a[1]};
  PointType u {b[0], 0., b[1]};

  return OctType(p, q, r, s, t, u);
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

OctType new_inscribed_prism(OctType &old_oct,
                            int p_off,
                            int s_off,
                            int t_off,
                            int u_off,
                            TwoDPointType pa,
                            TwoDPointType pb)
{
  OctType retval;
  retval[P] = old_oct[p_off];
  PointType new_q = PointType::lerp(old_oct[u_off], old_oct[s_off], 0.5);
  retval[Q] = rescale_YZ(new_q, pb[1]);
  PointType new_r = PointType::lerp(old_oct[p_off], old_oct[t_off], 0.5);
  retval[R] = rescale_YZ(new_r, pa[1]);
  retval[S] = old_oct[s_off];
  retval[T] = old_oct[t_off];
  retval[U] = old_oct[u_off];
  return retval;
}

/* Given a sphere and level of refinement, place the list of octahedra that
 * approximate that sphere at that LOR into the output argument.  Return true for
 * valid input and no errors, false for bad input or an error.  If we return
 * false, no octahedra are placed in the output argument.
 *
 * As noted in the doxygen, this function chops a sphere into O(4^levels)
 * octahedra.  That is exponential growth.  Use appropriate caution.
 */
bool discretize(const SphereType &sphere, int levels, std::vector<OctType> &out)
{
  // Check input.  Negative radius: return false.
  if (sphere.getRadius() < 0) { return false; }
  // Zero radius: return true without generating octahedra.
  if (sphere.getRadius() < PTINY) { return true; }

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
      // Octahedra are defined by points P, Q, R, S, T, U (indexes 0--5; see
      // previous enum).  Point oct[i] is opposite point oct[(i+3)%6].
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

  return true;
}

/* ------------------------------------------------------------ */
/* Discretize a cylinder (or truncated cone) into a hierarchy of
 * triangular prisms (or truncated tetrahedra) stored as octahedra.
 * Each level of refinement places a new prism/tet on all exposed faces.
 *
 * Input:  2D points a and b, number of levels of refinement, index into
 * the output vector.  The routine assumes a.x <= b.x, a.y >= 0, b.y >= 0.
 *
 * Output: the output vector of prisms (placed in the output std::vector
 * starting at the index).
 *
 * Return value: the number of prisms generated: zero for degenerate segments,
 * or a value of order 2^(levels) calculated within the routine.
 *
 * Conceptually, end points a and b are revolved around the X-axis, describing
 * circles that are the truncated cone's end-caps.  The segment ab revolved
 * about the X-axis becomes the side-wall.
 *
 * The level-zero prism is constructed by inscribing a triangle in each of the
 * end-cap circles.  In one circle, put the points PRT; in the other circle, put
 * the points UQS.  The edges UP, QR, and ST lie in the side-wall described by
 * rotating the segment ab about the x-axis and are each co-planar with the
 * X-axis.  (This five-sided prism, with two end-caps, three side-walls, and six
 * vertices, is stored in an Octahedron data record.  The edges PQ, RS, and TU
 * split the quadrilateral prism side-walls into pairs of coplanar triangles.)
 *
 * Each subsequent level of refinement adds a prism to each exposed
 * quadrilateral side-wall.
 */
int discrSeg(const TwoDPointType &a,
             const TwoDPointType &b,
             int levels,
             std::vector<OctType> &out,
             int idx)
{
  // Deal with degenerate segments
  if (abs(a[0] - b[0]) < PTINY) { return 0; }
  if (a[1] < PTINY && b[1] < PTINY) { return 0; }

  // Establish a prism (in an octahedron record) with one triangular
  // end lying on the circle described by rotating point a around the
  // x-axis and the other lying on circle from rotating b.
  out[idx + 0] = from_segment(a, b);

  // curr_lvl indexes to the first prism in the level we're currently refining
  int curr_lvl = idx;
  int curr_lvl_count = 1;
  int next_lvl = curr_lvl + curr_lvl_count;
  int total_count = 0;

  TwoDPointType pa, pb;

  // Refine: add an octahedron to each exposed face.
  for(int level = 0; level < levels; ++level)
  {
    // Each level of refinement generates a prism for each exposed
    // face, so twice the prisms in the preceding level.  Refining the
    // initial prism is the only different step, since all three of its
    // side-faces are exposed.
    int lvl_factor = 2;
    if(level == 0)
    {
      lvl_factor = 3;
    }

    // The ends of the prisms switch each level.
    if(level & 1)
    {
      pa = a;
      pb = b;
    }
    else
    {
      pa = b;
      pb = a;
    }

    // This loop generates the prisms of the next level of refinement.
    // The specified vertices ensure that the new prisms always have
    // triangular end-caps QSU and RTP, and that side-face PTSU faces
    // the parent-level prism.
    // Of note, the child-level end-cap QSU is coplanar with parent-
    // level cap RTP, and vice versa.  Hence the preceding if-statement
    // with comment "the ends switch each level."
    for(int i = 0; i < curr_lvl_count; ++i)
    {
      out[next_lvl + i * lvl_factor + 0] =
        new_inscribed_prism(out[curr_lvl + i], Q, T, S, R, pa, pb);
      out[next_lvl + i * lvl_factor + 1] =
        new_inscribed_prism(out[curr_lvl + i], U, R, Q, P, pa, pb);
      if(curr_lvl == 0)
      {
        out[next_lvl + i * lvl_factor + 2] =
          new_inscribed_prism(out[curr_lvl + i], S, P, U, T, pa, pb);
      }
    }

    curr_lvl = next_lvl;
    curr_lvl_count *= lvl_factor;
    next_lvl = curr_lvl + curr_lvl_count;
    total_count += curr_lvl_count;
  }

  return total_count;
}

/* Given a surface of revolution and level of refinement, place the list of
 * octahedra that approximate that shape at that LOR in an output argument.
 * Return true for valid input and lack of errors, false otherwise.  If we
 * return false, put nothing in the output argument.
 *
 * As noted in the doxygen, this function chops a surface of revolution into
 * n*O(2^levels) octahedra, where n is the number of segments in the SoR (one
 * less than the polyline's length).  That is exponential growth.  Use
 * appropriate caution.
 */
bool discretize(std::vector<TwoDPointType> &polyline,
                int levels,
                std::vector<OctType> &out)
{
  // Check for invalid input.  If any segment is invalid, exit returning false.
  bool stillValid = true;
  int segmentcount = polyline.size() - 1;
  for (int seg = 0; seg < segmentcount && stillValid; ++seg)
  {
    TwoDPointType & a = polyline[seg];
    TwoDPointType & b = polyline[seg + 1];
    // invalid if a.x > b.x
    if (a[0] > b[0]) { stillValid = false; }
    if (a[1] < 0 || b[1] < 0) { stillValid = false; }
  }
  if (!stillValid) { return false; }

  // How many octahedra will we generate in each segment of the polyline?
  // This is done in discretizeSegment() (which see for further details).
  // - One oct for the central octahedron (level 0),
  // - three more for its side faces (level 1),
  // - and for level i > 1, two times the octahedron count in level i-1,
  //   to refine each exposed face.
  //
  // Total = 1 + 3*sum[i=0 to levels-1](2^i)
  int octcount = 1;
  for(int level = levels; level > 0; --level)
  {
    if(level == 1)
    {
      octcount *= 3;
    }
    else
    {
      octcount *= 2;
    }
    octcount += 1;
  }

  // That was the octahedron count for one segment.  Multiply by the number
  // of segments we will compute.
  int totaloctcount = octcount * segmentcount;
  out.resize(totaloctcount);
  int total_prism_count = 0;

  for(int seg = 0; seg < segmentcount; ++seg)
  {
    int segment_prism_count = discrSeg(polyline[seg], polyline[seg + 1], levels, out, total_prism_count);
    total_prism_count += segment_prism_count;
  }

  out.resize(total_prism_count);

  // TODO check for errors in each segment's computation
  return true;
}

}  // end namespace quest
}  // end namespace axom

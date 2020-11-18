// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#define _USE_MATH_DEFINES
#include <cmath>

#include "axom/quest/Discretize.hpp"
#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/operators/squared_distance.hpp"

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
PointType project_to_shape(const PointType & p, const SphereType &sphere)
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
OctType from_sphere(const SphereType & sphere)
{
   NAType center(sphere.getCenter());
   NAType ihat({1., 0., 0.});
   NAType jhat({0., 1., 0.});
   NAType khat({0., 0., 1.});

   NAType dp = center; dp += ihat;
   NAType dq = center; dq += jhat;
   NAType dr = center; dr += khat;
   NAType ds = center; ds -= ihat;
   NAType dt = center; dt -= jhat;
   NAType du = center; du -= khat;

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
OctType new_inscribed_oct(const SphereType & sphere,
                          OctType & o,
                          int s,
                          int t,
                          int u)
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
 * As noted in the doxygen, this function chops a shape into O(levels^4)
 * octahedra.  That is exponential growth.  Use appropriate caution.
 */
void discretize(const SphereType & sphere,
                int levels,
                std::vector<OctType> & out)
{
   // How many octahedra will we generate?  One oct for the central octahedron
   // (level 0), eight more for its faces (level 1), and four more to refine
   // each exposed face of the last generation.
   // 
   // Total = 1 + 8*sum[i=0 to levels-1](4^i)
   int octcount = 1;
   for (int level = levels; level > 0; --level)
   {
      octcount *= 4;
      if (level == 1) { octcount *= 2; }
      octcount += 1;
   }
   out.reserve(octcount);
   
   // Establish an octahedron with all of its points lying on the sphere (level 0)
   out.push_back(from_sphere(sphere));

   // last_gen indexes to an octahedron of the last generation.
   int last_gen = 0;

   // max_last_gen indexes to the last oct of the last generation.
   int max_last_gen = 0;

   // Refine: add an octahedron to each exposed face.
   for (int level = 0; level < levels; ++level)
   {
      max_last_gen = out.size();
      while (last_gen <= max_last_gen)
      {
         // Octahedra are defined by points P, Q, R, S, T, U (indexes 0--5).
         // Point oct[i] is opposite point oct[(i+3)%6].
         // Convention for new octahedra: P, Q, R are new points from midpoints,
         // S, T, U are inherited from the last gen.
         /* newoct[0] uses (P,Q,R)-old. */
         /* newoct[1] uses (T,P,R)-old. */
         /* newoct[2] uses (P,U,Q)-old. */
         /* newoct[3] uses (R,Q,S)-old. */
         out.push_back(new_inscribed_oct(sphere, out[last_gen], 0, 1, 2));
         out.push_back(new_inscribed_oct(sphere, out[last_gen], 4, 0, 2));
         out.push_back(new_inscribed_oct(sphere, out[last_gen], 0, 5, 1));
         out.push_back(new_inscribed_oct(sphere, out[last_gen], 2, 1, 3));
         if (last_gen == 0)
         {
            out.push_back(new_inscribed_oct(sphere, out[last_gen], 0, 4, 5));
            out.push_back(new_inscribed_oct(sphere, out[last_gen], 1, 5, 3));
            out.push_back(new_inscribed_oct(sphere, out[last_gen], 4, 2, 3));
            out.push_back(new_inscribed_oct(sphere, out[last_gen], 5, 4, 3));
         }
         
         last_gen += 1;
      }
   }
}

enum ReflectDimension
{
   X = 0,
   Y,
   Z
};

/* Reflect an octahedron in a given dimension.
 */
OctType reflect(ReflectDimension d, OctType o)
{
   OctType out(o);
   for (int i = 0; i < OctType::NUM_OCT_VERTS; ++i)
   {
      PointType & pt = out[i];
      pt[d] *= -1;
   }
   return out;
}

/* Return a handwritten list of the octahedra discretizing the unit sphere.
 */
void discretized_sphere(std::vector<OctType> & out)
{
   // We're going to return three generations in the out-vector:
   // one in the first generation, eight in the second (covering each
   // of the faces of the first generation), 32 in the third (covering
   // the four exposed faces in each of the second-gen octs).
   constexpr int FIRST_GEN_COUNT = 1;
   constexpr int SECOND_GEN_COUNT = 8;
   constexpr int THIRD_GEN_COUNT = 32;
   out.resize(FIRST_GEN_COUNT + SECOND_GEN_COUNT + THIRD_GEN_COUNT);
   
   // First generation: one octahedron, with vertices on the unit vectors.
   NAType ihat({1., 0., 0.});
   NAType jhat({0., 1., 0.});
   NAType khat({0., 0., 1.});

   OctType center(PointType(ihat), PointType(jhat), PointType(khat),
                  PointType(-1 * ihat), PointType(-1 * jhat), 
                  PointType(-1 * khat));
   out[0] = center;
   
   // Second generation: eight octs, one for each face of the unit oct.
   // Point ij is halfway between (1, 0, 0) and (0, 1, 0),
   // point jk is halfway between (0, 1, 0) and (0, 0, 1),
   // point ki is halfway between (0, 0, 1) and (1, 0, 0).
   PointType ij(NAType({M_SQRT1_2, M_SQRT1_2, 0.}));
   PointType jk(NAType({0., M_SQRT1_2, M_SQRT1_2}));
   PointType ki(NAType({M_SQRT1_2, 0., M_SQRT1_2}));

   OctType second_gen(PointType(ihat), PointType(jhat), PointType(khat),
                      jk, ki, ij);
   out[1] = second_gen;
   out[2] = reflect(X, second_gen);
   out[3] = reflect(Y, second_gen);
   out[4] = reflect(Z, second_gen);
   out[5] = reflect(Y, reflect(X, second_gen));
   out[6] = reflect(Z, reflect(Y, second_gen));
   out[7] = reflect(X, reflect(Z, second_gen));
   out[8] = reflect(X, reflect(Y, reflect(Z, second_gen)));
   
   // Third generation: 32 new octs, one for each exposed face of the previous
   // generation.
   double SQRT1_6 = 1./sqrt(6.);
   // There are three interior points, derived from ij, jk, and ki.
   // Point a is halfway between ij and ki, at (1/sqrt(6))(2, 1, 1).
   PointType a(SQRT1_6*NAType({2, 1, 1}));
   // Point b is halfway between ij and jk, at (1/sqrt(6))(1, 2, 1).
   PointType b(SQRT1_6*NAType({1, 2, 1}));
   // Point c is halfway between jk and ki, at (1/sqrt(6))(1, 1, 2).
   PointType c(SQRT1_6*NAType({1, 1, 2}));
   
   // There are six edge points, derived from the original corner points and
   // ij, jk, and ki.
   // Point d is halfway between ihat and ij, at
   // (1/sqrt(4 + 2 sqrt(2)))(1+sqrt(2), 1, 0)
   double FACTOR_3G = 1./sqrt(2.*M_SQRT2 + 4);
   PointType d(FACTOR_3G*NAType({1. + M_SQRT2, 1., 0}));
   // Point e splits jhat and ij, at (1/sqrt(2 sqrt(2) + 2))(1, 1+sqrt(2), 0)
   PointType e(FACTOR_3G*NAType({1., 1. + M_SQRT2, 0}));
   // Point f splits jhat and jk, at (1/sqrt(2 sqrt(2) + 2))(0, 1+sqrt(2), 1)
   PointType f(FACTOR_3G*NAType({0, 1. + M_SQRT2, 1.}));
   // Point g splits khat and jk, at (1/sqrt(2 sqrt(2) + 2))(0, 1, 1+sqrt(2))
   PointType g(FACTOR_3G*NAType({0, 1., 1. + M_SQRT2}));
   // Point m splits khat and ki, at (1/sqrt(2 sqrt(2) + 2))(1, 0, 1+sqrt(2))
   PointType m(FACTOR_3G*NAType({1., 0, 1. + M_SQRT2}));
   // Point n splits ihat and ki, at (1/sqrt(2 sqrt(2) + 2))(1+sqrt(2), 0, 1)
   PointType n(FACTOR_3G*NAType({1. + M_SQRT2, 0, 1.}));

   int offset = FIRST_GEN_COUNT + SECOND_GEN_COUNT;
   // Here's the first octant of third-generation octahedra (octant 0).
   int octant = 0;
   // First, the interior oct
   out[offset+0+octant*4] = OctType(ij, jk, ki, c, a, b);
   // The one next to ihat
   out[offset+1+octant*4] = OctType(PointType(ihat), ij, ki, a, n, d);
   // The one next to jhat
   out[offset+2+octant*4] = OctType(PointType(jhat), jk, ij, b, e, f);
   // The one next to khat
   out[offset+3+octant*4] = OctType(PointType(khat), ki, jk, c, g, m);

   // Now we get to transform these into all seven remaining octants.
   // Reflect in X
   octant = 1;
   ReflectDimension rd0 = X;
   out[offset+0+octant*4] = reflect(rd0, out[offset+0]);
   out[offset+1+octant*4] = reflect(rd0, out[offset+1]);
   out[offset+2+octant*4] = reflect(rd0, out[offset+2]);
   out[offset+3+octant*4] = reflect(rd0, out[offset+3]);
   // Reflect in Y
   octant = 2;
   rd0 = Y;
   out[offset+0+octant*4] = reflect(rd0, out[offset+0]);
   out[offset+1+octant*4] = reflect(rd0, out[offset+1]);
   out[offset+2+octant*4] = reflect(rd0, out[offset+2]);
   out[offset+3+octant*4] = reflect(rd0, out[offset+3]);
   // Reflect in Z
   octant = 3;
   rd0 = Z;
   out[offset+0+octant*4] = reflect(rd0, out[offset+0]);
   out[offset+1+octant*4] = reflect(rd0, out[offset+1]);
   out[offset+2+octant*4] = reflect(rd0, out[offset+2]);
   out[offset+3+octant*4] = reflect(rd0, out[offset+3]);
   // Reflect in X, then Y
   octant = 4;
   rd0 = X;
   ReflectDimension rd1 = Y;
   out[offset+0+octant*4] = reflect(rd1, reflect(rd0, out[offset+0]));
   out[offset+1+octant*4] = reflect(rd1, reflect(rd0, out[offset+1]));
   out[offset+2+octant*4] = reflect(rd1, reflect(rd0, out[offset+2]));
   out[offset+3+octant*4] = reflect(rd1, reflect(rd0, out[offset+3]));
   // Reflect in Y, then Z
   octant = 5;
   rd0 = Y;
   rd1 = Z;
   out[offset+0+octant*4] = reflect(rd1, reflect(rd0, out[offset+0]));
   out[offset+1+octant*4] = reflect(rd1, reflect(rd0, out[offset+1]));
   out[offset+2+octant*4] = reflect(rd1, reflect(rd0, out[offset+2]));
   out[offset+3+octant*4] = reflect(rd1, reflect(rd0, out[offset+3]));
   // Reflect in Z, then X
   octant = 6;
   rd0 = Z;
   rd1 = X;
   out[offset+0+octant*4] = reflect(rd1, reflect(rd0, out[offset+0]));
   out[offset+1+octant*4] = reflect(rd1, reflect(rd0, out[offset+1]));
   out[offset+2+octant*4] = reflect(rd1, reflect(rd0, out[offset+2]));
   out[offset+3+octant*4] = reflect(rd1, reflect(rd0, out[offset+3]));
   // And the grand finale: reflect in Z, then Y, then X
   octant = 7;
   rd0 = Z;
   rd1 = Y;
   ReflectDimension rd2 = X;
   out[offset+0+octant*4] = reflect(rd2, reflect(rd1, reflect(rd0, out[offset+0])));
   out[offset+1+octant*4] = reflect(rd2, reflect(rd1, reflect(rd0, out[offset+1])));
   out[offset+2+octant*4] = reflect(rd2, reflect(rd1, reflect(rd0, out[offset+2])));
   out[offset+3+octant*4] = reflect(rd2, reflect(rd1, reflect(rd0, out[offset+3])));
}

}  // end namespace quest
}  // end namespace axom


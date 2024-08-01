// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file printers.hpp
 *
 * \brief Consists of jacob's neat little debugging print statements
 * 
 * \note If you see this in a repo, it reflects a failure on jacob's part.
 *       Reprimand him for his error.
 */

#ifndef PRIMAL_PRINTERS_HPP_
#define PRIMAL_PRINTERS_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/primal.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <regex>

namespace axom
{
namespace primal
{

template <typename T>
void python_print(const Point<T, 3>& point, const char* plot_lines = nullptr)
{
  printf("plot_query( fig, ax, (%.17f, %.17f, %.17f), ",
         point[0],
         point[1],
         point[2]);
  if(plot_lines) printf("plot_lines=\"%s\", ", plot_lines);
  printf("size=111)\n\n");
}

template <typename T>
void python_print(const Vector<T, 3>& vec,
                  const Point<T, 3>& point,
                  const char* color = "blue")
{
  printf("plot_vector( fig, ax, (%.17f, %.17f, %.17f), ", vec[0], vec[1], vec[2]);
  printf("origin=(%.17f, %.17f, %.17f),", point[0], point[1], point[2]);
  printf("color=\"%s\")\n", color);
}

template <typename T>
void desmos_print(const Point<T, 3>& point)
{
  printf("Q = [%.17f, %.17f, %.17f]\n", point[0], point[1], point[2]);
}

template <typename T>
void desmos_print(const Point<T, 2>& point)
{
  printf("Q = (%.17f, %.17f)\n", point[0], point[1]);
}

template <typename T>
void python_print(const BezierCurve<T, 3>& curve,
                  bool plot_tangent = false,
                  bool plot_endpoints = false,
                  char color = 'k')
{
  const int ord = curve.getOrder();

  printf("my_BezierCurve( [");
  printf("(%.17f,%.17f,%.17f)", curve[0][0], curve[0][1], curve[0][2]);
  for(int i = 1; i <= ord; ++i)
    printf(",(%.17f,%.17f,%.17f)", curve[i][0], curve[i][1], curve[i][2]);
  if(curve.isRational())
  {
    printf(" ], [%.17f", curve.getWeight(0));
    if(curve.isRational())
      for(int i = 1; i <= ord; ++i) printf(",%.17f", curve.getWeight(i));
  }
  printf("] ).plot( fig, ax, None, '%c'", color);
  if(plot_tangent) printf(", plot_tangent=True");
  if(plot_endpoints) printf(", plot_endpoints=True");
  printf(" )\n\n");
}

template <typename T>
void desmos_print(const BezierCurve<T, 2>& curve, int num = 0)
{
  printf("C(P_{%d}, W_{%d})\n", num, num);
  printf("P_{%d} = [(%.17f,%.17f)", num, curve[0][0], curve[0][1]);
  for(int p = 1; p <= curve.getOrder(); ++p)
    printf(",(%.17f,%.17f)", curve[p][0], curve[p][1]);
  printf("]\n");

  if(curve.isRational())
  {
    printf("W_{%d} = [%.17f", num, curve.getWeight(0));
    for(int p = 1; p <= curve.getOrder(); ++p)
      printf(",%.17f", curve.getWeight(p));
    printf("]\n");
  }
  else
  {
    printf("W_{%d} = [1", num);
    for(int p = 1; p <= curve.getOrder(); ++p) printf(",1");
    printf("]\n");
  }
}

template <typename T>
void desmos_print(const BezierCurve<T, 3>& curve)
{
  const int ord = curve.getOrder();
  for(int i = 0; i <= ord; ++i)
  {
    printf("P_{%d} = [%.17f, %.17f, %.17f]\n",
           i + 1,
           curve[i][0],
           curve[i][1],
           curve[i][2]);
  }

  if(curve.isRational())
  {
    printf("W = [%.17f", curve.getWeight(0));
    for(int i = 1; i <= ord; ++i)
    {
      printf(", %.17f", curve.getWeight(i));
    }
    printf("]\n");
  }
}

template <typename T>
void python_print(const BezierPatch<T, 3>& patch,
                  const char* cmap = "cm.Reds",
                  char cpcolor = '\0',
                  bool plot_normal = false)
{
  const int ord_u = patch.getOrder_u();
  const int ord_v = patch.getOrder_v();

  printf("my_BezierPatch(%d, %d, [", ord_u, ord_v);
  printf("(%.17f,%.17f,%.17f)", patch(0, 0)[0], patch(0, 0)[1], patch(0, 0)[2]);
  for(int i = 0; i <= ord_u; ++i)
    for(int j = (i == 0 ? 1 : 0); j <= ord_v; ++j)
      printf(",(%.17f,%.17f,%.17f)",
             patch(i, j)[0],
             patch(i, j)[1],
             patch(i, j)[2]);
  if(patch.isRational())
  {
    printf(" ], [%.17f", patch.getWeight(0, 0));
    if(patch.isRational())
      for(int i = 0; i <= ord_u; ++i)
        for(int j = (i == 0 ? 1 : 0); j <= ord_v; ++j)
          printf(",%.17f", patch.getWeight(i, j));
  }
  printf("] ).plot( fig, ax, %s", cmap);
  if(plot_normal) printf(", plot_normal=True");
  if(cpcolor != '\0') printf(", cpcolor='%c'", cpcolor);
  printf(")\n\n");
}

template <typename T>
void desmos_print(const BezierPatch<T, 3>& patch)
{
  // clang-format off
  const int ord_u = patch.getOrder_u();
  const int ord_v = patch.getOrder_v();

  printf("M = %d\n", ord_u);
  printf("N = %d\n", ord_v);

  // Print u basis functions
  for(int u = 0; u <= ord_u; ++u)
    printf("B_{%dM}(u) = nCr\\left(M, %d\\right)(1 - u)^{(M - %d)}u^{%d}\n", u + 1, u, u, u);

  // Print v basis functions
  for(int v = 0; v <= ord_v; ++v)
    printf("B_{%dN}(v) = nCr\\left(N, %d\\right)(1 - v)^{(N - %d)}v^{%d}\n", v + 1, v, v, v);

  // ===================== Print weight ======================
  printf("w(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("W_{%d}[%d]B_{%dM}(u)B_{%dN}(v)", v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nw_{u}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("W_{%d}[%d]B_{%dM}'(u)B_{%dN}(v)", v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nw_{v}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("W_{%d}[%d]B_{%dM}(u)B_{%dN}'(v)", v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nw_{uu}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("W_{%d}[%d]B_{%dM}''(u)B_{%dN}(v)", v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
    
  printf("\nw_{vv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("W_{%d}[%d]B_{%dM}(u)B_{%dN}''(v)", v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nw_{uv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("W_{%d}[%d]B_{%dM}'(u)B_{%dN}'(v)", v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  // ===================== Print X ======================
  printf("\nX(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("X_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}(v)",  v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
  
  printf("\nX_{u}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("X_{%d}[%d]W_{%d}[%d]B_{%dM}'(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nX_{v}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("X_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}'(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nX_{uu}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("X_{%d}[%d]W_{%d}[%d]B_{%dM}''(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
    
  printf("\nX_{vv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("X_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}''(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nX_{uv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("X_{%d}[%d]W_{%d}[%d]B_{%dM}'(u)B_{%dN}'(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
 
  // ===================== Print Y ======================
  printf("\nY(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Y_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
  
  printf("\nY_{u}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Y_{%d}[%d]W_{%d}[%d]B_{%dM}'(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nY_{v}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Y_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}'(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nY_{uu}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Y_{%d}[%d]W_{%d}[%d]B_{%dM}''(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
    
  printf("\nY_{vv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Y_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}''(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nY_{uv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Y_{%d}[%d]W_{%d}[%d]B_{%dM}'(u)B_{%dN}'(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  // ===================== Print Z ======================
  printf("\nZ(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Z_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
  
  printf("\nZ_{u}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Z_{%d}[%d]W_{%d}[%d]B_{%dM}'(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nZ_{v}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Z_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}'(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nZ_{uu}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Z_{%d}[%d]W_{%d}[%d]B_{%dM}''(u)B_{%dN}(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }
    
  printf("\nZ_{vv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Z_{%d}[%d]W_{%d}[%d]B_{%dM}(u)B_{%dN}''(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nZ_{uv}(u, v) = ");
  for(int u = 0; u <= ord_u; ++u)
    for(int v = 0; v <= ord_v; ++v)
    {
      printf("Z_{%d}[%d]W_{%d}[%d]B_{%dM}'(u)B_{%dN}'(v)", v + 1, u + 1, v + 1, u + 1, u + 1, v + 1);
      if(u != ord_u || v != ord_v) printf(" + ");
    }

  printf("\nP(u, v) = [X(u, v), Y(u, v), Z(u, v)]");
  printf("\nP_{u}(u, v) = [X_{u}(u, v), Y_{u}(u, v), Z_{u}(u, v)]");
  printf("\nP_{v}(u, v) = [X_{v}(u, v), Y_{v}(u, v), Z_{v}(u, v)]");
  printf("\nP_{uu}(u, v) = [X_{uu}(u, v), Y_{uu}(u, v), Z_{uu}(u, v)]");
  printf("\nP_{vv}(u, v) = [X_{vv}(u, v), Y_{vv}(u, v), Z_{vv}(u, v)]");
  printf("\nP_{uv}(u, v) = [X_{uv}(u, v), Y_{uv}(u, v), Z_{uv}(u, v)]");

  printf("\nS(u, v) = \\frac{P(u, v)}{w(u, v)}");
  printf("\nS_u(u, v) = \\frac{P_{u}(u, v) - S(u, v)w_{u}(u, v)}{w(u, v)}");
  printf("\nS_v(u, v) = \\frac{P_{v}(u, v) - S(u, v)w_{v}(u, v)}{w(u, v)}");
  printf("\nS_{uu}(u, v) = \\frac{P_{uu}(u, v) - 2S_u(u, v)w_u(u, v) - S(u, v)w_{uu}(u, v)}{w(u, v)}");
  printf("\nS_{vv}(u, v) = \\frac{P_{vv}(u, v) - 2S_v(u, v)w_v(u, v) - S(u, v)w_{vv}(u, v)}{w(u, v)}");
  printf("\nS_{uv}(u, v) = \\frac{P_{uv}(u, v) - S_u(u, v)w_v(u, v) - S_v(u, v)w_u(u, v) - S(u, v)w_{uv}(u, v)}{w(u, v)}");

  printf("\ns = 0.0000000001");
  printf("\nu_0 = ");
  printf("\nv_0 = ");
  
  printf("\nS_{uu}(u_0, v_0)[1]");
  printf("\nS_{uu}(u_0, v_0)[2]");
  printf("\nS_{uu}(u_0, v_0)[3]");
  
  printf("\nS_{vv}(u_0, v_0)[1]");
  printf("\nS_{vv}(u_0, v_0)[2]");
  printf("\nS_{vv}(u_0, v_0)[3]");

  printf("\nS_{uv}(u_0, v_0)[1]");
  printf("\nS_{uv}(u_0, v_0)[2]");
  printf("\nS_{uv}(u_0, v_0)[3]");
  printf("\n\n");
  for(int v = 0; v <= ord_v; ++v)
  {
    printf("W_{%d} = [%f", v + 1, patch.getWeight(0, v));
    for(int u = 1; u <= ord_u; ++u)
    {
      printf(", %f", patch.getWeight(u, v));
    }
    printf("]\nX_{%d} = [%f", v + 1, patch(0, v)[0]);
    for(int u = 1; u <= ord_u; ++u)
    {
      printf(", %f", patch(u, v)[0]);
    }
    printf("]\nY_{%d} = [%f", v + 1, patch(0, v)[1]);
    for(int u = 1; u <= ord_u; ++u)
    {
      printf(", %f", patch(u, v)[1]);
    }
    printf("]\nZ_{%d} = [%f", v + 1, patch(0, v)[2]);
    for(int u = 1; u <= ord_u; ++u)
    {
      printf(", %f", patch(u, v)[2]);
    }
    printf("]\n");
  }
  // clang-format on
}

template <typename T>
void python_print(const OrientedBoundingBox<T, 3> oBox,
                  const char* color = "blue")
{
  auto verts = oBox.vertices();
  printf("plot_box( fig, ax, [");
  printf("(%.17f,%.17f,%.17f)", verts[0][0], verts[0][1], verts[0][2]);
  for(int i = 1; i < verts.size(); ++i)
    printf(",(%.17f,%.17f,%.17f)", verts[i][0], verts[i][1], verts[i][2]);
  printf("], color='%s' )\n\n", color);
}

template <class Lambda, typename T>
void python_print(const OrientedBoundingBox<T, 3> oBox,
                  Lambda rotate_point,
                  numerics::Matrix<T> rotator,
                  const char* color = "blue")
{
  auto verts = oBox.vertices();
  printf("plot_box( fig, ax, [");
  auto rot_point = rotate_point(rotator, verts[0]);
  printf("(%.17f,%.17f,%.17f)", rot_point[0], rot_point[1], rot_point[2]);
  for(int i = 1; i < verts.size(); ++i)
  {
    rot_point = rotate_point(rotator, verts[i]);
    printf(",(%.17f,%.17f,%.17f)", rot_point[0], rot_point[1], rot_point[2]);
  }
  printf("], color='%s' )\n\n", color);
}

template <typename T>
void python_print(const BoundingBox<T, 3>& bBox, const char* color = "blue")
{
  std::vector<Point<T, 3>> verts;
  BoundingBox<T, 3>::getPoints(bBox, verts);
  printf("plot_bounding_box( fig, ax, [");
  printf("(%.17f,%.17f,%.17f)", verts[0][0], verts[0][1], verts[0][2]);
  for(int i = 1; i < verts.size(); ++i)
    printf(",(%.17f,%.17f,%.17f)", verts[i][0], verts[i][1], verts[i][2]);
  printf("], color='%s' )\n", color);
}

template <typename T>
void python_print(const Polygon<T, 3>& poly, bool closed = true)
{
  printf("plot_polygon(fig, ax, [");
  printf("(%.17f,%.17f,%.17f)", poly[0][0], poly[0][1], poly[0][2]);
  for(int i = 1; i < poly.numVertices(); ++i)
    printf(",(%.17f,%.17f,%.17f)", poly[i][0], poly[i][1], poly[i][2]);
  printf("]");
  if(!closed) printf(", closed=False");
  printf(")\n");
}

template <typename T>
void desmos_print(const Polygon<T, 2>& poly)
{
  printf("polygon(");
  printf("(%.17f,%.17f)", poly[0][0], poly[0][1]);
  for(int i = 1; i < poly.numVertices(); ++i)
    printf(",(%.17f,%.17f)", poly[i][0], poly[i][1]);
  printf(")\n");
}

template <typename T>
void convert_from_svg(std::string filename, axom::Array<BezierCurve<T, 2>>& curves)
{
  // Get regex pattern for all path elements
  std::regex pattern("<\\s*path[^>]*\\bd\\s*=\\s*[\"']([^\"']+)\"[^>]*>");

  // Read in file
  std::ifstream svg_file(filename);

  if(!svg_file.is_open())
  {
    std::cerr << "Failed to open the SVG file." << std::endl;
    return;  // Exit with an error code
  }

  std::string input((std::istreambuf_iterator<char>(svg_file)),
                    std::istreambuf_iterator<char>());

  std::sregex_iterator iter(input.begin(), input.end(), pattern);
  std::sregex_iterator end;

  while(iter != end)
  {
    std::smatch match = *iter;
    std::string dAttribute = match[1].str();

    // Replace all ',' with ' '
    dAttribute =
      std::regex_replace(dAttribute, std::regex(","), std::string(" "));

    // Put a space before and after every letter
    dAttribute = std::regex_replace(dAttribute,
                                    std::regex("([a-zA-Z])"),
                                    std::string(" $1 "));

    char command = '\0';
    double init_x = 0.0, init_y = 0.0;
    double curr_x = 0.0, curr_y = 0.0;
    double x1, y1, x2, y2, x, y;
    double dx1, dy1, dx2, dy2, dx, dy;

    // Iterate over words in the string
    std::istringstream iss(dAttribute);
    std::string word;

    while(iss >> word)
    {
      // Check if the word is a command
      if(isalpha(word[0]))
      {
        command = word[0];
        if(command == 'Z' || command == 'z')
        {
          axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                            Point<T, 2> {init_x, init_y}};

          // Add a line segment
          curves.push_back(BezierCurve<T, 2>(nodes, 1));

          // Update the current position
          curr_x = init_x;
          curr_y = init_y;
        }
      }
      else
      {
        // Check if the command is a move command
        if(command == 'M' || command == 'm')
        {
          // If the command is relative, add the value to the current position
          if(command == 'm')
          {
            dx = std::stod(word);
            iss >> dy;
            curr_x += dx;
            curr_y += dy;
            command = 'l';
          }
          else
          {
            x = std::stod(word);
            iss >> y;
            curr_x = x;
            curr_y = y;
            command = 'L';
          }

          // Set the initial position
          init_x = curr_x;
          init_y = curr_y;
        }
        else if(command == 'L' || command == 'l')
        {
          // If the command is relative, add the value to the current position
          if(command == 'l')
          {
            dx = std::stod(word);
            iss >> dy;
            axom::Array<Point<T, 2>> nodes = {
              Point<T, 2> {curr_x, curr_y},
              Point<T, 2> {curr_x + dx, curr_y + dy}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 1));

            // Update the current position
            curr_x += dx;
            curr_y += dy;
          }
          else
          {
            x = std::stod(word);
            iss >> y;
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                              Point<T, 2> {x, y}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 1));

            // Update the current position
            curr_x = x;
            curr_y = y;
          }
        }
        else if(command == 'C' || command == 'c')
        {
          // If the command is relative, add the value to the current position
          if(command == 'c')
          {
            dx1 = std::stod(word);
            iss >> dy1 >> dx2 >> dy2 >> dx >> dy;
            axom::Array<Point<T, 2>> nodes = {
              Point<T, 2> {curr_x, curr_y},
              Point<T, 2> {curr_x + dx1, curr_y + dy1},
              Point<T, 2> {curr_x + dx2, curr_y + dy2},
              Point<T, 2> {curr_x + dx, curr_y + dy}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 3));

            // Update the current position
            curr_x += dx;
            curr_y += dy;
          }
          else
          {
            x1 = std::stod(word);
            iss >> y1 >> x2 >> y2 >> x >> y;
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                              Point<T, 2> {x1, y1},
                                              Point<T, 2> {x2, y2},
                                              Point<T, 2> {x, y}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 3));

            // Update the current position
            curr_x = x;
            curr_y = y;
          }
        }
        else if(command == 'H' || command == 'h')
        {
          // If the command is relative, add the value to the current position
          if(command == 'h')
          {
            dx = std::stod(word);
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                              Point<T, 2> {curr_x + dx, curr_y}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 1));

            // Update the current position
            curr_x += dx;
          }
          else
          {
            x = std::stod(word);
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                              Point<T, 2> {x, curr_y}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 1));

            // Update the current position
            curr_x = x;
          }
        }
        else if(command == 'V' || command == 'v')
        {
          // If the command is relative, add the value to the current position
          if(command == 'v')
          {
            dy = std::stod(word);
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                              Point<T, 2> {curr_x, curr_y + dy}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 1));

            // Update the current position
            curr_y += dy;
          }
          else
          {
            y = std::stod(word);
            axom::Array<Point<T, 2>> nodes = {Point<T, 2> {curr_x, curr_y},
                                              Point<T, 2> {curr_x, y}};

            // Add a line segment
            curves.push_back(BezierCurve<T, 2>(nodes, 1));

            // Update the current position
            curr_y = y;
          }
        }
        else
        {
          std::cout << word << " is an unrecognized command. Whoopsie!"
                    << std::endl;
        }
      }
    }
    ++iter;
  }

  svg_file.close();
}

template <typename T>
BoundingBox<T, 2> curves_bbox(axom::Array<BezierCurve<T, 2>>& curves,
                              double scale_factor = 1.0,
                              bool make_square = false)
{
  // Find maximum and minimum x and y values
  double min_x = curves[0][0][0], min_y = curves[0][0][1];
  double max_x = curves[0][0][0], max_y = curves[0][0][1];

  for(auto& curve : curves)
  {
    for(int p = 0; p <= curve.getOrder(); ++p)
    {
      min_x = axom::utilities::min(min_x, curve[p][0]);
      max_x = axom::utilities::max(max_x, curve[p][0]);

      min_y = axom::utilities::min(min_y, curve[p][1]);
      max_y = axom::utilities::max(max_y, curve[p][1]);
    }
  }

  // Put those in a bounding box, then expand it
  primal::Point<double, 2> bounds[2] = {primal::Point<double, 2> {min_x, min_y},
                                        primal::Point<double, 2> {max_x, max_y}};
  primal::BoundingBox<double, 2> bb(bounds, 2);
  bb.scale(scale_factor);

  if(make_square)
  {
    int max_dim = bb.getLongestDimension();
    double max_len = bb.getMax()[max_dim] - bb.getMin()[max_dim];

    primal::Point<double, 2> centroid = bb.getCentroid();

    primal::Point<double, 2> new_min {centroid[0] - max_len / 2.0,
                                      centroid[1] - max_len / 2.0};
    primal::Point<double, 2> new_max {centroid[0] + max_len / 2.0,
                                      centroid[1] + max_len / 2.0};

    bb.addPoint(new_min);
    bb.addPoint(new_max);
  }

  return bb;
}

template <typename T>
void simple_grid_test(axom::Array<BezierCurve<T, 2>>& curves,
                      const BoundingBox<T, 2>& bb,
                      int npts_x,
                      int npts_y,
                      std::ofstream& wn_out)
{
  for(int xi = 0; xi < npts_x; ++xi)
  {
    double x = bb.getMin()[0] + xi * (bb.getMax()[0] - bb.getMin()[0]) / npts_x;

    //std::cout << x << std::endl;
    printLoadingBar((x - bb.getMin()[0]) / (bb.getMax()[0] - bb.getMin()[0]) * 100,
                    100);
    for(int yi = 0; yi < npts_y; ++yi)
	  {
			double y = bb.getMin()[1] + yi * (bb.getMax()[1] - bb.getMin()[1]) / npts_y;
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double wn = 0.0;
      for(int i = 0; i < curves.size(); ++i)
      {
        wn += primal::winding_number(query, curves[i], nevals, 1e-16, 1e-16);
      }

      wn_out << x << "," << y << "," << wn << std::endl;
    }
  }
}

template <typename T>
void simple_grid_test_memoized(
  axom::Array<std::unordered_map<std::pair<int, int>,
                                       primal::detail::BezierCurveMemo<double>,
                                       primal::detail::PairHash>>& marray,
  const BoundingBox<T, 2>& bb,
  int npts_x,
  int npts_y,
  std::ofstream& wn_out)
{
  primal::Polygon<double, 2> temp_approxogon(20);

  for(double x = bb.getMin()[0]; x <= bb.getMax()[0];
      x += (bb.getMax()[0] - bb.getMin()[0]) / npts_x)
  {
    //std::cout << x << std::endl;
    printLoadingBar((x - bb.getMin()[0]) / (bb.getMax()[0] - bb.getMin()[0]) * 100,
                    100);
    for(double y = bb.getMin()[1]; y <= bb.getMax()[1];
        y += (bb.getMax()[1] - bb.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double wn = 0.0;

      wn += primal::winding_number_approxogon_memoized(query,
                                                       marray,
                                                       temp_approxogon,
                                                       1e-16);

      wn_out << x << "," << y << "," << wn << std::endl;
    }
  }
}

template <typename T>
void simple_timing_test(axom::Array<BezierCurve<T, 2>>& curves,
                        const BoundingBox<T, 2>& bb,
                        int npts_x,
                        int npts_y,
                        std::ofstream& wn_out)
{
  for(double x = bb.getMin()[0]; x <= bb.getMax()[0];
      x += (bb.getMax()[0] - bb.getMin()[0]) / npts_x)
  {
    //std::cout << x << std::endl;
    printLoadingBar((x - bb.getMin()[0]) / (bb.getMax()[0] - bb.getMin()[0]) * 100,
                    100);
    for(double y = bb.getMin()[1]; y <= bb.getMax()[1];
        y += (bb.getMax()[1] - bb.getMin()[1]) / npts_y)
    {
      primal::Point<double, 2> query({x, y});

      int nevals = 0;

      double wn = 0.0;
      for(int i = 0; i < curves.size(); ++i)
      {
        wn += primal::winding_number(query, curves[i], nevals, 1e-16, 1e-16);
      }

      wn_out << x << "," << y << "," << wn << std::endl;
    }
  }
}
inline void printLoadingBar(int progress, int total, int barWidth = 40)
{
  float percentage = static_cast<float>(progress) / total;
  int progressWidth = static_cast<int>(percentage * barWidth);

  std::cout << "[";
  for(int i = 0; i < progressWidth; ++i) std::cout << "=";
  for(int i = progressWidth; i < barWidth; ++i) std::cout << "-";
  std::cout << "] " << std::setw(3) << static_cast<int>(percentage * 100.0)
            << "%\r";
  std::cout.flush();
}

}  // namespace primal

}  // namespace axom

#endif

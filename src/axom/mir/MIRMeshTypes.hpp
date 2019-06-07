// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef __MIR_MESH_TYPES_H__
#define __MIR_MESH_TYPES_H__


#include "axom/core.hpp"  // for axom macros
// #include "axom/mir.hpp"  // for Mir classes & functions
#include "axom/slam.hpp"


namespace numerics = axom::numerics;
namespace slam = axom::slam;

namespace axom
{
namespace mir
{
  /**
   * \brief Simple 2D Point class for example
   */
  struct Point2
  {
    Point2(double x = 0., double y = 0.) : m_x(x), m_y(y) {}

    Point2(const Point2& other) : m_x(other.m_x), m_y(other.m_y) {}

    Point2& operator=(const Point2& other)
    { m_x = other.m_x; m_y = other.m_y; return *this; }

    Point2& operator+=(const Point2& other)
    { m_x += other.m_x; m_y += other.m_y; return *this; }

    Point2& operator/=(double val)
    { m_x /= val; m_y += val; return *this; }

    double& operator[] (int i) { return (i==0) ? m_x : m_y; }
    const double& operator[] (int i) const { return (i==0) ? m_x : m_y; }

    friend std::ostream& operator<<(std::ostream& os, const Point2& pt)
    { return os << "{x:" << pt.m_x << ", y:" << pt.m_y <<"}"; }

    double m_x, m_y;
  };

  // SET TYPE ALIASES
  using PosType = slam::DefaultPositionType;
  using ElemType = slam::DefaultElementType;

  using ArrayIndir = slam::policies::ArrayIndirection< PosType, ElemType >;

  using VertSet = slam::PositionSet< PosType, ElemType >;
  using ElemSet = slam::PositionSet< PosType, ElemType >;

  // RELATION TYPE ALIASES
  using VarCard = slam::policies::VariableCardinality< PosType, ArrayIndir >;

  // Note: This is the actual relation type, which takes in a bunch of policies and data entries to relate to each other.
  // Note: It is the relation of the elements to the vertices.

  using ElemToVertRelation = slam::StaticRelation< PosType, ElemType, VarCard, ArrayIndir, ElemSet, VertSet >;
  using VertToElemRelation = slam::StaticRelation< PosType, ElemType, VarCard, ArrayIndir, VertSet, ElemSet >;

  // MAP TYPE ALIASES
  using BaseSet = slam::Set< PosType, ElemType >;
  using ScalarMap = slam::Map< BaseSet, axom::float64 >;  // Note: Documentation has a typo in it.
  using PointMap = slam::Map< BaseSet, Point2 >;
  
}
}
#endif
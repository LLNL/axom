// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Part class holds the material data for a single material.

#ifndef __TINY_HYDRO_TYPES_H__
#define __TINY_HYDRO_TYPES_H__


#include "VectorXY.hpp"

#include "fmt/fmt.hpp"

#include "axom/slam.hpp"

namespace slam = axom::slam;

namespace tinyHydro {

  typedef slam::Map<int>              IndexField;
  typedef int                         PositionType;
  typedef PositionType                IndexType;

  typedef slam::Map<double>           ScalarField;
  typedef ScalarField                 NodalScalarField;
  typedef ScalarField                 ZonalScalarField;

  typedef slam::Map<VectorXY>         VectorField;
  typedef VectorField                 NodalVectorField;
  typedef VectorField                 ZonalVectorField;
  typedef VectorField                 FaceVectorField;

  typedef VectorField                 BoundaryEdgeVectorField;


  typedef slam::PositionSet           ZoneSet;
  typedef slam::PositionSet           NodeSet;
  typedef slam::PositionSet           FaceSet;
  typedef slam::PositionSet           CornerSet;

  typedef slam::VectorIndirectionSet  ZoneSubset;
  typedef slam::VectorIndirectionSet  NodeSubset;



  enum
  {
    NODES_PER_ZONE        = 4,
    FACES_PER_ZONE        = 4,
    BD_BOTTOM             = 0,         // lower boundary nodes
    BD_RIGHT              = 1,         // right boundary nodes
    BD_TOP                = 2,         // top boundary nodes
    BD_LEFT               = 3,         // left boundary nodes
    NUM_DOMAIN_BOUNDARIES = 4
  };

  typedef slam::policies::STLVectorIndirection<PositionType, PositionType>              STLIndirection;

  typedef slam::policies::CompileTimeStride<PositionType, NODES_PER_ZONE>               ZNStride;
  typedef slam::policies::ConstantCardinality<PositionType, ZNStride>                   ZNCard;
  typedef slam::StaticRelation<ZNCard, STLIndirection, ZoneSet, NodeSet>                ZoneToNodeRelation;
  typedef ZoneToNodeRelation::RelationSubset                                            ZNodeSet;


  typedef slam::policies::CompileTimeStride<PositionType, FACES_PER_ZONE>               ZFStride;
  typedef slam::policies::ConstantCardinality<PositionType, ZFStride>                   ZFCard;
  typedef slam::StaticRelation<ZFCard, STLIndirection, ZoneSet, FaceSet>                ZoneToFaceRelation;
  typedef ZoneToFaceRelation::RelationSubset                                            ZFaceSet;


  typedef slam::policies::CompileTimeSize<ZoneSet::PositionType, NUM_DOMAIN_BOUNDARIES> NUM_BD_SZ;
  typedef slam::OrderedSet< PositionType, NUM_BD_SZ>                                    BoundaryEdgeSet;


  typedef slam::Map<IndexType>                                                          IndexMap;

  typedef slam::FieldRegistry<ZoneSet::IndexType>                                       IndexRegistry;
  typedef IndexRegistry::BufferType                                                     IndexBuffer;

  struct DataRegistry
  {
    static IndexRegistry setRegistry;
  };

} // end namespace tinyHydro

#endif

/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


// Part class holds the material data for a single material.

#ifndef __TINY_HYDRO_TYPES_H__
#define __TINY_HYDRO_TYPES_H__


#include "VectorXY.hpp"
#include "slam/OrderedSet.hpp"
#include "slam/RangeSet.hpp"
#include "slam/IndirectionSet.hpp"

#include "slam/StaticConstantRelation.hpp"

#include "slam/Map.hpp"

#include "slam/FieldRegistry.hpp"


namespace tinyHydro {

  typedef axom::slam::Map<int>              IndexField;


  typedef axom::slam::Map<double>           ScalarField;
  typedef ScalarField                             NodalScalarField;
  typedef ScalarField                             ZonalScalarField;

  typedef axom::slam::Map<VectorXY>         VectorField;
  typedef VectorField                             NodalVectorField;
  typedef VectorField                             ZonalVectorField;
  typedef VectorField                             FaceVectorField;

  typedef VectorField                             BoundaryEdgeVectorField;


  typedef axom::slam::PositionSet           ZoneSet;
  typedef axom::slam::PositionSet           NodeSet;
  typedef axom::slam::PositionSet           FaceSet;
  typedef axom::slam::PositionSet           CornerSet;

  typedef axom::slam::VectorIndirectionSet  ZoneSubset;
  typedef axom::slam::VectorIndirectionSet  NodeSubset;



  enum
  {
    NODES_PER_ZONE        = 4
    , FACES_PER_ZONE        = 4
    , BD_BOTTOM             = 0         // lower boundary nodes
    , BD_RIGHT              = 1         // right boundary nodes
    , BD_TOP                = 2         // top boundary nodes
    , BD_LEFT               = 3         // left boundary nodes
    , NUM_DOMAIN_BOUNDARIES = 4
  };
  typedef axom::slam::policies::CompileTimeStrideHolder<ZoneSet::PositionType, NODES_PER_ZONE>      ZNStride;
  typedef axom::slam::policies::CompileTimeStrideHolder<ZoneSet::PositionType, FACES_PER_ZONE>      ZFStride;

  typedef axom::slam::StaticConstantRelation<ZNStride, ZoneSet, NodeSet>                            ZoneToNodeRelation;
  typedef ZoneToNodeRelation::RelationSet                                                                 ZNodeSet;

  typedef axom::slam::StaticConstantRelation<ZFStride, ZoneSet, FaceSet>                            ZoneToFaceRelation;
  typedef ZoneToFaceRelation::RelationSet                                                                 ZFaceSet;


  typedef axom::slam::policies::CompileTimeSizeHolder<ZoneSet::PositionType, NUM_DOMAIN_BOUNDARIES> NUM_BD_SZ;
  typedef axom::slam::OrderedSet< NUM_BD_SZ>                                                        BoundaryEdgeSet;


  typedef ZoneSet::PositionType                                                                           IndexType;
  typedef axom::slam::Map<IndexType>                                                                IndexMap;

  typedef axom::slam::FieldRegistry<ZoneSet::IndexType>                                             SubsetRegistry;

  struct DataRegistry
  {
    static SubsetRegistry setRegistry;
  };

} // end namespace tinyHydro

#endif

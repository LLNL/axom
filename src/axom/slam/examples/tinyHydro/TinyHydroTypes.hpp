/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// Part class holds the material data for a single material.

#ifndef __TINY_HYDRO_TYPES_H__
#define __TINY_HYDRO_TYPES_H__


#include "VectorXY.hpp"

#include "fmt/fmt.hpp"

#include "axom/slam.hpp"


namespace tinyHydro {

  typedef axom::slam::Map<int>              IndexField;
  typedef int                               PositionType;
  typedef PositionType                      IndexType;

  typedef axom::slam::Map<double>           ScalarField;
  typedef ScalarField                       NodalScalarField;
  typedef ScalarField                       ZonalScalarField;

  typedef axom::slam::Map<VectorXY>         VectorField;
  typedef VectorField                       NodalVectorField;
  typedef VectorField                       ZonalVectorField;
  typedef VectorField                       FaceVectorField;

  typedef VectorField                       BoundaryEdgeVectorField;


  typedef axom::slam::PositionSet           ZoneSet;
  typedef axom::slam::PositionSet           NodeSet;
  typedef axom::slam::PositionSet           FaceSet;
  typedef axom::slam::PositionSet           CornerSet;

  typedef axom::slam::VectorIndirectionSet  ZoneSubset;
  typedef axom::slam::VectorIndirectionSet  NodeSubset;



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

  typedef axom::slam::policies::STLVectorIndirection<PositionType, PositionType>              STLIndirection;

  typedef axom::slam::policies::CompileTimeStride<PositionType, NODES_PER_ZONE>               ZNStride;
  typedef axom::slam::policies::ConstantCardinality<PositionType, ZNStride>                   ZNCard;
  typedef axom::slam::StaticRelation<ZNCard, STLIndirection, ZoneSet, NodeSet>                ZoneToNodeRelation;
  typedef ZoneToNodeRelation::RelationSubset                                                  ZNodeSet;


  typedef axom::slam::policies::CompileTimeStride<PositionType, FACES_PER_ZONE>               ZFStride;
  typedef axom::slam::policies::ConstantCardinality<PositionType, ZFStride>                   ZFCard;
  typedef axom::slam::StaticRelation<ZFCard, STLIndirection, ZoneSet, FaceSet>                ZoneToFaceRelation;
  typedef ZoneToFaceRelation::RelationSubset                                                  ZFaceSet;


  typedef axom::slam::policies::CompileTimeSize<ZoneSet::PositionType, NUM_DOMAIN_BOUNDARIES> NUM_BD_SZ;
  typedef axom::slam::OrderedSet< NUM_BD_SZ>                                                  BoundaryEdgeSet;


  typedef axom::slam::Map<IndexType>                                                          IndexMap;

  typedef axom::slam::FieldRegistry<ZoneSet::IndexType>                                       IndexRegistry;
  typedef IndexRegistry::BufferType                                                           IndexBuffer;

  struct DataRegistry
  {
    static IndexRegistry setRegistry;
  };

} // end namespace tinyHydro

#endif

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

  typedef asctoolkit::slam::Map<int>              IndexField;


  typedef asctoolkit::slam::Map<double>           ScalarField;
  typedef ScalarField                             NodalScalarField;
  typedef ScalarField                             ZonalScalarField;

  typedef asctoolkit::slam::Map<VectorXY>         VectorField;
  typedef VectorField                             NodalVectorField;
  typedef VectorField                             ZonalVectorField;
  typedef VectorField                             FaceVectorField;

  typedef VectorField                             BoundaryEdgeVectorField;


  typedef asctoolkit::slam::PositionSet           ZoneSet;
  typedef asctoolkit::slam::PositionSet           NodeSet;
  typedef asctoolkit::slam::PositionSet           FaceSet;
  typedef asctoolkit::slam::PositionSet           CornerSet;

  typedef asctoolkit::slam::VectorIndirectionSet  ZoneSubset;
  typedef asctoolkit::slam::VectorIndirectionSet  NodeSubset;



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
  typedef asctoolkit::slam::policies::CompileTimeStrideHolder<ZoneSet::PositionType, NODES_PER_ZONE>      ZNStride;
  typedef asctoolkit::slam::policies::CompileTimeStrideHolder<ZoneSet::PositionType, FACES_PER_ZONE>      ZFStride;

  typedef asctoolkit::slam::StaticConstantRelation<ZNStride, ZoneSet, NodeSet>                            ZoneToNodeRelation;
  typedef ZoneToNodeRelation::RelationSet                                                                 ZNodeSet;

  typedef asctoolkit::slam::StaticConstantRelation<ZFStride, ZoneSet, FaceSet>                            ZoneToFaceRelation;
  typedef ZoneToFaceRelation::RelationSet                                                                 ZFaceSet;


  typedef asctoolkit::slam::policies::CompileTimeSizeHolder<ZoneSet::PositionType, NUM_DOMAIN_BOUNDARIES> NUM_BD_SZ;
  typedef asctoolkit::slam::OrderedSet< NUM_BD_SZ>                                                        BoundaryEdgeSet;


  typedef ZoneSet::PositionType                                                                           IndexType;
  typedef asctoolkit::slam::Map<IndexType>                                                                IndexMap;

  typedef asctoolkit::slam::FieldRegistry<ZoneSet::IndexType>                                             SubsetRegistry;

  struct DataRegistry
  {
    static SubsetRegistry setRegistry;
  };

} // end namespace tinyHydro

#endif

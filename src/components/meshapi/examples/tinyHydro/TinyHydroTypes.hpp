// Part class holds the material data for a single material.

#ifndef __TINY_HYDRO_TYPES_H__
#define __TINY_HYDRO_TYPES_H__


#include "VectorXY.hpp"
#include "meshapi/OrderedSet.hpp"
#include "meshapi/RangeSet.hpp"
#include "meshapi/IndirectionSet.hpp"

#include "meshapi/StaticConstantRelation.hpp"

#include "meshapi/Map.hpp"

#include "meshapi/FieldRegistry.hpp"


namespace tinyHydro {

    typedef asctoolkit::meshapi::Map<int> IndexField;


    typedef asctoolkit::meshapi::Map<double> ScalarField;
    typedef ScalarField NodalScalarField;
    typedef ScalarField ZonalScalarField;

    typedef asctoolkit::meshapi::Map<VectorXY> VectorField;
    typedef VectorField NodalVectorField;
    typedef VectorField ZonalVectorField;
    typedef VectorField FaceVectorField;

    typedef VectorField BoundaryEdgeVectorField;


    typedef asctoolkit::meshapi::PositionSet ZoneSet;
    typedef asctoolkit::meshapi::PositionSet NodeSet;
    typedef asctoolkit::meshapi::PositionSet FaceSet;
    typedef asctoolkit::meshapi::PositionSet CornerSet;

    typedef asctoolkit::meshapi::VectorIndirectionSet ZoneSubset;
    typedef asctoolkit::meshapi::VectorIndirectionSet NodeSubset;



    enum { NODES_PER_ZONE       = 4
        , FACES_PER_ZONE        = 4
        , BD_BOTTOM             = 0     // lower boundary nodes
        , BD_RIGHT              = 1     // right boundary nodes
        , BD_TOP                = 2     // top boundary nodes
        , BD_LEFT               = 3     // left boundary nodes
        , NUM_DOMAIN_BOUNDARIES = 4
    };
    typedef asctoolkit::meshapi::policies::CompileTimeStrideHolder<ZoneSet::PositionType, NODES_PER_ZONE> ZNStride;
    typedef asctoolkit::meshapi::policies::CompileTimeStrideHolder<ZoneSet::PositionType, FACES_PER_ZONE> ZFStride;

    typedef asctoolkit::meshapi::StaticConstantRelation<ZNStride, ZoneSet, NodeSet> ZoneToNodeRelation;
    typedef ZoneToNodeRelation::RelationSet ZNodeSet;

    typedef asctoolkit::meshapi::StaticConstantRelation<ZFStride, ZoneSet, FaceSet> ZoneToFaceRelation;
    typedef ZoneToFaceRelation::RelationSet ZFaceSet;


    typedef asctoolkit::meshapi::policies::CompileTimeSizeHolder<ZoneSet::PositionType, NUM_DOMAIN_BOUNDARIES> NUM_BD_SZ;
    typedef asctoolkit::meshapi::OrderedSet< NUM_BD_SZ> BoundaryEdgeSet;


    typedef ZoneSet::PositionType IndexType;
    typedef asctoolkit::meshapi::Map<IndexType> IndexMap;

    typedef asctoolkit::meshapi::FieldRegistry<ZoneSet::IndexType> SubsetRegistry;

    struct DataRegistry
    {
        static SubsetRegistry setRegistry;
    };

} // end namespace tinyHydro

#endif

// arbitrary mesh of polygons in XY geom
// Thu Mar 26 09:38:50 PDT 2015

#ifndef __PolygonMeshXY_hh__
#define __PolygonMeshXY_hh__

#include "VectorXY.hpp"
#include <string>
#include <stdio.h>
#include <assert.h>

#include "meshapi/RangeSet.hpp"
#include "meshapi/StaticConstantRelation.hpp"
#include "meshapi/Map.hpp"

class PolygonMeshXY
{
 public:
    typedef asctoolkit::meshapi::PositionSet ZoneSet;
    typedef asctoolkit::meshapi::PositionSet NodeSet;
    typedef asctoolkit::meshapi::PositionSet FaceSet;
    typedef asctoolkit::meshapi::PositionSet CornerSet;

    enum { NODES_PER_ZONE = 4, FACES_PER_ZONE = 4};
    typedef asctoolkit::meshapi::policies::CompileTimeStrideHolder<ZoneSet::PositionType, NODES_PER_ZONE> ZNStride;
    typedef asctoolkit::meshapi::policies::CompileTimeStrideHolder<ZoneSet::PositionType, FACES_PER_ZONE> ZFStride;

    typedef asctoolkit::meshapi::StaticConstantRelation<ZNStride, ZoneSet, NodeSet> ZoneToNodeRelation;
    typedef ZoneToNodeRelation::RelationSet ZNodeSet;

    typedef asctoolkit::meshapi::StaticConstantRelation<ZFStride, ZoneSet, FaceSet> ZoneToFaceRelation;
    typedef ZoneToFaceRelation::RelationSet ZFaceSet;

    typedef ZoneSet::PositionType IndexType;
    typedef asctoolkit::meshapi::Map<IndexType> IndexMap;

    typedef asctoolkit::meshapi::Map<double>   ScalarField;
    typedef asctoolkit::meshapi::Map<VectorXY> VectorField;

    typedef ScalarField NodalScalarField;
    typedef ScalarField ZonalScalarField;

    typedef VectorField NodalVectorField;
    typedef VectorField ZonalVectorField;
    typedef VectorField FaceVectorField;

 public:
   // c'tor for doing quads
   PolygonMeshXY(int kmax, int lmax,
                double xmin = 0.0, double xmax = 1.0,
                double ymin = 0.0, double ymax = 1.0);

   // c'tor for reading in a mesh file
   enum MeshFileType {MT0, SILO, PDB, TRIANGLE};
   PolygonMeshXY(std::string /*filename*/, MeshFileType /*meshfiletype = MT0*/)
   {
      printf("not implemented yet\n");
      assert(false);
   }

   // destructor
   ~PolygonMeshXY();
   
   ZoneSet zones;                   // MeshAPI -- Set of zones -- PositionSet (wrapper around an int)
   NodeSet nodes;                   // MeshAPI -- Set of nodes -- PositionSet (wrapper around an int)
   FaceSet faces;                   // MeshAPI -- Set of faces -- missing (wrapper around an int)
   CornerSet corners;               // MeshAPI -- Set of corners -- missing (wrapper around an int)

   int numNodes() const { return nodes.size(); }
   int numZones() const { return zones.size(); }

   // connectivity arrays
   ZoneToFaceRelation zoneToFaces;
//   int * z2firstFace;           // MeshAPI -- StaticConstantRelation ( compile time cardinality of 4 )
//   int * zFaces;

   ZoneToNodeRelation zoneToNodes;
//   int * z2firstNode;           // MeshAPI -- StaticConstantRelation ( compile time cardinality of 4 )
//   int * zNodes;
//   int * zNumNodes;

   // geometry information
   NodalVectorField nodePos;          // MeshAPI -- Map: Node -> VectorXY -- node positions
   ZonalVectorField zonePos;          // MeshAPI -- Map: Zone -> VectorXY -- zone positions (barycenters)
   FaceVectorField faceArea;         // MeshAPI -- Map: Face -> VectorXY -- face areas -- outward normals
   ZonalScalarField zoneVolume;         // MeshAPI -- Map: Zone -> area (scalar)

   // functions that modify the mesh data
   void moveNodesToPosition(const VectorXY * newPos);
   void computeNewGeometry(void);
   
   // accessors, mostly used from Python
   VectorXY getPos(int i) const {return nodePos[i];}
   VectorXY getZonePos(int i) const {return zonePos[i];}
   double zoneVol(int i) const {return zoneVolume[i];}
   void setPos(int i, const VectorXY & v) {nodePos[i] = v;}

   int zoneNumNodes(int i) const {return zoneToNodes.size(i);}
   int zoneNode(int iz, int in) const {return zoneToNodes[iz][in];}
   VectorXY zoneNodePos(int iz, int in) const {return nodePos[ zoneNode(iz,in) ];}

   VectorXY meshAverageKLZMemOrderA();
   double timeElapsed;


   void dumpMesh();
} ;





#endif      // __PolygonMeshXY_hh__

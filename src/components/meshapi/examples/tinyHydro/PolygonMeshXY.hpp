// arbitrary mesh of polygons in XY geom
// Thu Mar 26 09:38:50 PDT 2015

#ifndef __PolygonMeshXY_hh__
#define __PolygonMeshXY_hh__

#include "VectorXY.hpp"
#include "TinyHydroTypes.hpp"

#include <string>
#include <stdio.h>

namespace tinyHydro {

class PolygonMeshXY
{

 public:
   // c'tor for doing quads
   PolygonMeshXY(int kmax, int lmax,
                double xmin = 0.0, double xmax = 1.0,
                double ymin = 0.0, double ymax = 1.0);

   // c'tor for reading in a mesh file
   enum MeshFileType {MT0, SILO, PDB, TRIANGLE};
   PolygonMeshXY(std::string /*filename*/, MeshFileType /*meshfiletype = MT0*/)
   {
      SLIC_ERROR("not implemented yet\n");
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
   ZoneToNodeRelation zoneToNodes;

   // geometry information
   NodalVectorField nodePos;          // MeshAPI -- Map: Node -> VectorXY -- node positions
   ZonalVectorField zonePos;          // MeshAPI -- Map: Zone -> VectorXY -- zone positions (barycenters)
   FaceVectorField faceArea;         // MeshAPI -- Map: Face -> VectorXY -- face areas -- outward normals
   ZonalScalarField zoneVolume;         // MeshAPI -- Map: Zone -> area (scalar)

   // functions that modify the mesh data
   void moveNodesToPosition(const NodalVectorField& newPos);
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

} // end namespace tinyHydro

#endif      // __PolygonMeshXY_hh__

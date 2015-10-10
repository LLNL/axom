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


   // functions that access the mesh topology
   int numNodes()                           const { return nodes.size(); }
   int numZones()                           const { return zones.size(); }
   int zoneNumNodes(int i)                  const { return zoneToNodes.size(i);}
   int zoneNode(int iz, int in)             const { return zoneToNodes[iz][in];}
   VectorXY zoneNodePos(int iz, int in)     const { return nodePos[ zoneNode(iz,in) ];}

   // functions that access and modify the mesh geometry
   void moveNodesToPosition(const NodalVectorField& newPos);
   void computeNewGeometry(void);
   
   // accessors, mostly used from Python
   VectorXY getPos(int i)                   const { return nodePos[i];}
   VectorXY getZonePos(int i)               const { return zonePos[i];}
   double zoneVol(int i)                    const { return zoneVolume[i];}
   void setPos(int i, const VectorXY & v)         { nodePos[i] = v;}

   VectorXY meshAverageKLZMemOrderA();

   void dumpMesh();

 public:
   ZoneSet zones;                   // Set of zones in the mesh -- PositionSet (wrapper around an int)
   NodeSet nodes;                   // Set of nodes in the mesh -- PositionSet (wrapper around an int)
   FaceSet faces;                   // Set of faces in the mesh -- missing (wrapper around an int)
   CornerSet corners;               // Set of corners in the mesh -- missing (wrapper around an int)

   // Mesh topology
   ZoneToFaceRelation zoneToFaces;
   ZoneToNodeRelation zoneToNodes;

   // Mesh geometry
   NodalVectorField nodePos;          // node positions
   ZonalVectorField zonePos;          // zone positions (barycenters)
   FaceVectorField faceArea;          // face areas -- scaled outward-facing normals
   ZonalScalarField zoneVolume;       // zone volumes (scalar)

   double timeElapsed;
} ;

} // end namespace tinyHydro

#endif      // __PolygonMeshXY_hh__

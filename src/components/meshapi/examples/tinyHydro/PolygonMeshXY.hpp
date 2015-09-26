// arbitrary mesh of polygons in XY geom
// Thu Mar 26 09:38:50 PDT 2015

#ifndef __PolygonMeshXY_hh__
#define __PolygonMeshXY_hh__

#include "VectorXY.hpp"
#include <string>
#include <stdio.h>
#include <assert.h>

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
      printf("not implemented yet\n");
      assert(false);
   }

   // destructor
   ~PolygonMeshXY();
   
   int nzones;
   int nnodes;

   // connectivity arrays
   int * z2firstFace;
   int * z2firstNode;
   int * zNodes;
   int * zFaces;
   int * zNumNodes;

   // geometry information
   VectorXY * nodePos;
   VectorXY * zonePos;
   VectorXY * faceArea;
   double * zoneVolume;

   // functions that modify the mesh data
   void moveNodesToPosition(const VectorXY * newPos);
   void computeNewGeometry(void);
   
   // accessors, mostly used from Python
   VectorXY getPos(int i) const {return nodePos[i];}
   VectorXY getZonePos(int i) const {return zonePos[i];}
   double zoneVol(int i) const {return zoneVolume[i];}
   void setPos(int i, const VectorXY & v) {nodePos[i] = v;}
   int zoneNumNodes(int i) const {return zNumNodes[i];}
   int zoneNode(int iz, int in) const {return zNodes[z2firstNode[iz]+in];}
   VectorXY zoneNodePos(int iz, int in) const {return nodePos[zNodes[z2firstNode[iz]+in]];}

   VectorXY meshAverageKLZMemOrderA();
   double timeElapsed;


   void dumpMesh();
} ;





#endif      // __PolygonMeshXY_hh__

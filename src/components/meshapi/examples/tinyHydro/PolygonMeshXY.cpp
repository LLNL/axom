#include "PolygonMeshXY.hpp"
#include "slic/slic.hpp"

//----------------------------------------------
PolygonMeshXY::PolygonMeshXY(int kmax, int lmax,
                           double xmin, double xmax,
                           double ymin, double ymax)
{
   // This c'tor constructs a mesh that looks like a logically regular
   // mesh of quads, but in fact it is arbitrary polygons.  It's
   // intended for making quick test meshes and to do direct speed
   // comparisons between the polygon and regular quad mesh.
   
   SLIC_ASSERT(kmax > 0 && lmax > 0);
   
   nodes = NodeSet(kmax*lmax);
   zones = ZoneSet((kmax-1)*(lmax-1));
   faces = FaceSet( zones.size() * 4);
   corners = CornerSet( zones.size() * 4);

   // make space for connectivity data
   zoneToNodes = ZoneToNodeRelation(&zones, &nodes);
   zoneToFaces = ZoneToFaceRelation(&zones, &faces);

   IndexMap faceIndices(&faces);
   IndexMap nodeIndices(&corners);



//   zNumNodes = new int[nzones];
//   z2firstNode = new int[nzones];
//   z2firstFace = new int[nzones];
//   zNodes = new int[4*nzones]; // ok because we know we have quads
//   zFaces = new int[4*nzones]; // ok because we know we have quads
   
   // make space for geometry info
   nodePos = NodalVectorField(&nodes);
   zoneVolume = ZonalScalarField(&zones);
   zonePos = ZonalVectorField(&zones);
   faceArea = FaceVectorField(&faces);

   // set up connectivity data, knowing it's logical orthogonal quads
   for (int i = 0; i < numZones(); i++)
   {
      int k = i % (kmax-1);
      int l = i / (kmax-1);

      nodeIndices[4*i+0] = k     + l*kmax;
      nodeIndices[4*i+1] = k + 1 + l*kmax;
      nodeIndices[4*i+2] = k + 1 + (l+1)*kmax;;
      nodeIndices[4*i+3] = k     + (l+1)*kmax;;

      faceIndices[4*i]   = 4*i;
      faceIndices[4*i+1] = 4*i+1;
      faceIndices[4*i+2] = 4*i+2;
      faceIndices[4*i+3] = 4*i+3;
   }
   zoneToNodes.bindRelationData( nodeIndices.data() );
   zoneToFaces.bindRelationData( faceIndices.data() );

   
   // set up initial mesh as regular orthonormal grid
   double dx = (xmax - xmin)/(kmax-1);
   double dy = (ymax - ymin)/(lmax-1);
   
   for (int l = 0; l < lmax; l++)
   {
      for (int k = 0; k < kmax; k++)
      {
         nodePos[k+l*kmax].x = xmin + k*dx;
         nodePos[k+l*kmax].y = ymin + l*dy;
      }
   }

   computeNewGeometry();

   //dumpMesh();
}



//----------------------------------------------
PolygonMeshXY::~PolygonMeshXY(void)
{
}

//----------------------------------------------
void PolygonMeshXY::computeNewGeometry(void)
{
   SLIC_ASSERT(numZones() > 0 && numNodes() > 2);
   
   VectorXY tmp; // we will need this
   
   // zone volumes and positions
   for (int iz = 0; iz < numZones(); iz++)
   {
      zoneVolume[iz] = 0.0;
      zonePos[iz] = VectorXY(0.0, 0.0);

      ZNodeSet zNodes = zoneToNodes[iz];
      ZFaceSet zFaces = zoneToFaces[iz];

      const int numZNodes = zNodes.size();
      for (int in = 0; in < numZNodes; in++)
      {
         // assumption is that zNodes lists nodes in counter-clockwise order around the zone
         IndexType nextInd =  (in+1) %  numZNodes;

         const VectorXY n0Pos = nodePos[ zNodes[in] ];
         const VectorXY n1Pos = nodePos[ zNodes[nextInd] ];

         zoneVolume[iz] += n0Pos.cross(n1Pos); // sum volume
         zonePos[iz] += n1Pos;  // sum position

         // Face areas are outward normal vectors
         tmp = n1Pos - n0Pos;
         faceArea[zFaces[in] ] = VectorXY(tmp.y, -tmp.x);
      }
      zoneVolume[iz] *= 0.5; // correct 2x volume from cross product
      // normalize position; if we did it by volume it would be more
      // correct (centroid), but arithmetic average usually works fine
      zonePos[iz] *= 1.0/numZNodes;
   }

   // check our final volumes are positive
   for (int z = 0; z < numZones(); z++)
       SLIC_ASSERT(zoneVolume[z] > 0.0);
}
//----------------------------------------------
void PolygonMeshXY::moveNodesToPosition(const VectorXY *newPos)
{
   for (int i=0; i < numNodes(); i++)
   {
      nodePos[i] = newPos[i];
   }
}

//----------------------------------------------
#include "myTimer.hpp"
#include <stdio.h>
int clock_getres(clockid_t clk_id, struct timespec *res);
int clock_gettime(clockid_t clk_id, struct timespec *tp);
int clock_settime(clockid_t clk_id, const struct timespec *tp);

//----------------------------------------------

VectorXY PolygonMeshXY::meshAverageKLZMemOrderA()
{
   timespec time1, time2; // TIMER

   clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);  // TIMER
   VectorXY ret;
   for (int i = 0; i < zones.size(); i++)
   {
      ret.accum(zonePos[i]);
   }
   ret.x = ret.x / zones.size();
   ret.y = ret.y / zones.size();

   clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); // TIMER
   timeElapsed = diffSeconds(time1, time2);
          
   return ret;
}
//----------------------------------------------



void PolygonMeshXY::dumpMesh()
{
    printf("Mesh has %i nodes and %i zones", nodes.size(), zones.size());

    printf("\n\nNodes");
    for(int i=0; i< nodes.size(); ++i)
    {
        VectorXY p = getPos(i);
        printf("\n\t Node %i -- pos (%g,%g) ", i, p.x, p.y);
    }

    printf("\n\nZones");
    for(int i=0; i< zones.size(); ++i)
    {
        VectorXY p = getZonePos(i);
        ZoneToNodeRelation::RelationSet zNodes = zoneToNodes[i];
        ZoneToNodeRelation::RelationSet zFaces = zoneToFaces[i];
        printf("\n\t Zone %i -- pos (%g,%g) -- vol %g -- zNumNodes %i -- zoneNodes %i %i %i %i -- zoneFaces %i %i %i %i"
                , i
                , p.x, p.y
                , zoneVol(i)
                , zNodes.size()
                , zNodes[0], zNodes[1], zNodes[2], zNodes[3]
                , zFaces[0], zFaces[1], zFaces[2], zFaces[3]
                );
    }

    printf("\n\n--\n\n");

}


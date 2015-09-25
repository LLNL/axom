#include "PolygonMeshXY.hpp"
#include <assert.h>
//----------------------------------------------
PolygonMeshXY::PolygonMeshXY(int kmax, int lmax,
                           double xmin, double xmax,
                           double ymin, double ymax)
{
   // This c'tor constructs a mesh that looks like a logically regular
   // mesh of quads, but in fact it is arbitrary polygons.  It's
   // intended for making quick test meshes and to do direct speed
   // comparisons between the polygon and regular quad mesh.
   
   assert(kmax > 0 and lmax > 0);
   
   nnodes = kmax*lmax;
   nzones = (kmax-1)*(lmax-1);

   // make space for connectivity data
   zNumNodes = new int[nzones];
   z2firstNode = new int[nzones];
   z2firstFace = new int[nzones];
   zNodes = new int[4*nzones]; // ok because we know we have quads
   zFaces = new int[4*nzones]; // ok because we know we have quads
   
   // make space for geometry info
   nodePos = new VectorXY[nnodes];
   zoneVolume = new double[nzones];
   zonePos = new VectorXY[nzones];
   faceArea = new VectorXY[4*nzones]; // ok because we know we have quads

   // set up connectivity data, knowing it's logical orthogonal quads
   for (int i = 0; i < nzones; i++)
   {
      zNumNodes[i] = 4;
      z2firstFace[i] = 4*i;
      z2firstNode[i] = 4*i;
      int k = i % (kmax-1);
      int l = i / (kmax-1);
      zNodes[4*i+0] = k     + l*kmax;
      zNodes[4*i+1] = k + 1 + l*kmax;
      zNodes[4*i+2] = k + 1 + (l+1)*kmax;;
      zNodes[4*i+3] = k     + (l+1)*kmax;;
      zFaces[4*i]   = 4*i; 
      zFaces[4*i+1] = 4*i+1;
      zFaces[4*i+2] = 4*i+2;
      zFaces[4*i+3] = 4*i+3;
   }
   
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
}

//----------------------------------------------
PolygonMeshXY::~PolygonMeshXY(void)
{
   delete [] z2firstFace;
   delete [] z2firstNode;
   delete [] zNodes;
   delete [] zFaces;
   delete [] zNumNodes;

   delete [] nodePos;
   delete [] zonePos;
   delete [] zoneVolume;
   delete [] faceArea;
}

//----------------------------------------------
void PolygonMeshXY::computeNewGeometry(void)
{
   assert(nzones > 0 and nnodes > 2);
   
   VectorXY tmp; // we will need this
   
   // zone volumes and positions
   for (int iz = 0; iz < nzones; iz++)
   {
      zoneVolume[iz] = 0.0;
      zonePos[iz] = VectorXY(0.0, 0.0);

      for (int in = 0; in < zNumNodes[iz]; in++)
      {
         // assumption is that zNodes lists nodes in counter-clockwise
         // order around the zone
         const VectorXY n0Pos = nodePos[zNodes[z2firstNode[iz]+in]];
         const VectorXY n1Pos = nodePos[zNodes[z2firstNode[iz]+(in+1) % zNumNodes[iz]]];
         zoneVolume[iz] += n0Pos.cross(n1Pos); // sum volume
         zonePos[iz] += n1Pos;  // sum position
         // Face areas are outward normal vectors
         tmp = n1Pos - n0Pos;
         faceArea[z2firstFace[iz] + in] = VectorXY(tmp.y, -tmp.x);
      }
      zoneVolume[iz] *= 0.5; // correct 2x volume from cross product
      // normalize position; if we did it by volume it would be more
      // correct (centroid), but arithmetic average usually works fine
      zonePos[iz] *= 1.0/zNumNodes[iz];
   }


   for (int z = 0; z < nzones; z++)
   {
       printf("Element %i -- vol %g -- pos (%g,%g) \n", z, zoneVolume[z], zonePos[z].x, zonePos[z].y);
   }


   // check our final volumes are positive
   for (int z = 0; z < nzones; z++) assert(zoneVolume[z] > 0.0);
}
//----------------------------------------------
void PolygonMeshXY::moveNodesToPosition(const VectorXY *newPos)
{
   for (int i=0; i < nnodes; i++)
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
   for (int i = 0; i < nzones; i++)
   {
      ret.accum(zonePos[i]);
   }
   ret.x = ret.x / nzones;
   ret.y = ret.y / nzones;

   clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); // TIMER
   timeElapsed = diffSeconds(time1, time2);
          
   return ret;
}
//----------------------------------------------

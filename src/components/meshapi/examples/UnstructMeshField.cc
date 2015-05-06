
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <set>
#include <vector>
#include <cmath>

typedef double DataType;

struct Point
{
   DataType x;
   DataType y;
   DataType z;

   Point(): x(0.0), y(0.0), z(0.0){}
};

int main(int argc, char** argv)
{
   //--------------------------------------------------------------

   typedef size_t IndexType;

   //--------------------------------------------------------------

   std::stringstream filePath;
   filePath << "../src/components/meshapi/examples/"
            << "ball_8.vtk";

   std::ifstream vtkMesh( filePath.str().c_str() );

   //--------------------------------------------------------------

   // Read some initial header stuff.  There's no error checking
   // like there should be, but this is a memory stressing routine,
   // not a robust mesh reader.

   std::string junk;
   while( junk != "POINTS" )
   {
      vtkMesh >> junk;
   }

   //--------------------------------------------------------------

   // Read in the POINT data, (that we'll call the nodes)
   IndexType numNodes;
   vtkMesh >> numNodes;
   vtkMesh >> junk;

   std::cout << "Number of nodes = " << numNodes << "\n";

   std::vector< Point > nodePosition( numNodes );
   for(IndexType i=0; i != numNodes; ++i)
   {
      vtkMesh >> nodePosition[i].x;
      vtkMesh >> nodePosition[i].y;
      vtkMesh >> nodePosition[i].z;
   }

   //--------------------------------------------------------------

   // Read in the CELL data, that we'll call zones.  We're going to assume
   // a VTK type 12 (Hexahedral zones).

   IndexType numZones;
   IndexType listSize;

   vtkMesh >> junk;
   vtkMesh >> numZones;
   vtkMesh >> listSize;
   IndexType numNodeZoneIndices = listSize - numZones;

   // This is only because we're assuming Hex's.  General meshes can be different.
   if( numZones * 8 != numNodeZoneIndices )
   {
      std::cout << "Error in reading mesh!\n";
      std::cout << "  numZones = " << numZones << "\n";
      std::cout << "  numZones*8 = " << numZones*8 << "\n";
      std::cout << "  numNodeZoneIndices = " << numNodeZoneIndices << "\n";
      return -1;
   }

   std::cout << "Number of zones = " << numZones << "\n";

   // This stores in one vector a list of nodes in each zone.
   std::vector< IndexType > nodesOfZone( numNodeZoneIndices );
   // This stores the index into nodesOfZone for the first node in a given
   // zone.  There's an extra one for how we'll actually use it later.
   // You could notice that this will be multiples of 8 into nodesOfZone,
   // and do an optimization, but won't here.
   typedef std::vector< IndexType >::iterator ZoneNodeIndex;
   std::vector< ZoneNodeIndex > startNodeOfZone( numZones + 1 );

   ZoneNodeIndex zni = nodesOfZone.begin();
   IndexType nodeCount;
   for(IndexType z=0; z != numZones; ++z )
   {
      vtkMesh >> nodeCount;
      if( nodeCount != 8 )
      {
         std::cout << "Unsupported mesh type with zone = " << z << ", nodeCount = " << nodeCount << "\n";
         return -1;
      }

      startNodeOfZone[z] = zni;
      for( IndexType n = 0; n != nodeCount; ++n )
      {
         vtkMesh >> *zni;
         ++zni;
      }
   }

   if( zni != nodesOfZone.end() )
   {
      std::cout << "Error reading nodes of zones!\n";
      std::cout << nodesOfZone.end() - zni << "\n";
      return -1;
   }

   // We specify one past the end to make iterating later easier.
   startNodeOfZone[numZones] = nodesOfZone.end();

   IndexType numNodesOfZone = nodesOfZone.size();
   std::cout << "numNodesOfZone = " << numNodesOfZone << "\n";

   //--------------------------------------------------------------

   // Close the file.
   vtkMesh.close();

   //--------------------------------------------------------------

   // Now build the zonesOfNode list.  This is tricky because we'll first
   // do it into a flexible data structure, then into an "optimized" one.

   std::vector< std::set< IndexType > > tmpZonesOfNode( numNodes );
   IndexType numZonesOfNode = 0;
   for(IndexType z=0; z != numZones; ++z)
   {
      for( ZoneNodeIndex zni = startNodeOfZone[z]; zni != startNodeOfZone[z+1]; ++zni )
      {
         tmpZonesOfNode[ *zni ].insert( z );
         ++numZonesOfNode;
      }
   }

   typedef std::vector< IndexType >::iterator NodeZoneIndex;
   std::vector< ZoneNodeIndex > startZoneOfNode( numNodes + 1 );
   std::vector< IndexType > zonesOfNode( numZonesOfNode );

   NodeZoneIndex nzi = zonesOfNode.begin();
   for(IndexType n=0; n != numNodes; ++n)
   {
      startZoneOfNode[n] = nzi;
      for( std::set<IndexType>::iterator z = tmpZonesOfNode[n].begin();
           z != tmpZonesOfNode[n].end();
           ++z )
      {
         *nzi = *z;
         ++nzi;
      }
   }
   startZoneOfNode[numNodes] = zonesOfNode.end();

   // Now free the temporary stuff. Yes, ugly, but what C++ wanted.
   std::vector< std::set< IndexType > >().swap( tmpZonesOfNode );

   std::cout << "numZonesOfNode = " << numZonesOfNode << "\n";
   if( nzi != zonesOfNode.end() )
   {
      std::cout << "Error creating zones of Node list!\n";
      return -1;
   }

   //--------------------------------------------------------------
   //--------------------------------------------------------------

   // Now that we have the mesh in memory, we can start to do things with it.

   // Compute the zone position, or the average of the nodes around the zone
   std::vector< Point > zonePosition( numZones );
   for( IndexType z = 0; z != numZones; ++z )
   {
      // Instead of counting, we could subtract two consecutive startNodeOfZones.
      int count = 0;
      Point& zp = zonePosition[z];
      for( ZoneNodeIndex n = startNodeOfZone[z];
           n != startNodeOfZone[z+1];
           ++n)
      {
         zp.x += nodePosition[*n].x;
         zp.y += nodePosition[*n].y;
         zp.z += nodePosition[*n].z;
         ++count;
      }
      zp.x /= count;
      zp.y /= count;
      zp.z /= count;
   }

   //--------------------------------------------------------------

   // Compute a zone field based on the position
   std::vector< DataType > zoneField( numZones );
   for( IndexType z = 0; z != numZones; ++z )
   {
      const Point& zp = zonePosition[z];
      // What's the radius?
      zoneField[z] = std::sqrt( zp.x*zp.x + zp.y*zp.y + zp.z*zp.z);
   }

   //--------------------------------------------------------------

   // Compute the node average version, and the maximum relative error
   // Relying on zeroing out the average field by the vector constructor.
   std::vector< DataType > nodeFieldAvg( numNodes, 0.0 );
   std::vector< DataType > nodeFieldExact( numNodes );
   double errSqSum = 0.0;
   for( IndexType n = 0; n != numNodes; ++n )
   {
      const Point& np = nodePosition[n];
      // What's the radius?
      nodeFieldExact[n] = std::sqrt( np.x*np.x + np.y*np.y + np.z*np.z);

      IndexType numZones = startZoneOfNode[n+1]-startZoneOfNode[n];
      for( NodeZoneIndex z = startZoneOfNode[n];
           z != startZoneOfNode[n+1];
           ++z)
      {
         //std::cout << "node, zone = " << n << " " <<  *z << "\n";
         nodeFieldAvg[n] += zoneField[*z];
      }
      nodeFieldAvg[n] /= numZones;
      const double err = nodeFieldAvg[n] - nodeFieldExact[n];
      errSqSum += err*err;
   }

   std::cout << "The L2-ish error in the node average radius was " << std::sqrt( errSqSum / numNodes ) << "\n";

   //--------------------------------------------------------------
   return 0;
}

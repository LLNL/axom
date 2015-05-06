/**
 * \file UnstructMeshField.cc
 *
 * \brief Simple user of the mesh api class.  Loads a hex mesh from a VTK file, generates the Node to Zone relation and does simple mesh processing
 *
 * \author T. Brunner (original)
 * \author K. Weiss (modified to use mesh API)
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <set>
#include <vector>
#include <cmath>
#include <cstdlib>

#include <cstdio>  /* defines FILENAME_MAX */
#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif


namespace asctoolkit {
namespace meshapi {
namespace examples {


typedef double DataType;

struct Point
{
   DataType x;
   DataType y;
   DataType z;

   Point(): x(0.0), y(0.0), z(0.0){}
};

typedef size_t IndexType;

struct HexMesh {

    enum { NODES_PER_ZONE = 8 };

    typedef std::vector< Point >                    PositionsVec;
    typedef std::vector< DataType >                 NodeField;
    typedef std::vector< DataType >                 ZoneField;

    typedef std::vector< IndexType >                NodesOfZone;
    typedef NodesOfZone::iterator                   ZoneNodeIndex;
    typedef std::vector< ZoneNodeIndex >            StartNodeOfZones;

    typedef std::vector< IndexType >                ZonesOfNode;
    typedef ZonesOfNode::iterator                   NodeZoneIndex;
    typedef std::vector< NodeZoneIndex >            StartZoneOfNodes;

    /** \brief Simple accessor for the number of nodes in the mesh  */
    IndexType  numNodes() const { return nodePosition.size(); }

    /** \brief Simple accessor for the number of zones in the mesh */
    IndexType  numZones() const { return startNodeOfZone.size() -1; }


    // Node-centered position field
    PositionsVec                                    nodePosition;
    PositionsVec                                    zonePosition;

    // storage for relation_(3,0) -- zones -> nodes
    NodesOfZone                                     nodesOfZone;
    StartNodeOfZones                                startNodeOfZone;

    // storage for relation_(0,3) -- nodes -> zones
    ZonesOfNode                                     zonesOfNode;
    StartZoneOfNodes                                startZoneOfNode;

    // A zone field and two node fields
    ZoneField                                       zoneField;
    NodeField                                       nodeFieldExact;
    NodeField                                       nodeFieldAvg;
};

void printCWD()
{
    char cCurrentPath[FILENAME_MAX];

    if (!GetCurrentDir(cCurrentPath, FILENAME_MAX))
    {
        std::abort();
    }

    std::cout << "The current working directory is " << std::string(cCurrentPath);

}

HexMesh readHexMesh(std::string fileName)
{
    HexMesh mesh;

    std::ifstream vtkMesh( fileName.c_str() );

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

    mesh.nodePosition = HexMesh::PositionsVec( numNodes );
    for(IndexType i=0; i != numNodes; ++i)
    {
       vtkMesh >> mesh.nodePosition[i].x;
       vtkMesh >> mesh.nodePosition[i].y;
       vtkMesh >> mesh.nodePosition[i].z;
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
    if( numZones * (HexMesh::NODES_PER_ZONE) != numNodeZoneIndices )
    {
       std::cerr << "Error in reading mesh!\n";
       std::cerr << "  numZones = " << numZones << "\n";
       std::cerr << "  numZones*8 = " << numZones* (HexMesh::NODES_PER_ZONE) << "\n";
       std::cerr << "  numNodeZoneIndices = " << numNodeZoneIndices << "\n";
       std::abort();
    }

    std::cout << "Number of zones = " << numZones << "\n";

    // This stores in one vector a list of nodes in each zone.
    mesh.nodesOfZone = HexMesh::NodesOfZone( numNodeZoneIndices );

    // This stores the index into nodesOfZone for the first node in a given
    // zone.  There's an extra one for how we'll actually use it later.
    // You could notice that this will be multiples of 8 into nodesOfZone,
    // and do an optimization, but won't here.

    mesh.startNodeOfZone = HexMesh::StartNodeOfZones( numZones + 1 );

    HexMesh::ZoneNodeIndex zni = mesh.nodesOfZone.begin();
    IndexType nodeCount;
    for(IndexType z=0; z != numZones; ++z )
    {
       vtkMesh >> nodeCount;
       if( nodeCount != 8 )
       {
          std::cerr << "Unsupported mesh type with zone = " << z << ", nodeCount = " << nodeCount << "\n";
          std::abort();
       }

       mesh.startNodeOfZone[z] = zni;
       for( IndexType n = 0; n != nodeCount; ++n )
       {
          vtkMesh >> *zni;
          ++zni;
       }
    }

    if( zni != mesh.nodesOfZone.end() )
    {
       std::cerr << "Error reading nodes of zones!\n";
       std::cerr << mesh.nodesOfZone.end() - zni << "\n";
       std::abort();
    }

    // We specify one past the end to make iterating later easier.
    mesh.startNodeOfZone[numZones] = mesh.nodesOfZone.end();

    IndexType numNodesOfZone = mesh.nodesOfZone.size();
    std::cout << "numNodesOfZone = " << numNodesOfZone << "\n";

    //--------------------------------------------------------------

    // Close the file.
    vtkMesh.close();

    return mesh;
}

void generateNodeZoneRelation(HexMesh & mesh)
{
    // Now build the zonesOfNode list.
    // This is tricky because we'll first
    // do it into a flexible data structure, then into an "optimized" one.

    typedef std::set< IndexType> NodeIndexSet;
    typedef NodeIndexSet::iterator NodeIndexSetIterator;


    std::vector< NodeIndexSet > tmpZonesOfNode( mesh.numNodes() );
    IndexType numZonesOfNode = 0;
    for(IndexType z=0; z != mesh.numZones(); ++z)
    {
      for( HexMesh::ZoneNodeIndex zni = mesh.startNodeOfZone[z]; zni != mesh.startNodeOfZone[z+1]; ++zni )
      {
         tmpZonesOfNode[ *zni ].insert( z );
         ++numZonesOfNode;
      }
    }

    // Create the static-variable node->zone relation
    mesh.startZoneOfNode = HexMesh::StartZoneOfNodes( mesh.numNodes() + 1 );
    mesh.zonesOfNode = HexMesh::ZonesOfNode( numZonesOfNode );
    HexMesh::NodeZoneIndex nzi = mesh.zonesOfNode.begin();
    for(IndexType n=0; n != mesh.numNodes(); ++n)
    {
      mesh.startZoneOfNode[n] = nzi;
      for( NodeIndexSetIterator z = tmpZonesOfNode[n].begin(); z != tmpZonesOfNode[n].end(); ++z )
      {
         *nzi = *z;
         ++nzi;
      }
    }
    mesh.startZoneOfNode[mesh.numNodes()] = mesh.zonesOfNode.end();

    // Now free the temporary stuff. Yes, ugly, but what C++ wanted.
    std::vector< NodeIndexSet >().swap( tmpZonesOfNode );

    std::cout << "\n numZonesOfNode = " << numZonesOfNode << "\n";
    if( nzi != mesh.zonesOfNode.end() )
    {
      std::cerr << "Error creating zones of Node list!\n";
      std::abort();
    }

}

void computeZoneBarycenters(HexMesh & mesh)
{
    // Compute the zone position, or the average of the nodes around the zone
    mesh.zonePosition = HexMesh::PositionsVec( mesh.numZones() );
    for( IndexType z = 0; z != mesh.numZones(); ++z )
    {
       // Instead of counting, we could subtract two consecutive startNodeOfZones.
       int count = 0;
       Point& zp = mesh.zonePosition[z];
       for( HexMesh::ZoneNodeIndex n = mesh.startNodeOfZone[z]; n != mesh.startNodeOfZone[z+1]; ++n)
       {
          zp.x += mesh.nodePosition[*n].x;
          zp.y += mesh.nodePosition[*n].y;
          zp.z += mesh.nodePosition[*n].z;
          ++count;
       }
       zp.x /= count;
       zp.y /= count;
       zp.z /= count;
    }
}

void createZoneRadiusField (HexMesh & mesh)
{
    // Compute a zone field based on the position
    mesh.zoneField = HexMesh::ZoneField ( mesh.numZones() );
    for (IndexType z = 0; z != mesh.numZones(); ++z)
    {
        const Point& zp = mesh.zonePosition[z];
        // What's the radius?
        mesh.zoneField[z] = std::sqrt (zp.x * zp.x + zp.y * zp.y + zp.z * zp.z);
    }
}

void computeNodalErrors(HexMesh & mesh)
{
    // Compute the node average version, and the maximum relative error
    // Relying on zeroing out the average field by the vector constructor.
    mesh.nodeFieldAvg = HexMesh::NodeField( mesh.numNodes(), 0.0 );
    mesh.nodeFieldExact = HexMesh::NodeField( mesh.numNodes() );
    double errSqSum = 0.0;
    for( IndexType n = 0; n != mesh.numNodes(); ++n )
    {
       const Point& np = mesh.nodePosition[n];
       // What's the radius?
       mesh.nodeFieldExact[n] = std::sqrt( np.x*np.x + np.y*np.y + np.z*np.z);

       IndexType numZones = mesh.startZoneOfNode[n+1]- mesh.startZoneOfNode[n];
       for( HexMesh::NodeZoneIndex z = mesh.startZoneOfNode[n]; z != mesh.startZoneOfNode[n+1]; ++z)
       {
          //std::cout << "node, zone = " << n << " " <<  *z << "\n";
          mesh.nodeFieldAvg[n] += mesh.zoneField[*z];
       }
       mesh.nodeFieldAvg[n] /= numZones;
       const double err = mesh.nodeFieldAvg[n] - mesh.nodeFieldExact[n];
       errSqSum += err*err;
    }

    std::cout << "\n\tThe L2-ish error in the node average radius was " << std::sqrt( errSqSum / mesh.numNodes() ) << "\n";
}

}   // end namespace examples
}   // end namespace meshapi
}   // end namespace asctoolkit


int main()
{
   using namespace asctoolkit::meshapi::examples;

   //--------------------------------------------------------------
   // Load the hexmesh from the vtk file
   //--------------------------------------------------------------

   printCWD();

   int const NUM_RESOLUTIONS = 4;
   int fileResolutions[] = {1,2,4,8};

   for(int res = 0; res < NUM_RESOLUTIONS; ++res)
   {
       std::stringstream filePath;
       filePath << "../src/components/meshapi/examples/"
                 << "ball_"<< fileResolutions[res] << ".vtk";
       std::string meshName = filePath.str();

       std::cout<<"\n** Loading mesh file '" << meshName << "' and generating zone-> node relation...\n";

       HexMesh hexMesh = readHexMesh( meshName );

       //--------------------------------------------------------------

       // Now build the node to zone relation
       std::cout<<"\n** Generating node->zone relation...";
       generateNodeZoneRelation(hexMesh);

       //--------------------------------------------------------------

       // Now that we have the mesh in memory, we can start to do things with it.

       std::cout<<"\n** Computing zone barycenters using zone->node relation...";
       computeZoneBarycenters(hexMesh);

       std::cout<<"\n** Generating a zone-centered radius field...";
       createZoneRadiusField(hexMesh);

       std::cout<<"\n** Computing node-based errors using node->zone relation...";
       computeNodalErrors(hexMesh);

       std::cout<<"\ndone." << std::endl;
   }

   //--------------------------------------------------------------
   return 0;
}




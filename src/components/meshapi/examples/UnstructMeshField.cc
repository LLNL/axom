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

#include "meshapi/FileUtilities.hpp"

#include "meshapi/OrderedSet.hpp"
#include "meshapi/StaticVariableRelation.hpp"




namespace asctoolkit {
namespace meshapi {
namespace examples {


typedef size_t IndexType;
typedef double DataType;


struct Point
{
   Point() : x(DataType()), y(DataType()), z(DataType()){}

   DataType radius() const { return std::sqrt( x*x + y*y + z*z); }

   DataType x;
   DataType y;
   DataType z;
};

struct HexMesh {

public:
    enum { NODES_PER_ZONE = 8 };

    // types for sets
    typedef asctoolkit::meshapi::OrderedSet NodeSet;
    typedef asctoolkit::meshapi::OrderedSet ZoneSet;

    // types for relations
    typedef asctoolkit::meshapi::StaticVariableRelation NodeZoneRelation;
    typedef NodeZoneRelation::RelationVecConstIterator  NodeZoneIterator;

    typedef asctoolkit::meshapi::StaticVariableRelation ZoneNodeRelation;
    typedef ZoneNodeRelation::RelationVecConstIterator  ZoneNodeIterator;

    typedef NodeZoneRelation::Index         IndexType;
    typedef NodeZoneRelation::size_type     SizeType;

    // types for maps
    // TODO: Convert to meshapi::Map
    typedef std::vector< Point >                    PositionsVec;
    typedef std::vector< DataType >                 NodeField;
    typedef std::vector< DataType >                 ZoneField;

public:
    /** \brief Simple accessor for the number of nodes in the mesh  */
    SizeType  numNodes() const { return nodes.size(); }

    /** \brief Simple accessor for the number of zones in the mesh */
    SizeType  numZones() const { return zones.size(); }

public:
    /// Sets in the mesh
    NodeSet nodes;
    ZoneSet zones;

    /// Relations in the mesh
    ZoneNodeRelation    relationZoneNode;   // storage for relation_(3,0) -- zones -> nodes
    NodeZoneRelation    relationNodeZone;   // storage for relation_(0,3) -- nodes -> zones


    /// Maps (fields) in the mesh -- positions and scalar fields

    // Node-centered position field
    PositionsVec                                    nodePosition;
    PositionsVec                                    zonePosition;

    // A zone field and two node fields
    ZoneField                                       zoneField;
    NodeField                                       nodeFieldExact;
    NodeField                                       nodeFieldAvg;
};


void readHexMesh(std::string fileName, HexMesh& mesh)
{
    std::ifstream vtkMesh( fileName.c_str() );
    if(!vtkMesh)
    {
        std::cerr <<"fstream error -- problem opening file: '" << fileName <<"'"
                  << "\nThe current working directory is: '" << asctoolkit::meshapi::util::getCWD() <<"'"
                  << std::endl;
        std::abort();
    }

    //--------------------------------------------------------------

    // Read some initial header stuff.  Note: this is not a robust vtkreader
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

    // Create the set of Nodes
    mesh.nodes = HexMesh::NodeSet(numNodes);

    // Create the nodal position field, and read in file file
    mesh.nodePosition = HexMesh::PositionsVec( numNodes );
    for(IndexType i=0; i != numNodes; ++i)
    {
       vtkMesh >> mesh.nodePosition[i].x;
       vtkMesh >> mesh.nodePosition[i].y;
       vtkMesh >> mesh.nodePosition[i].z;
    }

    //--------------------------------------------------------------

    // Read in the CELL data, that we'll call zones.  We're going to assume hexahedra (VTK type 12)

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

    // Create the set of Zones
    mesh.zones = HexMesh::ZoneSet(numZones);

    // Create the topological incidence relation from zones to nodes
    mesh.relationZoneNode = HexMesh::ZoneNodeRelation(&mesh.zones, &mesh.nodes);

    // Read in and encode this as a static variable relation
    // TODO: Replace with a static constant relation when that class is developed
    typedef HexMesh::ZoneNodeRelation::RelationVec  RelationVec;
    typedef RelationVec::iterator                   RelationVecIterator;

    // Setup the 'begins' vector
    //  -- exploit the fact that the relation is constant
    //  -- note that for a constant relation, this array is not really necessary
    RelationVec beginsVec( mesh.zones.size() + 1 );
    for(HexMesh::IndexType idx=0; idx <= mesh.zones.size(); ++idx)
    {
        beginsVec[idx] = idx * HexMesh::NODES_PER_ZONE;
    }

    // Setup the 'offsets' vector
    RelationVec offsetsVec ( numNodeZoneIndices );
    RelationVecIterator oIt = offsetsVec.begin();
    IndexType nodeCount;
    typedef HexMesh::ZoneSet::iterator SetIterator;
    for(SetIterator zIt = mesh.zones.begin(); zIt != mesh.zones.end(); ++zIt)
    {
       vtkMesh >> nodeCount;
       if( nodeCount != HexMesh::NODES_PER_ZONE )
       {
          std::cerr << "Unsupported mesh type with zone = " << *zIt
                    << ", nodeCount = " << nodeCount
                    << " (expected " << HexMesh::NODES_PER_ZONE << ")\n";
          std::abort();
       }

       for( IndexType n = 0; n != nodeCount; ++n )
       {
          vtkMesh >> *oIt++;
       }
    }

    // Close the file.
    vtkMesh.close();

    // Set the relation here by copying over the data buffers
    // NOTE: This should be a lot cleaner once we hook up to the datastore
    mesh.relationZoneNode.setRelation(beginsVec, offsetsVec);

    // Check that the relation is valid
    mesh.relationZoneNode.isValid();

    //--------------------------------------------------------------

    if( oIt != offsetsVec.end() )
    {
       std::cerr << "Error reading nodes of zones!\n";
       std::cerr << offsetsVec.end() - oIt << "\n";
       std::abort();
    }

}

void generateNodeZoneRelation(HexMesh & mesh)
{
    // Now build the zonesOfNode list.

    // We will first build this as a vector of sets
    // and then linearize this into a static variable relation
    // TODO: Replace this with a dynamic variable relation once it is developed

    typedef std::set< HexMesh::IndexType> NodeIndexSet;
    typedef NodeIndexSet::iterator NodeIndexSetIterator;

    std::vector< NodeIndexSet > tmpZonesOfNode( mesh.numNodes() );
    IndexType numZonesOfNode = 0;

    typedef HexMesh::ZoneSet::iterator ZoneIterator;
    for(ZoneIterator zIt = mesh.zones.begin(); zIt != mesh.zones.end(); ++zIt)
    {
        typedef HexMesh::ZoneNodeRelation::RelationVecConstIterator RelVecIt;
        RelVecIt znIt = mesh.relationZoneNode.begin(*zIt);
        RelVecIt znEnd = mesh.relationZoneNode.end(*zIt);
        for( ; znIt != znEnd; ++znIt )
        {
            tmpZonesOfNode[ *znIt ].insert( *zIt );
            ++numZonesOfNode;
        }
    }

    // Create the relation here
    mesh.relationNodeZone = HexMesh::NodeZoneRelation( &mesh.nodes, &mesh.zones);

    // Now, linearize the dynamic relation into a static relation here
    typedef HexMesh::NodeZoneRelation::RelationVec  RelationVec;
    typedef RelationVec::iterator                   RelationVecIterator;

    RelationVec beginsVec( mesh.nodes.size() + 1 );
    RelationVec offsetsVec( numZonesOfNode );
    HexMesh::IndexType count = 0;
    for(HexMesh::NodeSet::iterator nIt = mesh.nodes.begin(), nEnd = mesh.nodes.end(); nIt != nEnd; ++nIt)
    {
        beginsVec[*nIt] = count;
        NodeIndexSetIterator sIt = tmpZonesOfNode[*nIt].begin();
        NodeIndexSetIterator sEnd = tmpZonesOfNode[*nIt].end();
        for(; sIt != sEnd; ++sIt)
        {
            offsetsVec[count++] = *sIt;
        }
    }
    beginsVec[mesh.nodes.size()] = count;

    // Send the data buffers over to the relation now and check the validity
    mesh.relationNodeZone.setRelation(beginsVec, offsetsVec);
    mesh.relationNodeZone.isValid();

    std::cout << "\n numZonesOfNode = " << numZonesOfNode << "\n";
    if( count != numZonesOfNode )
    {
      std::cerr << "Error creating zones of Node list!\n";
      std::abort();
    }
}

void computeZoneBarycenters(HexMesh & mesh)
{
    typedef HexMesh::ZoneSet::iterator ZoneIter;
    typedef HexMesh::ZoneNodeIterator ZNIterator;

    // Compute the zone positions as the the averages of the positions of the nodes around each zone
    mesh.zonePosition = HexMesh::PositionsVec( mesh.numZones() );

    // Outer loop over each zone in the set
    for(ZoneIter zIt = mesh.zones.begin(); zIt != mesh.zones.end(); ++zIt )
    {
       Point& zonePos = mesh.zonePosition[ *zIt ];

       // Inner loop over each node of the zone
       ZNIterator  nIt  = mesh.relationZoneNode.begin(*zIt);
       ZNIterator  nEnd = mesh.relationZoneNode.end(*zIt);
       DataType numNodes = static_cast<DataType>(mesh.relationZoneNode.size(*zIt));
       for( ; nIt != nEnd; ++nIt)
       {
           Point& nodePos = mesh.nodePosition[*nIt];
           zonePos.x += nodePos.x;
           zonePos.y += nodePos.y;
           zonePos.z += nodePos.z;
       }
       zonePos.x /= numNodes;
       zonePos.y /= numNodes;
       zonePos.z /= numNodes;
    }
}

void createZoneRadiusField (HexMesh & mesh)
{
    // Compute a zone field based on the L2 norm of their position vectors
    mesh.zoneField = HexMesh::ZoneField ( mesh.numZones() );

    typedef HexMesh::ZoneSet::iterator ZoneIter;
    for (ZoneIter zIt = mesh.zones.begin(); zIt != mesh.zones.end();++zIt)
    {
        Point const& zp = mesh.zonePosition[*zIt];
        // What's the radius?
        mesh.zoneField[*zIt] = zp.radius();
    }
}

void computeNodalErrors(HexMesh & mesh)
{
    // Compute the node average version, and the maximum relative error
    typedef HexMesh::NodeSet::iterator NodeIter;
    typedef HexMesh::NodeZoneIterator NZIterator;

    // Relying on zeroing out the average field by the vector constructor.
    mesh.nodeFieldAvg = HexMesh::NodeField( mesh.numNodes(), 0.0 );
    mesh.nodeFieldExact = HexMesh::NodeField( mesh.numNodes() );
    double errSqSum = 0.0;

    // Outer loop over each node
    for (NodeIter nIt = mesh.nodes.begin(); nIt != mesh.nodes.end();++nIt)
    {
       const Point& np = mesh.nodePosition[*nIt];

       // What's the radius?
       mesh.nodeFieldExact[*nIt] = np.radius();

       // Inner loop over each zone of the node
       NZIterator  zIt  = mesh.relationNodeZone.begin(*nIt);
       NZIterator  zEnd = mesh.relationNodeZone.end(*nIt);
       DataType numZones = static_cast<DataType>(mesh.relationNodeZone.size(*nIt));
       for( ; zIt != zEnd; ++zIt)
       {
           mesh.nodeFieldAvg[*nIt] += mesh.zoneField[*zIt];
       }

       mesh.nodeFieldAvg[*nIt] /= numZones;

       double const err = mesh.nodeFieldAvg[*nIt] - mesh.nodeFieldExact[*nIt];
       errSqSum += err*err;
    }

    std::cout   << "\n\tThe L2-ish error in the node average radius was "
                << std::sqrt( errSqSum / mesh.numNodes() )
                << std::endl;
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



#ifndef USE_ONE
   int const NUM_RESOLUTIONS = 4;
   int fileResolutions[] = {1,2,4,8};
#else
   int const NUM_RESOLUTIONS = 1;
   int fileResolutions[] = {1};
#endif

   for(int res = 0; res < NUM_RESOLUTIONS; ++res)
   {
       std::stringstream filePath;
       filePath  << "../src/components/meshapi/examples/"
                 << "ball_"<< fileResolutions[res] << ".vtk";
       std::string meshName = filePath.str();

       std::cout<<"\n** Loading mesh file '" << meshName << "' and generating zone-> node relation...\n";

       HexMesh hexMesh;
       readHexMesh( meshName, hexMesh );

       //--------------------------------------------------------------

       // Now build the node to zone relation
       std::cout<<"\n** Generating node->zone relation...";
       generateNodeZoneRelation( hexMesh );

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


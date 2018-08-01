/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * \file
 *
 * \brief Simple example that uses Slam for generating and processing a simple
 *  3D mesh.
 *
 * \details Loads a hex mesh from a VTK file, generates the Node to Zone
 *  relation and does simple mesh processing.
 *
 * \author T. Brunner (original)
 * \author K. Weiss (modified to use axom's Slam component)
 */

#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/utilities/FileUtilities.hpp"

#include "axom/slic/interface/slic.hpp"
#include "axom/slic/core/UnitTestLogger.hpp"

#include "axom/slam/policies/IndirectionPolicies.hpp"
#include "axom/slam/policies/CardinalityPolicies.hpp"
#include "axom/slam/Set.hpp"
#include "axom/slam/RangeSet.hpp"
#include "axom/slam/StaticRelation.hpp"

#include "axom/slam/DynamicVariableRelation.hpp"

#include "axom/slam/Map.hpp"
#include "axom/slam/Utilities.hpp"
#include "axom/slam/FieldRegistry.hpp"

#include "fmt/fmt.hpp"


namespace slamUnstructuredHex
{

typedef axom::slam::MeshIndexType IndexType;
typedef double DataType;


/** Simple point class for this example */
struct Point
{
  Point(const DataType& x, const DataType& y, const DataType& z)
    : m_x(x), m_y(y), m_z(z){}
  Point() : m_x(DataType()), m_y(DataType()), m_z(DataType()){}

  DataType        radius() const
  {
    return std::sqrt( m_x * m_x + m_y * m_y + m_z * m_z);
  }

  Point& operator +=(const Point& pt)
  {
    m_x += pt.m_x;
    m_y += pt.m_y;
    m_z += pt.m_z;
    return *this;
  }

  Point& operator *=(const DataType& sc)
  {
    m_x *= sc;
    m_y *= sc;
    m_z *= sc;
    return *this;
  }

  template<typename T>
  Point& operator /=(const T& sc)
  {
    return operator*=( 1. / sc);
  }

  DataType m_x, m_y, m_z;
};

/// Some operations on Points
Point operator  +(const Point& pt1, const Point& pt2)
{
  Point pt(pt1);
  pt += pt2;
  return pt;
}

Point operator  *(const Point& pt1, const DataType& sc)
{
  Point pt(pt1);
  pt *= sc;
  return pt;
}

Point operator  *(const DataType& sc,const Point& pt1)
{
  Point pt(pt1);
  pt *= sc;
  return pt;
}

template<typename T>
Point operator  /(const Point& pt1, const T& sc)
{
  Point pt(pt1);
  pt *= (1. / sc);
  return pt;
}


/**
 * \brief Simple hex mesh for this example.
 *  Contains the necessary types, for a set of nodes and zones,
 *  the incidence relations between them and some fields defined on them.
 */
struct HexMesh
{

public:
  enum
  {
    COORDS_PER_NODE = 3,
    NODES_PER_ZONE = 8
  };

  /// types for sets
  typedef axom::slam::PositionSet NodeSet;
  typedef axom::slam::PositionSet ZoneSet;
  typedef ZoneSet::PositionType PositionType;
  typedef ZoneSet::IndexType IndexType;

  /// types for relations
  typedef axom::slam::policies::
    STLVectorIndirection<PositionType, PositionType>  STLIndirection;
  typedef axom::slam::policies::
    VariableCardinality<PositionType, STLIndirection> VariableCardinality;
  typedef axom::slam::
    StaticRelation<VariableCardinality,STLIndirection,NodeSet,ZoneSet>
    NodeToZoneRelation;
  typedef NodeToZoneRelation::RelationConstIterator NodeZoneIterator;

  typedef axom::slam::policies::
    CompileTimeStride<PositionType, NODES_PER_ZONE>   ZNStride;
  typedef axom::slam::policies::
    ConstantCardinality<PositionType, ZNStride>       ConstantCardinality;
  typedef axom::slam::
    StaticRelation<ConstantCardinality,STLIndirection,ZoneSet,NodeSet>
    ZoneToNodeRelation;
  typedef ZoneToNodeRelation::RelationConstIterator ZoneNodeIterator;

  /// types for maps
  typedef axom::slam::Map< Point >                      PositionsVec;
  typedef axom::slam::Map< DataType >                   NodeField;
  typedef axom::slam::Map< DataType >                   ZoneField;

public:
  /** \brief Simple accessor for the number of nodes in the mesh  */
  PositionType  numNodes() const { return nodes.size(); }

  /** \brief Simple accessor for the number of zones in the mesh */
  PositionType  numZones() const { return zones.size(); }

public:
  /// Sets in the mesh
  NodeSet nodes;
  ZoneSet zones;

  /// Relations in the mesh
  ZoneToNodeRelation zoneToNodeRelation;        // storage for relation_(3,0) --
                                                // zones -> nodes
  NodeToZoneRelation nodeToZoneRelation;        // storage for relation_(0,3) --
                                                // nodes -> zones


  /// Maps (fields) defined on the mesh -- nodal and zonal positions and scalar
  // fields
  PositionsVec nodePosition;
  PositionsVec zonePosition;
  ZoneField zoneField;
  NodeField nodeFieldExact;
  NodeField nodeFieldAvg;
};

/// The repository is a proxy for a data allocator/manager
struct Repository
{
  // Define the explicit instances of our local (key/value) datastore for int
  // and double
  typedef axom::slam::FieldRegistry<int>    IntsRegistry;
  typedef axom::slam::FieldRegistry<double> RealsRegistry;
  typedef axom::slam::Map<int>              IntField;
  typedef axom::slam::Map<double>           RealField;

  static IntsRegistry intsRegistry;
  static RealsRegistry realsRegistry;
};

Repository::IntsRegistry Repository::intsRegistry;
Repository::RealsRegistry Repository::realsRegistry;


/** A simple class to read a VTK hex mesh */
class SimpleVTKHexMeshReader
{
public:
  // uses RAII to open/close the file
  SimpleVTKHexMeshReader(const std::string & fileName)
    : vtkMesh( fileName.c_str() )
  {
    if(!vtkMesh)
    {
      using namespace axom::slam::util;
      std::string ancesFile = findFileInAncestorDirs( fileName);
      SLIC_ERROR_IF( !axom::utilities::filesystem::pathExists( ancesFile),
                     fmt::format(
                       "Tried opening file '{}', but it does not exist.",
                       ancesFile) );

      SLIC_INFO("Opening file " << ancesFile);
      vtkMesh.open( ancesFile.c_str() );
    }

  }
  ~SimpleVTKHexMeshReader()
  {
    vtkMesh.close();        // Close the file.
  }

  void parseMeshFile()
  {
    typedef Repository::RealsRegistry::BufferType RealBuf;
    typedef Repository::IntsRegistry::BufferType IndexBuf;

    // Read some initial header stuff.  Note: this is not a robust vtkreader
    std::string junk;

    IndexType numNodes;

    /// Read in POINT data, (which we'll call the nodes)
    while( junk != "POINTS" ) { vtkMesh >> junk; }
    vtkMesh >> numNodes >> junk;

    Repository::intsRegistry.addScalar("num_nodes", numNodes);
    SLIC_INFO("-- Number of nodes: " << numNodes);

    const IndexType numCoords = HexMesh::COORDS_PER_NODE * numNodes;
    RealBuf& pointData = Repository::realsRegistry.addBuffer("node_positions",
                                                             numCoords);
    for(IndexType idx = 0 ; idx < numCoords ; ++idx)
    {
      vtkMesh >> pointData[idx];
    }

    /// Read in the CELL data (which we'll call zones).
    /// We're going to assume hexahedra (VTK type 12)
    IndexType numZones, listSize, nodeCount;
    vtkMesh >> junk >> numZones >> listSize;
    Repository::intsRegistry.addScalar("num_zones", numZones);
    const IndexType numNodeZoneIndices = HexMesh::NODES_PER_ZONE  * numZones;

    SLIC_INFO("-- Number of zones: " << numZones );

    // Note: The VTK format has an extra value per zone for the number of
    // indices
    // This is constant since we're assuming a Hex mesh.  General meshes can be
    // different.
    SLIC_ASSERT_MSG( (listSize - numZones) == numNodeZoneIndices,
                     fmt::format(
                       "Error while reading mesh!\n "
                       "numZones = {0}; numZones*{1} = {2}; indices in file = {3}",
                       numZones, static_cast<int>(HexMesh::NODES_PER_ZONE),
                       numNodeZoneIndices, listSize - numZones));

    IndexBuf& zn_indices = Repository::intsRegistry
                           .addBuffer("zone_node_indices", numNodeZoneIndices );
    IndexType idx = 0;
    for(IndexType zoneIdx = 0 ; zoneIdx < numZones ; ++zoneIdx)
    {
      vtkMesh >> nodeCount;
      SLIC_ASSERT( nodeCount == HexMesh::NODES_PER_ZONE);

      for( IndexType n = 0 ; n < (HexMesh::NODES_PER_ZONE) ; ++n )
      {
        vtkMesh >> zn_indices[idx++];
      }
    }
  }

private:
  std::ifstream vtkMesh;
};

void readHexMesh(std::string fileName, HexMesh* mesh)
{
  {
    SimpleVTKHexMeshReader vtkMeshReader(fileName);
    vtkMeshReader.parseMeshFile();
  }
  typedef Repository::RealsRegistry::BufferType RealBuf;
  typedef Repository::IntsRegistry::BufferType IndexBuf;

  /// Check that the mesh has been loaded properly
  if(   !(Repository::intsRegistry.hasScalar("num_nodes")
          && Repository::intsRegistry.hasScalar("num_zones")
          && Repository::realsRegistry.hasBuffer("node_positions")
          && Repository::intsRegistry.hasBuffer("zone_node_indices") ))
  {
    SLIC_ERROR("Hex mesh not loaded properly from file " << fileName);
  }

  /// Create the sets of nodes and zones in the mesh
  IndexType numNodes = Repository::intsRegistry.getScalar("num_nodes");
  mesh->nodes = HexMesh::NodeSet(numNodes);

  IndexType numZones = Repository::intsRegistry.getScalar("num_zones");
  mesh->zones = HexMesh::ZoneSet(numZones);

  /// Create the nodal position field
  mesh->nodePosition = HexMesh::PositionsVec( &mesh->nodes );
  RealBuf::iterator ptIt
    = Repository::realsRegistry.getBuffer("node_positions").begin();
  for(IndexType idx = 0 ; idx < mesh->numNodes() ; ++idx)
  {
    mesh->nodePosition[idx] = Point(*ptIt++, *ptIt++, *ptIt++);
  }

  /// Create the topological incidence relation from zones to nodes
  IndexBuf& zn_indices
    = Repository::intsRegistry.getBuffer("zone_node_indices");
  mesh->zoneToNodeRelation
    = HexMesh::ZoneToNodeRelation(&mesh->zones, &mesh->nodes);
  mesh->zoneToNodeRelation.bindIndices(zn_indices.size(), &zn_indices);

  // Check that the relation is valid
  SLIC_ASSERT_MSG( mesh->zoneToNodeRelation.isValid(),
                   "Error creating (static) relation from zones to nodes!");
  SLIC_INFO("-- numNodesOfZones: " << zn_indices.size());

}

void generateNodeZoneRelation(HexMesh* mesh)
{
  // Create NodeToZone relation by inverting the ZoneToZone relation
  // TODO: This function to invert a relation should be moved into Slam

  typedef Repository::IntsRegistry::BufferType IndexBuf;
  typedef HexMesh::ZoneToNodeRelation::RelationSubset RelationSubset;

  /// Step 1: Compute the cardinalities of each node by looping through zone to
  // node relation
  IndexBuf& nzBegins = Repository::intsRegistry
                       .addBuffer("node_zone_begins", mesh->nodes.size() + 1 );
  for(IndexType zIdx = 0 ; zIdx < mesh->numZones() ; ++zIdx)
  {
    RelationSubset nSet = mesh->zoneToNodeRelation[zIdx];
    for(IndexType idx = 0 ; idx < nSet.size() ; ++idx)
    {
      ++nzBegins[ nSet[idx] ];
    }
  }

  /// Step 2: Compute begin offsets for each node based on cardinalities
  // Strategy: perform (inplace) exclusive prefix sum of cardinalities in
  // nzBegins
  IndexType prevVal = nzBegins[0];
  nzBegins[0] = 0;
  for(int i = 1 ; i <= mesh->numNodes() ; ++i)
  {
    IndexType nextVal = nzBegins[i];
    nzBegins[i] = nzBegins[i - 1] + prevVal;
    prevVal = nextVal;
  }


  /// Step 3: Invert the zone_node relation, use nzBegins[node_index] as offset
  // for next zone
  IndexBuf& zIndices = Repository::intsRegistry
                       .addBuffer("node_zone_indices",
                                  nzBegins[mesh->numNodes()] );
  for(IndexType zIdx = 0 ; zIdx < mesh->numZones() ; ++zIdx)
  {
    RelationSubset nSet = mesh->zoneToNodeRelation[zIdx];
    for(IndexType idx = 0 ; idx < nSet.size() ; ++idx)
    {
      const IndexType nIdx = nSet[idx];
      const IndexType offset = nzBegins[nIdx]++;
      zIndices[offset] = zIdx;
    }
  }

  /// Step 4: Fix begin offsets by shifting back by one index
  for(int i = mesh->numNodes() ; i > 0 ; --i)
  {
    nzBegins[i] = nzBegins[i - 1];
  }
  nzBegins[0] = 0;


  /// We can finally create the node to zone relation
  mesh->nodeToZoneRelation
    = HexMesh::NodeToZoneRelation( &mesh->nodes, &mesh->zones);
  mesh->nodeToZoneRelation.bindBeginOffsets(mesh->nodes.size(), &nzBegins);
  mesh->nodeToZoneRelation.bindIndices(zIndices.size(), &zIndices);

  SLIC_ASSERT_MSG( mesh->nodeToZoneRelation.isValid(true),
                   "Error creating (static) relation from nodes to zones!\n");

  SLIC_INFO("-- numZonesOfNode: " << zIndices.size());

}

void computeZoneBarycenters(HexMesh* mesh)
{
  typedef HexMesh::ZoneToNodeRelation::RelationSubset NodeSet;

  // Compute the zone positions as the the averages
  // of the positions of the nodes around each zone
  mesh->zonePosition = HexMesh::PositionsVec( &mesh->zones );

  // Outer loop over each zone in the mesh
  for(IndexType zIdx = 0 ; zIdx < mesh->numZones() ; ++zIdx )
  {
    Point zonePos;

    // Inner loop over each node of the zone
    const NodeSet& nodeSet = mesh->zoneToNodeRelation[zIdx];
    for(IndexType idx = 0 ; idx < nodeSet.size() ; ++idx)
    {
      zonePos += mesh->nodePosition[ nodeSet[idx] ];
    }
    zonePos /= nodeSet.size();

    mesh->zonePosition[ zIdx ] = zonePos;
  }
}

void createZoneRadiusField (HexMesh* mesh)
{
  // Compute a zone field based on the L2 norm of their position vectors
  mesh->zoneField = HexMesh::ZoneField ( &mesh->zones );

  for (IndexType zIdx = 0 ; zIdx < mesh->numZones() ; ++zIdx )
  {
    mesh->zoneField[zIdx] = mesh->zonePosition[zIdx].radius();
  }
}

DataType computeNodalErrors(HexMesh* mesh)
{
  // Compute the node average version, and the approximation error
  typedef HexMesh::NodeToZoneRelation::RelationSubset ZoneSet;


  mesh->nodeFieldAvg   = HexMesh::NodeField(&mesh->nodes);
  mesh->nodeFieldExact = HexMesh::NodeField(&mesh->nodes);
  double errSqSum = 0.0;

  // Outer loop over each node
  for (IndexType nIdx = 0 ; nIdx < mesh->numNodes() ; ++nIdx)
  {
    // What's the radius?
    const DataType nodalValExact = mesh->nodePosition[nIdx].radius();

    // Inner loop over each zone of the node to find average value
    const ZoneSet& zoneSet = mesh->nodeToZoneRelation[nIdx];
    DataType nodalAvg = DataType();
    for(IndexType idx = 0 ; idx < zoneSet.size() ; ++idx)
    {
      nodalAvg += mesh->zoneField[ zoneSet[idx] ];
    }
    nodalAvg /= zoneSet.size();

    // Update fields and compute error
    mesh->nodeFieldAvg[nIdx] = nodalAvg;
    mesh->nodeFieldExact[nIdx] = nodalValExact;

    const double err = nodalAvg - nodalValExact;
    errSqSum += err * err;
  }

  DataType err = std::sqrt( errSqSum / mesh->numNodes() );
  SLIC_INFO("-> The L2-ish error in the node average radius was " << err);

  return err;
}

}   // end namespace slamUnstructuredHex


int main(int argc, char** argv)
{
  using namespace slamUnstructuredHex;

  axom::slic::UnitTestLogger logger;

#ifndef USE_ONE
  int const NUM_RESOLUTIONS = 4;
#else
  int const NUM_RESOLUTIONS = 1;
#endif

  int fileResolutions[] = {1,2,4,8};
  DataType expectedResults[]
    = {0.10736689892, 0.037977237476, 0.013251067479, 0.0046357167735};

  std::string dataDir;
  if(argc > 1)
  {
    dataDir = std::string(argv[1]);
  }
  else
  {
    // Parse command line for data directory, with fallback
    const std::string DEFAULT_DATA_DIR = "../src/components/slam/data";
    dataDir = DEFAULT_DATA_DIR;
    SLIC_INFO("Using default data directory "
              << DEFAULT_DATA_DIR
              << "\n First parameter can be a custom directory.");
  }

  int numFailedTests = 0;
  for(int res = 0 ; res < NUM_RESOLUTIONS ; ++res)
  {
    std::string meshName = fmt::format("{}/ball_{}.vtk",
                                       dataDir, fileResolutions[res]);

    SLIC_INFO("Loading mesh file '"
              << meshName
              << "' and generating zone-> node relation");

    HexMesh hexMesh;
    readHexMesh( meshName, &hexMesh );
    SLIC_ASSERT(hexMesh.zoneToNodeRelation.isValid());
    //--------------------------------------------------------------

    // Now build the node to zone relation
    SLIC_INFO("Generating node->zone relation");
    generateNodeZoneRelation( &hexMesh );
    SLIC_ASSERT(hexMesh.nodeToZoneRelation.isValid());

    //--------------------------------------------------------------
    // Now that we have the mesh in memory, we can start to do things with it.
    SLIC_INFO("Computing zone barycenters using zone->node relation");
    computeZoneBarycenters(&hexMesh);

    SLIC_INFO("Generating a zone-centered radius field");
    createZoneRadiusField(&hexMesh);

    DataType errVal = computeNodalErrors(&hexMesh);

    // Some error checking based on precomputed values
    if(!axom::utilities::isNearlyEqual(errVal, expectedResults[res]) )
    {
      SLIC_WARNING(
        "Error differed from expected value -- "
        << fmt::format("Expected {}, but got {} (difference: {}",
                       expectedResults[res], errVal,
                       errVal - expectedResults[res]));

      ++numFailedTests;
    }

    SLIC_INFO("-- done.\n");
  }

  //--------------------------------------------------------------
  SLIC_INFO(fmt::format("-- {} tests out of {} passed",
                        NUM_RESOLUTIONS - numFailedTests, NUM_RESOLUTIONS));

  return (numFailedTests == 0) ? 0 : 1;
}

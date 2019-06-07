// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"  // for axom macros
// #include "axom/mir.hpp"  // for Mir classes & functions
#include "axom/slam.hpp"



// C/C++ includes
#include <cmath>          // for definition of M_PI, exp()
#include <vector>

// namespace aliases
//namespace mir     = axom::mir;
namespace numerics = axom::numerics;
namespace slam = axom::slam;

#define EPSILON_ZERO = 0.00000000001f;

//--------------------------------------------------------------------------------

/**
 * \brief Simple 2D Point class for example
 */
struct Point2
{
  Point2(double x = 0., double y = 0.) : m_x(x), m_y(y) {}

  Point2(const Point2& other) : m_x(other.m_x), m_y(other.m_y) {}

  Point2& operator=(const Point2& other)
  { m_x = other.m_x; m_y = other.m_y; return *this; }

  Point2& operator+=(const Point2& other)
  { m_x += other.m_x; m_y += other.m_y; return *this; }

  Point2& operator/=(double val)
  { m_x /= val; m_y += val; return *this; }

  double& operator[] (int i) { return (i==0) ? m_x : m_y; }
  const double& operator[] (int i) const { return (i==0) ? m_x : m_y; }

  friend std::ostream& operator<<(std::ostream& os, const Point2& pt)
  { return os << "{x:" << pt.m_x << ", y:" << pt.m_y <<"}"; }

  double m_x, m_y;
};

//--------------------------------------------------------------------------------

struct EdgeClipInfo
{
  int vertexOne;
  int vertexTwo;
  int colorOne;
  int colorTwo;
  float t;        // The percent distance from vertexOne to VertexTwo where the clip occurs.
};

//--------------------------------------------------------------------------------

struct MIRMesh
{
  /****************************************************************
   *                        TYPE ALIASES 
   ****************************************************************/
  // SET TYPE ALIASES
  using PosType = slam::DefaultPositionType;
  using ElemType = slam::DefaultElementType;

  using ArrayIndir = slam::policies::ArrayIndirection< PosType, ElemType >;

  using VertSet = slam::PositionSet< PosType, ElemType >;
  using ElemSet = slam::PositionSet< PosType, ElemType >;

  // RELATION TYPE ALIASES
  using VarCard = slam::policies::VariableCardinality< PosType, ArrayIndir >;

  // Note: This is the actual relation type, which takes in a bunch of policies and data entries to relate to each other.
  // Note: It is the relation of the elements to the vertices.
  using ElemToVertRelation = slam::StaticRelation< PosType, ElemType, VarCard, ArrayIndir, ElemSet, VertSet >;
  using VertToElemRelation = slam::StaticRelation< PosType, ElemType, VarCard, ArrayIndir, VertSet, ElemSet >;

  // MAP TYPE ALIASES
  using BaseSet = slam::Set< PosType, ElemType >;
  using ScalarMap = slam::Map< BaseSet, axom::float64 >;  // Note: Documentation has a typo in it.
  using PointMap = slam::Map< BaseSet, Point2 >;

  enum { GREEN = 0, BLUE = 1 };

  /****************************************************************
   *                       MESH FUNCTIONS
   ****************************************************************/
  void InitializeMesh()
  {
    // Create the mesh connectivity information
    // data for element-vertex boundary relation
    evInds = {
      0,4,5,1,     // elem 0, card 4, start 0
      1,5,6,2,     // elem 1, card 4, start 4
      2,6,7,3,     // elem 2, card 4, start 8
      4,8,9,5,     // elem 3, card 4, start 12
      5,9,10,6,    // elem 4, card 4, start 16
      6,10,11,7,   // elem 5, card 4, start 20
      8,12,13,9,   // elem 6, card 4, start 24
      9,13,14,10,  // elem 7, card 4, start 28
      10,14,15,11  // elem 8, card 4, start 32, end 36
    };

    evBegins = {
      0,4,8,12,16,20,24,28,32,36
    };

    // data for vertex-element coboundary relation
    veInds = {
      0,          // vert  0, card 1, start 0
      0,1,        // vert  1, card 2, start 1
      1,2,        // vert  2, card 2, start 3
      2,          // vert  3, card 1, start 5
      0,3,        // vert  4, card 2, start 6
      0,1,3,4,    // vert  5, card 4, start 8
      1,2,4,5,    // vert  6, card 4, start 12
      2,5,        // vert  7, card 2, start 16
      3,6,        // vert  8, card 2, start 18
      3,4,6,7,    // vert  9, card 4, start 20
      4,5,7,8,    // vert  10, card 4, start 24
      5,8,        // vert  11, card 2, start 28
      6,          // vert  12, card 1, start 30
      6,7,        // vert  13, card 2, start 31
      7,8,        // vert  14, card 2, start 33
      8,          // vert  15, card 1, start 35, end 36
    };
    veBegins = {
      0,1,3,5,6,8,12,16,18,20,24,28,30,31,33,35,36
    };

    // Initialize the mesh sets
    verts = VertSet(16);  // Construct a vertex set with 11 vertices
    elems = ElemSet(9);   // Construct an element set with 5 elements

    // Check validity of the sets
    SLIC_ASSERT_MSG( verts.isValid(), "Vertex set is not valid.");
    SLIC_ASSERT_MSG( elems.isValid(), "Element set is not valid.");

    printf("Mesh has been initialized.\n");
  }

  /// Constructs the mesh boundary and coboundary relations
  void constructMeshRelations()
  {

    // construct boundary relation from elements to vertices using variable cardinality
    {
      using RelationBuilder = ElemToVertRelation::RelationBuilder;
      bdry = RelationBuilder()
              .fromSet( &elems )
              .toSet( &verts )
              .begins( RelationBuilder::BeginsSetBuilder()
                        .size( elems.size() ) 
                        .data( evBegins.data() )  )
              .indices ( RelationBuilder::IndicesSetBuilder() 
                          .size( evInds.size() ) 
                          .data( evInds.data() ) );
    }


    {
      // _quadmesh_example_construct_cobdry_relation_start
      // construct coboundary relation from vertices to elements
      using RelationBuilder = VertToElemRelation::RelationBuilder;
      cobdry = RelationBuilder()
               .fromSet( &verts )
               .toSet( &elems )
               .begins( RelationBuilder::BeginsSetBuilder()
                        .size( verts.size() )
                        .data( veBegins.data() ) )
               .indices( RelationBuilder::IndicesSetBuilder()
                         .size( veInds.size() )
                         .data( veInds.data() ) );
      // _quadmesh_example_construct_cobdry_relation_end
    }


    SLIC_ASSERT_MSG( bdry.isValid(), "Boundary relation is not valid.");
    SLIC_ASSERT_MSG( cobdry.isValid(), "Coboundary relation is not valid.");

    SLIC_INFO("Elem-Vert relation has size " << bdry.totalSize());
    SLIC_INFO("Vert-Elem relation has size " << cobdry.totalSize());

    printf("Finished constructing mesh relations.\n");
  }

  /// Constructs the volume fraction maps on the vertices
  void constructMeshVolumeFractionMaps()
  {
    // TODO: Determine better way to store volume fraction data than 1 map per material?
    // Construct the scalar map on the elements for the green material
    greenVolumeFractionsElements = ScalarMap( &elems );

    // Set the green volume fraction for each element 
    greenVolumeFractionsElements[0] = 1.0;
    greenVolumeFractionsElements[1] = 1.0;
    greenVolumeFractionsElements[2] = 1.0;
    greenVolumeFractionsElements[3] = 1.0;
    greenVolumeFractionsElements[4] = 0.5;
    greenVolumeFractionsElements[5] = 0.2;
    greenVolumeFractionsElements[6] = 0.2;
    greenVolumeFractionsElements[7] = 0.0;
    greenVolumeFractionsElements[8] = 0.0;
    
    // Construct the scalar map on the elements for the blue material
    blueVolumeFractionsElements = ScalarMap( &elems );

    blueVolumeFractionsElements[0] = 0.0;
    blueVolumeFractionsElements[1] = 0.0;
    blueVolumeFractionsElements[2] = 0.0;
    blueVolumeFractionsElements[3] = 0.0;
    blueVolumeFractionsElements[4] = 0.5;
    blueVolumeFractionsElements[5] = 0.8;
    blueVolumeFractionsElements[6] = 0.8;
    blueVolumeFractionsElements[7] = 1.0;
    blueVolumeFractionsElements[8] = 1.0;

    printf("Finished constructing volume fractions for mesh elements.\n");
  }

    /// Constucts the positions map on the vertices
  void constructVertexPositionMap()
  {
    // construct the position map on the vertices
    vertexPositions = PointMap( &verts );

    // first vertex is at origin
    vertexPositions[0] = Point2( 0.0, 3.0 );
    vertexPositions[1] = Point2( 1.0, 3.0 );
    vertexPositions[2] = Point2( 2.0, 3.0 );
    vertexPositions[3] = Point2( 3.0, 3.0 );

    vertexPositions[4] = Point2( 0.0, 2.0 );
    vertexPositions[5] = Point2( 1.0, 2.0 );
    vertexPositions[6] = Point2( 2.0, 2.0 );
    vertexPositions[7] = Point2( 3.0, 2.0 );

    vertexPositions[8] = Point2( 0.0, 1.0 );
    vertexPositions[9] = Point2( 1.0, 1.0 );
    vertexPositions[10] = Point2( 2.0, 1.0 );
    vertexPositions[11] = Point2( 3.0, 1.0 );

    vertexPositions[12] = Point2( 0.0, 0.0 );
    vertexPositions[13] = Point2( 1.0, 0.0 );
    vertexPositions[14] = Point2( 2.0, 0.0 );
    vertexPositions[15] = Point2( 3.0, 0.0 );

    SLIC_ASSERT_MSG( vertexPositions.isValid(), "Position map is not valid.");

    SLIC_INFO("-- Vertex positions:");
    for(int vID=0 ; vID < verts.size() ; ++vID)
    {
      SLIC_INFO("Position of vert " << vID << " is " << vertexPositions[vID]);
    }
  }

  // Computes the average volume fraction at each vertex of the mesh
  void computeVolumeFractionAverages()
  {
    // Initialize the vertex volume fraction maps
    greenVolumeFractionsVertices = ScalarMap( &verts );
    blueVolumeFractionsVertices = ScalarMap( &verts );

    for (int vID = 0; vID < verts.size(); ++vID)
    {
      // Compute the per vertex volume fractions for the green material
      axom::float64 sum = 0;
      auto vertexElements = cobdry[vID];
      for (int i = 0; i < vertexElements.size(); ++i)
      {
        auto eID = vertexElements[i];
        sum += greenVolumeFractionsElements[eID];
      }
      greenVolumeFractionsVertices[vID] = sum / vertexElements.size();

      // Compute the per vertex volume fractions for the blue material
      sum = 0;
      for (int i = 0; i < vertexElements.size(); ++i)
      {
        auto eID = vertexElements[i];
        sum += blueVolumeFractionsElements[eID];
      }
      blueVolumeFractionsVertices[vID] = sum / vertexElements.size();
    }

    printf("Finished computing volume fraction averages.\n");
  }

  /// Prints out the map values for each element
  void printElementScalarMap(ScalarMap& elements, std::string prefix)
  {
    std::cout << prefix;
    for (int i = 0; i < elems.size(); ++i)
    {
      printf("Element %d: %f\n", i, elements[i]);
    }
  }

  /// Prints out the map values for each vertex
  void printVertexScalarMap(ScalarMap& vertices, std::string prefix)
  {
    std::cout << prefix;
    for (int i = 0; i < verts.size(); ++i)
    {
      printf("Vertex %d: %f\n", i, vertices[i]);
    }
  }

  /// Use bilinear interpolation to compute the points where the given element should be clipped.
  /// Returns true if clip should occur, otherwise false.
  /// If a clip should occur, the clipping points will be given in p1 and p2.
  void computeClippingPoints(const int eID)//, Point2* p1, Point2* p2)
  {
    // Find the vertices associated with each element
    auto elementVertices = bdry[eID];

    // Check if one of the two currently considered materials is not present at the vertices of the current element
    if (blueVolumeFractionsElements[eID] != 0.0 && greenVolumeFractionsElements[eID] != 0.0)
    {
      // Triangle Case
      if (elementVertices.size() == 3)
      {
        computeTriangleClippingPoints(eID);
      }
      // Quad Case
      if (elementVertices.size() == 4)
      {
        computeQuadClippingPoints(eID);
      }
    }
    else
    {
      printf("No cuts in element %d.\n", eID);
    }
  }

  /// Computes the points where the given quad element should be clipped
  void computeQuadClippingPoints(int eID)
  {
    auto elementVertices = bdry[eID];
    EdgeClipInfo eci[4];

    // Initialize the edge clipping info
    for (int i = 0; i < 4; ++i)
    {
      eci[i].vertexOne = elementVertices[ i ];
      eci[i].vertexTwo = elementVertices[ (i + 1) % 4 ];

      if (blueVolumeFractionsVertices[eci[i].vertexOne] > greenVolumeFractionsVertices[eci[i].vertexOne])
        eci[i].colorOne = BLUE;
      else
        eci[i].colorOne = GREEN;

      if (blueVolumeFractionsVertices[eci[i].vertexTwo] > greenVolumeFractionsVertices[eci[i].vertexTwo])
        eci[i].colorTwo = BLUE;
      else
        eci[i].colorTwo = GREEN;

      eci[i].t = -1.0;  // default t value means one material dominates the other and no clip should occur on this edge
    }

    // Analyze the edge clipping info and determine where to clip the cell
    for (int i = 0; i < 4; ++i)
    {
      if (eci[i].colorOne != eci[i].colorTwo)
      {
        axom::float64 numerator = blueVolumeFractionsVertices[eci[i].vertexOne] - greenVolumeFractionsVertices[eci[i].vertexOne];
        axom::float64 denominator = -greenVolumeFractionsVertices[eci[i].vertexOne]
                                    + greenVolumeFractionsVertices[eci[i].vertexTwo]
                                    + blueVolumeFractionsVertices[eci[i].vertexOne]
                                    - blueVolumeFractionsVertices[eci[i].vertexTwo];
        
        if (denominator != 0.0)
          eci[i].t = numerator / denominator;
      }
    }

    // Print out the cutting information
    bool cut = false;
    for (int i = 0; i < 4; ++i)
    {
      if (eci[i].t > 0.0)
      {
        printf("Cutting element %d between vertices %d and %d. t = %f\n", eID, eci[i].vertexOne, eci[i].vertexTwo, eci[i].t);
        cut = true;
      }
    }
    if (!cut)
      printf("No cuts in element %d.\n", eID);
      
  }

  /// Computes the points where the given triangle element should be clipped
  void computeTriangleClippingPoints(int eID)
  {
    printf("Triangle clipping case not yet implemented.\n");
  }

  /****************************************************************
   *                        VARIABLES
   ****************************************************************/
  public:
    // Mesh Set Definitions
    VertSet verts;  // the set of vertices in the mesh
    ElemSet elems;  // the set of elements in the mesh

    // Mesh Relation Definitions
    ElemToVertRelation bdry;      // Boundary relation from elements to vertices
    VertToElemRelation cobdry;   // Coboundary relation from vertices to elements

    // Mesh Map Definitions
    PointMap vertexPositions;   // vertex position
    ScalarMap greenVolumeFractionsElements;   // the volume fractions of the green material for each element (NOT each vertex)
    ScalarMap blueVolumeFractionsElements;    // the volume fractions of the blue material for each element (NOT each vertex)

    ScalarMap greenVolumeFractionsVertices;
    ScalarMap blueVolumeFractionsVertices;

    // support data for mesh connectivity
    std::vector<PosType> evInds;
    std::vector<PosType> evBegins;
    std::vector<PosType> veInds;
    std::vector<PosType> veBegins;
};


//--------------------------------------------------------------------------------


/*!
 * \file
 *
 * \brief Various code snippets/examples used for the Mir tutorial section.
 *
 * \note These examples are designed to illustrate specific Mir
 *  concepts and capabilities. Consult the Tutorial section of Mir's
 *  User Guide for more details.
 *
 */


/*!
 * \brief Tutorial main
 */
int main( int AXOM_NOT_USED(argc), char** AXOM_NOT_USED(argv) )
{
  printf("Beginning the MIR Tutorial program.\n");

  MIRMesh testMesh;
  testMesh.InitializeMesh();
  testMesh.constructMeshRelations();
  // testMesh.constructMeshPositionMap();
  testMesh.constructMeshVolumeFractionMaps();
  testMesh.computeVolumeFractionAverages();

  // Print volume fraction debugging information
  // testMesh.printElementScalarMap(testMesh.greenVolumeFractionsElements, "Green Volume Fractions per Element: \n");
  // testMesh.printElementScalarMap(testMesh.blueVolumeFractionsElements, "Blue Volume Fractions per Element: \n");

  // testMesh.printVertexScalarMap(testMesh.greenVolumeFractionsVertices, "Green Volume Fractions per Vertex: \n");
  // testMesh.printVertexScalarMap(testMesh.blueVolumeFractionsVertices, "Blue Volume Fractions per Vertex: \n");

  // Determine where to clip the cells
  for (int i = 0; i < 9; ++i)
    testMesh.computeClippingPoints(i);


  return 0;
}

//--------------------------------------------------------------------------------
#ifndef QUEST_TEST_UTILITIES_HPP_
#define QUEST_TEST_UTILITIES_HPP_


#include "slic/slic.hpp"

#include "quest/Point.hpp"
#include "quest/Triangle.hpp"

#include "mint/Mesh.hpp"
#include "mint/UnstructuredMesh.hpp"
#include "mint/Field.hpp"
#include "mint/FieldData.hpp"


/**
 * \file
 *
 * This file contains several utility functions for testing the quest component,
 * e.g. generating random doubles and Points, creating a simple mesh of an octahedron...
 * We may later decide to move some of these into the actual component if they are deemed useful.
 */

namespace quest {
namespace utilities {


/**
 * \brief Simple utility function to generate a random double in the range of beg to end
 * \param beg The lower value of the range (default 0.)
 * \param end The upper value of the range (default 1.)
 * \pre If the range is zero, we set it to 1
 * \return A double in the range [beg,end]
 */
double randomDouble(double beg = 0., double end = 1.)
{
    double range = end-beg;

    if(range == 0)
        range = 1.;

    return beg + (rand() / ( RAND_MAX / range ) ) ;
}

 /**
  * \brief Simple utility to generate a Point whose entries are random values in the range [beg, end]
  */
template<int DIM>
quest::Point<double,DIM> randomSpacePt(double beg, double end)
{
    quest::Point<double,DIM> pt;
    for(int i=0; i< DIM; ++i)
        pt[i] = randomDouble(beg,end);

    return pt;
}


/**
 * \brief Simple utility to find the centroid of two points
 */
template<int DIM>
quest::Point<double,DIM> getCentroid(const quest::Point<double,DIM>& pt0, const quest::Point<double,DIM>& pt1)
{
    return (pt0.array() + pt1.array()) /2.;
}

/**
 * \brief Simple utility to find the centroid of three points
 */
template<int DIM>
quest::Point<double,DIM> getCentroid(const quest::Point<double,DIM>& pt0, const quest::Point<double,DIM>& pt1, const quest::Point<double,DIM>& pt2)
{
    return (pt0.array() + pt1.array() + pt2.array()) /3.;
}

/**
 * \brief Utility function to generate a triangle mesh of an octahedron
 * Vertices of the octahedron are at +-i, +-j and +-k.
 * \note The caller must delete the mesh
 */
mint::Mesh*  make_octahedron_mesh()
{
    typedef int VertexIndex;
    typedef quest::Point<double, 3> SpacePt;
    typedef quest::Triangle<double, 3> SpaceTriangle;

    enum { POS_X, NEG_X, POS_Y, NEG_Y, POS_Z, NEG_Z };

    // The six vertices of the octahedron
    const int NUM_VERTS = 6;
    SpacePt verts[NUM_VERTS]
                  =  { SpacePt::make_point( 1., 0., 0.)
                     , SpacePt::make_point(-1., 0., 0.)
                     , SpacePt::make_point( 0,  1., 0.)
                     , SpacePt::make_point( 0, -1., 0.)
                     , SpacePt::make_point( 0,  0,  1.)
                     , SpacePt::make_point( 0,  0, -1.)  };

    // The eight triangles of the octahedron
    // Explicit representation of triangle-vertex incidence relation
    // Note: We are orienting the triangles with normals pointing outside
    const int NUM_TRIS = 8;
    const int VERTS_PER_TRI = 3;
    VertexIndex tvRelation[NUM_TRIS*VERTS_PER_TRI]
                  = { POS_Z, POS_X, POS_Y
                    , POS_Z, POS_Y, NEG_X
                    , POS_Z, NEG_X, NEG_Y
                    , POS_Z, NEG_Y, POS_X
                    , NEG_Z, POS_Y, POS_X
                    , NEG_Z, NEG_X, POS_Y
                    , NEG_Z, NEG_Y, NEG_X
                    , NEG_Z, POS_X, NEG_Y };

      // Note (KW 3/2016) -- We are not currently using this
      // Explicit representation of edge-vertex incidence relation
      // Note: we don't care about the orientation here
      //const int NUM_EDGES = 12;
      //const int VERTS_PER_EDGE = 3;
      //VertexIndex evRelation[NUM_EDGES*VERTS_PER_EDGE]
      //              = { POS_Z, POS_X  // Four edges incident in +Z
      //                , POS_Z, POS_Y
      //                , POS_Z, NEG_X
      //                , POS_Z, NEG_Y
      //                , NEG_Z, POS_X  // Four edges incident in -Z
      //                , NEG_Z, POS_Y
      //                , NEG_Z, NEG_X
      //                , NEG_Z, NEG_Y
      //                , POS_Y, POS_X  // Four edges not incident in Z
      //                , NEG_Y, POS_Y
      //                , POS_Y, NEG_X
      //                , NEG_Y, NEG_Y };

    // First, confirm that all triangles have normals that point away from the origin
    for(int i =0; i < NUM_TRIS; ++i)
    {
        int baseIndex = i*VERTS_PER_TRI;
        SpaceTriangle tri( verts[ tvRelation[ baseIndex + 0]]
                         , verts[ tvRelation[ baseIndex + 1]]
                         , verts[ tvRelation[ baseIndex + 2]] );

        SLIC_ASSERT( quest::ON_NEGATIVE_SIDE == quest::orientation( SpacePt(), tri) );
    }

    // Now create an unstructured triangle mesh from the two arrays
    typedef mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;
    TriangleMesh* triMesh = new TriangleMesh(3);

    // insert verts
    for(int i=0; i< NUM_VERTS; ++i)
        triMesh->insertNode(verts[i][0], verts[i][1], verts[i][2]);

    // insert triangles
    for(int i=0; i< NUM_TRIS; ++i)
        triMesh->insertCell( &tvRelation[i*VERTS_PER_TRI], MINT_TRIANGLE, 3);

    SLIC_ASSERT( NUM_VERTS == triMesh->getMeshNumberOfNodes() );
    SLIC_ASSERT( NUM_TRIS == triMesh->getMeshNumberOfCells() );

    return triMesh;
}


/**
 * \brief Utility function to write a VTK file from a mint Mesh instance
 */
void write_vtk( mint::Mesh* mesh, const std::string& fileName )
{
  SLIC_ASSERT( mesh != ATK_NULLPTR );

  std::ofstream ofs;
  ofs.open( fileName.c_str() );

  // STEP 0: Write VTK header
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << " Unstructured Mesh (";
  ofs << mesh->getBlockId() << ", " << mesh->getPartitionId() << ")\n";

  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";

  // STEP 1: Write mesh nodes
  const int num_nodes = mesh->getMeshNumberOfNodes();
  ofs << "POINTS " << num_nodes << " double\n";
  for ( int node=0; node < num_nodes; ++node ) {

    ofs << mesh->getMeshNodeCoordinate( node, 0 ) << " ";
    ofs << mesh->getMeshNodeCoordinate( node, 1 ) << " ";
    if ( mesh->getDimension()==3 ) {

        ofs << mesh->getMeshNodeCoordinate( node, 2 );

    } else {

        ofs << "0.0";
    }

    ofs << std::endl;

  } // END for all nodes

  // STEP 2: Write mesh cell connectivity
  // TODO: note currently this does not work with mixed cell types
  const int ncells   = mesh->getMeshNumberOfCells();
  const int maxnodes = mesh->getMeshNumberOfCellNodes( 0 );
  ofs << "CELLS " << ncells << " " << ncells*(maxnodes+1) << std::endl;

  std::vector< int > cell;
  for ( int cellIdx=0; cellIdx < ncells; ++cellIdx ) {

    const int nnodes = mesh->getMeshNumberOfCellNodes( cellIdx );
    cell.resize( nnodes );
    mesh->getMeshCell( cellIdx, &cell[0] );

    ofs << nnodes << " ";
    for ( int i=0; i < nnodes; ++i ) {
      ofs << cell[ i ] << " ";
    } // END for all nodes
    ofs << std::endl;

  } // END for all cells

  // STEP 3: Write cell types
  ofs << "CELL_TYPES " << ncells << std::endl;
  for ( int cellIdx=0; cellIdx < ncells; ++cellIdx ) {
    int ctype    = mesh->getMeshCellType( cellIdx );
    int vtk_type = mint::cell::vtk_types[ ctype ];
    ofs << vtk_type << std::endl;
  } // END for all cells

  // STEP 4: Write Cell Data
  ofs << "CELL_DATA " << ncells << std::endl;
  mint::FieldData* CD = mesh->getCellFieldData();
  for ( int f=0; f < CD->getNumberOfFields(); ++f ) {

      mint::Field* field = CD->getField( f );

      ofs << "SCALARS " << field->getName() << " ";
      if ( field->getType() == mint::DOUBLE_FIELD_TYPE ) {

          double* dataPtr = field->getDoublePtr();
          SLIC_ASSERT( dataPtr != ATK_NULLPTR );

          ofs << "double\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < ncells; ++i) {
              ofs << dataPtr[ i ] << std::endl;
          }

      } else {

          int* dataPtr = field->getIntPtr();
          SLIC_ASSERT( dataPtr != ATK_NULLPTR );

          ofs << "int\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < ncells; ++i ) {
             ofs << dataPtr[ i ] << std::endl;
          }
      }



  }

  // STEP 5: Write Point Data
  const int nnodes = mesh->getMeshNumberOfNodes();
  ofs << "POINT_DATA " << nnodes << std::endl;
  mint::FieldData* PD = mesh->getNodeFieldData();
  for ( int f=0; f < PD->getNumberOfFields(); ++f ) {

      mint::Field* field = PD->getField( f );

      ofs << "SCALARS " << field->getName() << " ";
      if ( field->getType() == mint::DOUBLE_FIELD_TYPE ) {

          double* dataPtr = field->getDoublePtr();
          ofs << "double\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < nnodes; ++i) {
              ofs << dataPtr[ i ] << std::endl;
          }

      } else {

          int* dataPtr = field->getIntPtr();
          ofs << "int\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < nnodes; ++i) {
              ofs << dataPtr[ i ] << std::endl;
          }

      }
  }

  ofs.close();
}






}   // end namespace utilities
}   // end namespace quest


#endif // QUEST_TEST_UTILITIES_HPP_

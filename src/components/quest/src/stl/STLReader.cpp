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

#include "STLReader.hpp"

// axom includes
#include "axom/Macros.hpp"
#include "axom/Types.hpp"

#include "axom_utils/Utilities.hpp"  // For isLittleEndian() and swapEndian()

#include "slic/slic.hpp"


// C/C++ includes
#include <cstddef>   // for NULL
#include <fstream>   // for ifstream


namespace
{
const std::size_t BINARY_HEADER_SIZE = 80;      // bytes
const std::size_t BINARY_TRI_SIZE = 50;         // bytes
}


//------------------------------------------------------------------------------
//      STLReader Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace quest
{

STLReader::STLReader() :
  m_fileName(""),
  m_num_nodes(0),
  m_num_faces(0)
{}

//------------------------------------------------------------------------------
STLReader::~STLReader()
{
  this->clear();
}

//------------------------------------------------------------------------------
void STLReader::clear()
{
  m_num_nodes = 0;
  m_num_faces = 0;
  m_nodes.clear();
}

//------------------------------------------------------------------------------
bool STLReader::isAsciiFormat() const
{
  // The binary format consists of
  //    a header of size BINARY_HEADER_SIZE==80 bytes
  //    followed by a four byte int encoding the number of triangles
  //    followed by the triangle data (BINARY_TRI_SIZE == 50 bytes per triangle)

  // Open the file
  std::ifstream ifs( m_fileName.c_str(), std::ios::in| std::ios::binary);
  SLIC_ASSERT_MSG(
    ifs.is_open(),
    "There was a problem reading the provided STL file " << m_fileName);

  // Find out the file size
  ifs.seekg(0, ifs.end);
  axom::common::int32 fileSize = static_cast<axom::common::int32>(ifs.tellg());

  const int totalHeaderSize =
    (BINARY_HEADER_SIZE + sizeof(axom::common::int32));
  if(fileSize < totalHeaderSize)
    return true;

  // Find the number of triangles (if the file were binary)
  int numTris = 0;
  ifs.seekg(BINARY_HEADER_SIZE, ifs.beg);
  ifs.read( (char*)&numTris, sizeof(axom::common::int32));

  if(!axom::utilities::isLittleEndian() )
  {
    numTris = axom::utilities::swapEndian(numTris);
  }

  // Check if the size matches our expectation
  int expectedBinarySize = totalHeaderSize + (numTris * BINARY_TRI_SIZE);
  return (fileSize != expectedBinarySize);
}

//------------------------------------------------------------------------------
void STLReader::readAsciiSTL()
{
  std::ifstream ifs( m_fileName.c_str());
  SLIC_ASSERT_MSG(
    ifs.is_open(),
    "There was a problem reading the provided STL file " << m_fileName);

  std::string junk;
  double x,y,z;

  // In an STL  file, we only care about the vertex positions
  // Vertices are strings of the form: "vertex v_x v_y v_z"
  while(true)
  {
    do
    {
      ifs >> junk;
    }
    while( ifs.good() && junk != "vertex");

    if(ifs.fail())
      break;

    ifs >> x >> y >> z;
    m_nodes.push_back(x);
    m_nodes.push_back(y);
    m_nodes.push_back(z);
  }

  // Set the number of nodes and faces
  m_num_nodes = m_nodes.size() / 3;
  m_num_faces = m_num_nodes / 3;
}

void STLReader::readBinarySTL()
{
  // Binary STL format consists of
  //    an 80 byte header (BINARY_HEADER_SIZE)
  //    followed by a 32 bit int encoding the number of faces
  //    followed by the triangles, each of which is 50 bytes (BINARY_TRI_SIZE)

  // A local union data structure for triangles in a binary STL
  union BinarySTLTri
  {
    axom::common::int8 raw[BINARY_TRI_SIZE];
    struct
    {
      float normal[3];
      float vert[9];
      axom::common::uint16 attr;
    };
  } tri;

  bool const isLittleEndian = axom::utilities::isLittleEndian();

  // Open binary file, skip the header
  std::ifstream ifs( m_fileName.c_str(), std::ios::in| std::ios::binary);
  ifs.seekg(BINARY_HEADER_SIZE);

  // read the num faces and reserve room for the vertex positions
  ifs.read( (char*)&m_num_faces, sizeof(axom::common::int32));

  if(!isLittleEndian )
  {
    m_num_faces = axom::utilities::swapEndian(m_num_faces);
  }

  m_num_nodes = m_num_faces * 3;
  m_nodes.reserve( m_num_nodes * 3);

  // Read the triangles. Cast to doubles and ignore normals and attributes
  for(axom::mint::IndexType i=0 ; i < m_num_faces ; ++i)
  {
    ifs.read( (char*)tri.raw, BINARY_TRI_SIZE);

    for(int j=0 ; j<9 ; ++j)
    {
      float coord = isLittleEndian
                    ? tri.vert[j]
                    : axom::utilities::swapEndian(tri.vert[j]);

      m_nodes.push_back( static_cast<double>( coord) );
    }
  }

}


//------------------------------------------------------------------------------
void STLReader::read()
{
  SLIC_ASSERT( m_fileName != "" );

  // Clear internal data, check the format and load the data
  this->clear();

  if(isAsciiFormat())
    readAsciiSTL();
  else
    readBinarySTL();
}

//------------------------------------------------------------------------------
void STLReader::getMesh(
  axom::mint::UnstructuredMesh< mint::SINGLE_SHAPE >* mesh )
{
  /* Sanity checks */
  SLIC_ASSERT( mesh != AXOM_NULLPTR );
  SLIC_ASSERT(
    static_cast<axom::mint::IndexType>(m_nodes.size()) == 3* m_num_nodes );

  // Load the vertices into the mesh
  for ( axom::mint::IndexType i=0 ; i < m_num_nodes ; ++i )
  {
    mesh->appendNode( m_nodes[i*3], m_nodes[i*3+1], m_nodes[i*3+2] );
  }

  // Load the triangles.  Note that the indices are implicitly defined.
  for ( axom::mint::IndexType i=0 ; i < m_num_faces ; ++i )
  {
    axom::mint::IndexType tv[3] = {3*i, 3*i+1, 3*i+2};
    mesh->appendCell( tv );
  }

}

} // end namespace quest
} // end namespace axom

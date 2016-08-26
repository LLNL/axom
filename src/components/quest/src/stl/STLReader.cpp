/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file STLReader.cpp
 *******************************************************************************
 */

#include "STLReader.hpp"

// ATK includes
#include "common/ATKMacros.hpp"
#include "common/CommonTypes.hpp"
#include "slic/slic.hpp"


// C/C++ includes
#include <cstddef>   // for NULL
#include <fstream>   // for ifstream



//------------------------------------------------------------------------------
//      STLReader Implementation
//------------------------------------------------------------------------------
namespace quest
{

STLReader::STLReader() :
        m_fileName(""),
        m_num_nodes(0),
        m_num_faces(0)
{

}

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
    // An STL file is in ASCII format if the first word is 'solid'

    std::ifstream ifs( m_fileName.c_str());
    SLIC_ASSERT_MSG(ifs.is_open()
                   , "There was a problem reading the provided STL file " << m_fileName);

    std::string first;
    ifs >> first;

    return first == "solid";
}

//------------------------------------------------------------------------------
void STLReader::readAsciiSTL()
{
    std::ifstream ifs( m_fileName.c_str());
    SLIC_ASSERT_MSG(ifs.is_open()
                   , "There was a problem reading the provided STL file " << m_fileName);

    std::string junk;
    double x,y,z;

    // In an STL  file, we only care about the vertex positions
    // Vertices are strings of the form: "vertex v_x v_y v_z"
    while(true)
    {
        do { ifs >> junk;} while( ifs.good() && junk != "vertex");

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
  const std::size_t  BINARY_HEADER_SIZE = 80; // bytes
  const std::size_t  BINARY_TRI_SIZE = 50;    // bytes

  // Binary STL format consists of
  //    an 80 byte header
  //    followed by a 32 bit int encoding the number of faces
  //    followed by the triangles, each of which is 50 bytes


  // A local union data structure for triangles in a binary STL
  union BinarySTLTri {
    asctoolkit::common::int8 raw[BINARY_TRI_SIZE];
    struct {
      float normal[3];
      float vert[9];
      asctoolkit::common::uint16 attr;
    };
  } tri;

  std::ifstream ifs( m_fileName.c_str(), std::ios::in| std::ios::binary);

  // skip the header
  ifs.seekg(BINARY_HEADER_SIZE);

  // read the num faces and reserve room for the vertex positions
  ifs.read( (char*)&m_num_faces, 4);

  m_num_nodes = m_num_faces * 3;
  m_nodes.reserve( m_num_nodes * 3);

  // Read the triangles. Cast to doubles and ignore normals and attributes
  for(int i=0; i < m_num_faces; ++i)
  {
    ifs.read( (char*)tri.raw, BINARY_TRI_SIZE);

    for(int j=0; j<9; ++j)
    {
      m_nodes.push_back( static_cast<double>( tri.vert[j]));
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
        mint::UnstructuredMesh< mint::LINEAR_TRIANGLE >* mesh )
{
  /* Sanity checks */
  SLIC_ASSERT( mesh != ATK_NULLPTR );
  SLIC_ASSERT( static_cast<int>(m_nodes.size()) == 3* m_num_nodes );

  // Load the vertices into the mesh
  for ( int i=0; i < m_num_nodes; ++i ) {
      mesh->insertNode( m_nodes[i*3], m_nodes[i*3+1], m_nodes[i*3+2] );
  }

  // Load the triangles.  Note that the indices are implicitly defined.
  for ( int i=0; i < m_num_faces; ++i ) {
      int tv[3] = {3*i, 3*i+1, 3*i+2};
      mesh->insertCell( tv,mint::LINEAR_TRIANGLE,3);
  }

}

} /* namespace quest */


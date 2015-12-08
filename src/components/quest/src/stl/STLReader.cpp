/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file STLReader.cpp
 *
 * \date Dec 8, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "STLReader.hpp"

// ATK includes
#include "common/CommonTypes.hpp"
#include "quest/stla_io.hpp"
#include "slic/slic.hpp"


// C/C++ includes
#include <cstddef>  // for NULL

namespace quest
{

STLReader::STLReader() :
        m_fileName(""),
        m_num_nodes(0),
        m_num_faces(0),
        m_nodes(ATK_NULLPTR),
        m_face_normals(ATK_NULLPTR),
        m_face_connectivity(ATK_NULLPTR)

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

  if ( m_nodes != ATK_NULLPTR ) {
      delete [] m_nodes;
        m_nodes = ATK_NULLPTR;
  }

  if ( m_face_normals != ATK_NULLPTR ) {
      delete [] m_face_normals;
      m_face_normals = ATK_NULLPTR;
  }

  if ( m_face_connectivity != ATK_NULLPTR ) {
      delete [] m_face_connectivity;
      m_face_connectivity = ATK_NULLPTR;
  }

}

//------------------------------------------------------------------------------
void STLReader::read()
{
  SLIC_ASSERT( m_fileName != "" );

  // STEP 0: clear internal data-structures
  this->clear();

  // STEP 1: Query STL file for sizes
  int numSolids = 0;
  int numNodes  = 0;
  int numFaces  = 0;
  int numText   = 0;
  stla_size( m_fileName, &numSolids, &numNodes, &numFaces, &numText );

  m_num_nodes = numNodes;
  m_num_faces = numFaces;

  // STEP 2: Allocate internal data-structures
  m_nodes             = new double[ 3*numNodes ];
  m_face_connectivity = new int[ 3*numFaces ];
  m_face_normals      = new double[ 3*numFaces ];

  // STEP 3: Read in geometry to internal data-structures
  stla_read( m_fileName, numNodes, numFaces,
          m_nodes, m_face_connectivity, m_face_normals );

}

//------------------------------------------------------------------------------
void STLReader::getMesh(
        meshtk::UnstructuredMesh< meshtk::LINEAR_TRIANGLE >* mesh )
{
  /* Sanity checks */
  SLIC_ASSERT( mesh != ATK_NULLPTR );
  SLIC_ASSERT( m_nodes != ATK_NULLPTR );
  SLIC_ASSERT( m_face_connectivity != ATK_NULLPTR );

  for ( int i=0; i < m_num_nodes; ++i ) {
      mesh->insertNode( m_nodes[i*3], m_nodes[i*3+1], m_nodes[i*3+2] );
  }

  for ( int i=0; i < m_num_faces; ++i ) {
      mesh->insertCell( &m_face_connectivity[ i*3],meshtk::LINEAR_TRIANGLE,3);
  }

}

} /* namespace quest */

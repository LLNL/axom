/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file PSTLReader.cpp
 *
 * \date Mar 10, 2016
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "PSTLReader.hpp"

#include "common/CommonTypes.hpp"
#include "quest/stla_io.hpp"
#include "slic/slic.hpp"

namespace quest {

//------------------------------------------------------------------------------
PSTLReader::PSTLReader( MPI_Comm comm ) : m_comm( comm )
{
  MPI_Comm_rank( m_comm, &m_my_rank );
}

//------------------------------------------------------------------------------
PSTLReader::~PSTLReader()
{
    // TODO Auto-generated destructor stub
}

//------------------------------------------------------------------------------
void PSTLReader::read()
{
  SLIC_ASSERT( m_fileName != "" );
  SLIC_ASSERT( m_comm != MPI_COMM_NULL );

  // STEP 0: clear internal data-structures
  this->clear();

  unsigned char* buffer = ATK_NULLPTR;
  int nbytes            = 0;

  // STEP 1: Query STL file for sizes
  if ( m_my_rank==0 ) {

     int numNodes  = 0;
     int numFaces  = 0;
     int numSolids = 0;
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

     // STEP 3: calculate total number of bytes to send
     nbytes = 2*sizeof(int)+ 3*numNodes*sizeof(double)+ 3*numFaces*sizeof(int);

     // STEP 4: serialize data in a buffer
     buffer = new unsigned char[ nbytes ];
     unsigned char* ptr = buffer;
     memcpy(ptr, &m_num_nodes, sizeof(int)  );
     ptr += sizeof(int);

     memcpy(ptr, &m_num_faces, sizeof(int) );
     ptr += sizeof(int);

     memcpy(ptr, m_nodes, 3*numNodes*sizeof(double) );
     ptr += 3*numNodes*sizeof(double);

     memcpy(ptr, m_face_connectivity, 3*numFaces*sizeof(int) );
     ptr += 3*numFaces*sizeof(int);

     // STEP 4: broadcast stuff to other ranks
     MPI_Bcast( &nbytes, 1, MPI_INT, 0, m_comm );
     MPI_Bcast( buffer, nbytes, MPI_UNSIGNED_CHAR, 0, m_comm );

  } else {

     // receive nbytes broadcast from rank o
     MPI_Bcast( &nbytes, 1, MPI_INT, 0, m_comm );

     // allocate & receive
     buffer = new unsigned char[ nbytes ];
     MPI_Bcast( buffer, nbytes, MPI_UNSIGNED_CHAR, 0, m_comm );

     // unpack
     unsigned char* ptr = buffer;
     memcpy( &m_num_nodes, ptr, sizeof(int) );
     ptr += sizeof(int);

     memcpy( &m_num_faces, ptr, sizeof(int) );
     ptr += sizeof(int);

     m_nodes             = new double[ 3*m_num_nodes ];
     memcpy( m_nodes, ptr, 3*m_num_nodes*sizeof(double) );
     ptr += 3*m_num_nodes*sizeof(double);

     m_face_connectivity = new int[ 3*m_num_faces ];
     memcpy( m_face_connectivity, ptr, 3*m_num_faces*sizeof(int) );
     ptr += 3*m_num_faces*sizeof(int);

  } // END else

  if ( buffer != ATK_NULLPTR ) {
     delete [] buffer;
     buffer = ATK_NULLPTR;
  }

  MPI_Barrier( MPI_COMM_WORLD );
}

} /* namespace quest */

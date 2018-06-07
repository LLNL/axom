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

#include "quest/PSTLReader.hpp"

namespace axom
{
namespace quest
{

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
int PSTLReader::read()
{
  SLIC_ASSERT( m_fileName != "" );
  SLIC_ASSERT( m_comm != MPI_COMM_NULL );

  // Clear internal data-structures
  this->clear();

  int rc = -1; // return code

  if( m_my_rank == 0)
  {

    // Rank 0 reads the mesh and broadcasts vertex positions to the others
    rc = STLReader::read();
    if ( rc == 0 )
    {
      MPI_Bcast( &m_num_nodes, 1, MPI_INT, 0, m_comm );
      MPI_Bcast( &m_nodes[0], m_num_nodes * 3, MPI_DOUBLE, 0, m_comm);

    }
    else
    {

      MPI_Bcast( &rc, 1, MPI_INT, 0, m_comm );

    }

  }
  else
  {

    // Other ranks receive the mesh vertices from rank 0
    MPI_Bcast( &m_num_nodes, 1, MPI_INT, 0, m_comm );
    if ( m_num_nodes !=-1 )
    {
      m_nodes.resize( m_num_nodes * 3);
      MPI_Bcast( &m_nodes[0], m_num_nodes * 3, MPI_DOUBLE, 0, m_comm);
      m_num_faces = m_num_nodes / 3;
      rc = 0;
    }

  } // END else

  MPI_Barrier( MPI_COMM_WORLD );

  return ( rc );
}

} // end namespace quest
} // end namespace axom

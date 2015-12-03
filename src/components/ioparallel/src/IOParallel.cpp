/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Implementation file for IOParallel class.
 *
 ******************************************************************************
 */

// Associated header file
#include "IOParallel.hpp"

// Other toolkit component headers
#include "common/CommonTypes.hpp"

// SiDRe project headers
#include "sidre/DataGroup.hpp"
#include "sidre/SidreTypes.hpp"


namespace asctoolkit
{
namespace ioparallel
{


/*
 *************************************************************************
 *
 * Datastore ctor creates root group.
 *
 *************************************************************************
 */
IOParallel::IOParallel(MPI_Comm comm,
  std::vector<sidre::DataGroup *>& groups,
  int num_files)
: m_comm_size(1),
  m_my_rank(0),
  m_baton(comm, num_files),
  m_datagroups(groups)
//  m_mpi_tag(378)
{
  MPI_Comm_size(comm, &m_comm_size);
  MPI_Comm_rank(comm, &m_my_rank);
}


/*
 *************************************************************************
 *
 * Datastore dtor destroys all contents.
 *
 *************************************************************************
 */
IOParallel::~IOParallel()
{
}


/*
 *************************************************************************
 *
 * Return non-cost pointer to buffer with given index or null ptr.
 *
 *************************************************************************
 */
void IOParallel::write(const std::string& file_string, int cycle, const std::string& protocol)
{
  int group_id = m_baton.waitForMyTurn();
  std::ostringstream oss;
  oss << file_string << "_" << group_id << "_" << cycle;
  std::string file_name = oss.str();
  for (std::vector<sidre::DataGroup *>::const_iterator itr = m_datagroups.begin();
       itr != m_datagroups.end(); ++itr) {  
    (*itr)->save(file_name, protocol);
  }
  (void)m_baton.finishMyTurn();
}

/*
 *************************************************************************
 *
 * Create new data buffer and assign unique id.
 *
 *************************************************************************
 */
void IOParallel::read(const std::string& file_string, int cycle, const std::string& protocol)
{
  int group_id = m_baton.waitForMyTurn();
  std::ostringstream oss;
  oss << file_string << "_" <<  group_id << "_" << cycle;
  std::string file_name = oss.str();
  for (std::vector<sidre::DataGroup *>::iterator itr = m_datagroups.begin();
       itr != m_datagroups.end(); ++itr) {
    (*itr)->load(file_name, protocol);
  }
  (void)m_baton.finishMyTurn();
}

#if 0
void IOParallel::waitForMyTurn()
{
   if (m_rank_before_me != -1) {
      MPI_Status mpi_stat;
      int baton;
      int mpi_err = MPI_Recv(&baton, 1, MPI_INT, m_rank_before_me,
         m_mpi_tag, m_mpi_com, &mpi_stat);
      if (mpi_err == MPI_SUCCESS) {


      }
   }
}
#endif




} /* end namespace sidre */
} /* end namespace asctoolkit */

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
 * \brief   Implementation file for IOManager class.
 *
 ******************************************************************************
 */

// Associated header file
#include "IOManager.hpp"

// Other toolkit component headers
#include "common/CommonTypes.hpp"

// SiDRe project headers
#include "sidre/DataGroup.hpp"
#include "sidre/SidreTypes.hpp"

#include "conduit_mpi.hpp"

namespace asctoolkit
{
namespace spio
{


/*
 *************************************************************************
 *
 * Datastore ctor creates root group.
 *
 *************************************************************************
 */
IOManager::IOManager(MPI_Comm comm,
  sidre::DataGroup ** groups,
  int num_datagroups,
  int num_files)
: m_comm_size(1),
  m_my_rank(0),
  m_baton(comm, num_files),
  m_datagroups(groups),
  m_num_datagroups(num_datagroups)
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
IOManager::~IOManager()
{
}


/*
 *************************************************************************
 *
 * Return non-cost pointer to buffer with given index or null ptr.
 *
 *************************************************************************
 */
void IOManager::write(const std::string& file_string, int cycle, const std::string& protocol)
{
  int group_id = m_baton.waitForMyTurn();
  std::ostringstream oss;
  oss << file_string << "_" << group_id << "_" << cycle;
  std::string file_name = oss.str();
//  conduit::mpi::send();
  for (int i = 0; i < m_num_datagroups; ++i) {
    m_datagroups[i]->save(file_name, protocol);
  }
//  if (m_num_datagroups == 1) {
//    if (m_baton.groupSize() == 1) {
//      m_datagroups[0]->save(file_name, protocol);
//    } else if (m_baton.isLastInGroup()) {
//    } else {
//      conduit::Node n;
//      m_datagroups[0]->copyToNode(n);
//    }
//  }   
  (void)m_baton.finishMyTurn();
}

/*
 *************************************************************************
 *
 * Create new data buffer and assign unique id.
 *
 *************************************************************************
 */
void IOManager::read(const std::string& file_string, int cycle, const std::string& protocol)
{
  int group_id = m_baton.waitForMyTurn();
  std::ostringstream oss;
  oss << file_string << "_" <<  group_id << "_" << cycle;
  std::string file_name = oss.str();
  for (int i = 0; i < m_num_datagroups; ++i) {
    m_datagroups[i]->load(file_name, protocol);
  }
  (void)m_baton.finishMyTurn();
}




} /* end namespace spio */
} /* end namespace asctoolkit */

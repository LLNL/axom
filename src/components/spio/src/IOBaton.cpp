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
 * \brief   Implementation file for IOBaton class.
 *
 ******************************************************************************
 */

// Associated header file
#include "IOBaton.hpp"

// Other toolkit component headers
#include "common/CommonTypes.hpp"

// SiDRe project headers
#include "sidre/DataGroup.hpp"
#include "sidre/SidreTypes.hpp"


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
IOBaton::IOBaton(MPI_Comm comm,
  int num_files)
: m_mpi_comm(comm),
  m_comm_size(1),
  m_my_rank(0),
  m_rank_before_me(-1),
  m_rank_after_me(-1),
  m_mpi_tag(378)
{
  MPI_Comm_size(comm, &m_comm_size);
  MPI_Comm_rank(comm, &m_my_rank);
  m_num_groups = num_files;
  m_num_larger_groups = m_comm_size % num_files;
  m_group_size = m_comm_size / m_num_groups; // ?
  m_first_regular_group_rank = (m_group_size + 1) * m_num_larger_groups;
  if (m_my_rank < m_first_regular_group_rank) {
    m_group_id = m_my_rank / (m_group_size + 1);
    m_rank_within_group = m_my_rank % (m_group_size + 1);
    if (m_rank_within_group < m_group_size) {
      m_rank_after_me = m_my_rank + 1;
    } 
  } else {
    m_group_id = m_num_larger_groups + (m_my_rank - m_first_regular_group_rank) / m_group_size;
    m_rank_within_group = (m_my_rank - m_first_regular_group_rank) % m_group_size;
    if (m_rank_within_group < m_group_size - 1) {
      m_rank_after_me = m_my_rank + 1;
    } 
  }
  if (m_rank_within_group > 0) {
    m_rank_before_me = m_my_rank - 1;
  }
}


/*
 *************************************************************************
 *
 * Datastore dtor destroys all contents.
 *
 *************************************************************************
 */
IOBaton::~IOBaton()
{
}


int IOBaton::wait()
{
  int return_val = -1;
  if (m_rank_before_me != -1) {
    MPI_Status mpi_stat;
    int baton;
    int mpi_err = MPI_Recv(&baton, 1, MPI_INT, m_rank_before_me,
      m_mpi_tag, m_mpi_comm, &mpi_stat);
    if (mpi_err == MPI_SUCCESS) {
      return_val = m_group_id;
    }
  } else {
    return_val = m_group_id;
  }
  return return_val;
}

int IOBaton::pass()
{
  int return_val = 0;
  if (m_rank_after_me != -1) {
    int baton; 
    int mpi_err = MPI_Ssend(&baton, 1, MPI_INT, m_rank_after_me,
      m_mpi_tag, m_mpi_comm);
    if (mpi_err != MPI_SUCCESS) {
      return_val = -1;
    }
  }
  return return_val;
}




} /* end namespace spio */
} /* end namespace asctoolkit */

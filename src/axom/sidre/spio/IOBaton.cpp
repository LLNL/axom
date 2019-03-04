/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

// Associated header file
#include "IOBaton.hpp"

const int axom::sidre::IOBaton::s_invalid_rank_id = -1;

namespace axom
{
namespace sidre
{


/*
 *************************************************************************
 *
 * IOBaton constructor
 *
 *************************************************************************
 */
IOBaton::IOBaton(MPI_Comm comm,
                 int num_files,
                 int num_trees)
  : m_mpi_comm(comm),
  m_comm_size(1),
  m_my_rank(0),
  m_rank_before_me(s_invalid_rank_id),
  m_rank_after_me(s_invalid_rank_id),
  m_mpi_tag(MPI_ANY_TAG)
{
  SLIC_ASSERT(num_files > 0);
  SLIC_ASSERT(num_trees > 0);

  MPI_Comm_size(comm, &m_comm_size);
  MPI_Comm_rank(comm, &m_my_rank);
  m_num_files = num_files;
  m_num_trees = num_trees;
  int active_comm_size = m_comm_size;
  if (m_comm_size > m_num_trees)
  {
    active_comm_size = m_num_trees;
  }
  m_num_larger_groups = active_comm_size % num_files;
  if (m_my_rank < active_comm_size)
  {
    m_group_size = active_comm_size / m_num_files; // ?
  }
  else
  {
    m_group_size = 1;
  }
  m_first_regular_group_rank = (m_group_size + 1) * m_num_larger_groups;
  if (m_my_rank < m_first_regular_group_rank)
  {
    m_group_id = m_my_rank / (m_group_size + 1);
    m_rank_within_group = m_my_rank % (m_group_size + 1);
    if (m_rank_within_group < m_group_size)
    {
      m_rank_after_me = m_my_rank + 1;
    }
  }
  else if (m_my_rank < active_comm_size)
  {
    m_group_id = m_num_larger_groups +
                 (m_my_rank - m_first_regular_group_rank) / m_group_size;
    m_rank_within_group = (m_my_rank - m_first_regular_group_rank) %
                          m_group_size;
    if (m_rank_within_group < m_group_size - 1)
    {
      m_rank_after_me = m_my_rank + 1;
    }
  }
  else
  {
    m_group_id = m_my_rank;
    m_rank_within_group = 0;
  }
  if (m_rank_within_group > 0)
  {
    m_rank_before_me = m_my_rank - 1;
  }
}


/*
 *************************************************************************
 *
 * IOBaton dtor destroys all contents.
 *
 *************************************************************************
 */
IOBaton::~IOBaton()
{}


int IOBaton::wait()
{
  int return_val = -1;
  if (m_rank_before_me != s_invalid_rank_id)
  {
    MPI_Status mpi_stat;
    int baton;
    int mpi_err = MPI_Recv(&baton, 1, MPI_INT, m_rank_before_me,
                           m_mpi_tag, m_mpi_comm, &mpi_stat);
    if (mpi_err == MPI_SUCCESS)
    {
      return_val = m_group_id;
    }
  }
  else
  {
    return_val = m_group_id;
  }
  return return_val;
}

int IOBaton::pass()
{
  int return_val = 0;
  if (m_rank_after_me != s_invalid_rank_id)
  {
    int baton;
    int mpi_err = MPI_Ssend(&baton, 1, MPI_INT, m_rank_after_me,
                            0, m_mpi_comm);
    if (mpi_err != MPI_SUCCESS)
    {
      return_val = -1;
    }
  }
  return return_val;
}




} /* end namespace sidre */
} /* end namespace axom */

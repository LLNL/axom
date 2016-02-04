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
 * \brief   Header file containing definition of IOBaton class.
 *
 ******************************************************************************
 */

#ifndef IOBATON_HPP_
#define IOBATON_HPP_

// Standard C++ headers
#include <vector>
#include <stack>

#include "mpi.h"

// Other CS Toolkit headers
#include "common/CommonTypes.hpp"
#include "sidre/DataGroup.hpp"
#include "slic/slic.hpp"


namespace asctoolkit
{
namespace spio
{

// using directives to make Conduit usage easier and less visible
//using conduit::Node;

/*!
 * \class IOBaton
 *
 * \brief IOBaton
 *
 * It dumps and reads
 */
class IOBaton
{
public:

  /*!
   * \brief Default ctor initializes IOBaton
   */
  IOBaton(MPI_Comm com,
          int num_files);

  /*!
   * \brief Dtor destroys
   */
  ~IOBaton();

  int waitForMyTurn();
  int finishMyTurn(); 

  int groupSize()
  {
    return m_my_rank < m_first_regular_group_rank ? m_group_size + 1 : m_group_size; 
  }

  bool isFirstInGroup()
  {
     return (m_rank_within_group == 0); 
  }

  bool isLastInGroup()
  {
     return (m_rank_after_me == -1); 
  }

private:
  /*!
   *  Unimplemented ctors and copy-assignment operators.
   */
#ifdef USE_CXX11
  IOBaton( const IOBaton& ) = delete;
  IOBaton( IOBaton&& ) = delete;

  IOBaton& operator=( const IOBaton& ) = delete;
  IOBaton& operator=( IOBaton&& ) = delete;
#else
  IOBaton( const IOBaton& );
  IOBaton& operator=( const IOBaton& );
#endif

  MPI_Comm m_mpi_comm;

  int m_comm_size;  // num procs in the MPI communicator
  int m_my_rank;    // rank of this proc
  int m_num_groups; // number of groups (files)
  int m_num_larger_groups;  // some group have one extra
  int m_group_size; // regular group size (m_comm_size / m_num_groups) w/o remainder
  int m_group_id;
  int m_first_regular_group_rank;
  int m_rank_within_group;
  int m_rank_before_me;
  int m_rank_after_me;

  int m_mpi_tag;
};


} /* end namespace spio */
} /* end namespace asctoolkit */

#endif /* IOBATON_HPP_ */

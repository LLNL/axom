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
 * \brief   Header file containing definition of IOParallel class.
 *
 ******************************************************************************
 */

#ifndef IOPARALLEL_HPP_
#define IOPARALLEL_HPP_

// Standard C++ headers
#include <vector>
#include <stack>

#include "mpi.h"

// Other CS Toolkit headers
#include "common/CommonTypes.hpp"
#include "sidre/DataGroup.hpp"
#include "spio/IOBaton.hpp"
#include "slic/slic.hpp"


namespace asctoolkit
{
namespace spio
{

// using directives to make Conduit usage easier and less visible
//using conduit::Node;

/*!
 * \class IOParallel
 *
 * \brief IOParallel
 *
 * It dumps and reads
 */
class IOParallel
{
public:

  /*!
   * \brief Default ctor initializes IOParallel
   */
  IOParallel(MPI_Comm com,
             std::vector<sidre::DataGroup *>& groups,
             int num_files);

  /*!
   * \brief Dtor destroys
   */
  ~IOParallel();

  void write(const std::string& file_string, int cycle, const std::string& protocol);
  void read(const std::string& file_string, int cycle, const std::string& protocol);

private:
  /*!
   *  Unimplemented ctors and copy-assignment operators.
   */
#ifdef USE_CXX11
  IOParallel( const IOParallel& ) = delete;
  IOParallel( IOParallel&& ) = delete;

  IOParallel& operator=( const IOParallel& ) = delete;
  IOParallel& operator=( IOParallel&& ) = delete;
#else
  IOParallel( const IOParallel& );
  IOParallel& operator=( const IOParallel& );
#endif

  int m_comm_size;  // num procs in the MPI communicator
  int m_my_rank;    // rank of this proc

  IOBaton m_baton;

  std::vector<sidre::DataGroup *> m_datagroups;

};


} /* end namespace spio */
} /* end namespace asctoolkit */

#endif /* IOPARALLEL_HPP_ */

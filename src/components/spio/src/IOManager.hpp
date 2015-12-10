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
 * \brief   Header file containing definition of IOManager class.
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


/*!
 * \class IOManager
 *
 * \brief IOManager
 *
 * It dumps and reads
 */
class IOManager
{
public:

  /*!
   * \brief Default ctor initializes IOManager
   */
  IOManager(MPI_Comm com,
             std::vector<sidre::DataGroup *>& groups,
             int num_files);

  /*!
   * \brief Dtor destroys
   */
  ~IOManager();

  void write(const std::string& file_string, int cycle, const std::string& protocol);
  void read(const std::string& file_string, int cycle, const std::string& protocol);

private:
  /*!
   *  Unimplemented ctors and copy-assignment operators.
   */
#ifdef USE_CXX11
  IOManager( const IOManager& ) = delete;
  IOManager( IOManager&& ) = delete;

  IOManager& operator=( const IOManager& ) = delete;
  IOManager& operator=( IOManager&& ) = delete;
#else
  IOManager( const IOManager& );
  IOManager& operator=( const IOManager& );
#endif

  int m_comm_size;  // num procs in the MPI communicator
  int m_my_rank;    // rank of this proc

  IOBaton m_baton;

  std::vector<sidre::DataGroup *> m_datagroups;

};


} /* end namespace spio */
} /* end namespace asctoolkit */

#endif /* IOPARALLEL_HPP_ */

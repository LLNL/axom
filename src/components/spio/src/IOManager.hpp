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

#include "mpi.h"

#include "hdf5.h"

// Other CS Toolkit headers
#include "common/ATKMacros.hpp"
#include "common/CommonTypes.hpp"
#include "sidre/DataGroup.hpp"
#include "spio/IOBaton.hpp"


namespace asctoolkit
{
namespace spio
{


/*!
 * \class IOManager
 *
 * \brief IOManager manages and organizes the I/O operations.
 *
 * This class handles the bookkeeping and organizing tasks that must be done
 * before calling DataGroup's I/O methods.  It uses IOBaton to control the
 * parallel I/O operations, such that one rank at a time interacts with any
 * particular output file.
 */
class IOManager
{
public:

  /*!
   * \brief Constructor
   *
   * \param com               MPI communicator
   * \param group            Pointer to a DataGroup holding data for I/O
   * \param num_files         Number of files for I/O
   */
  IOManager(MPI_Comm com,
            sidre::DataGroup * group,
            int num_files);

  /*!
   * \brief Destructor
   */
  ~IOManager();

  /*!
   * \brief write
   *
   * \param file_string   base name for output file
   * \param cycle         cycle counter
   * \param protocol      identifies I/O protocol (format, e
   */
  void write(const std::string& file_string,
             int cycle,
             const std::string& protocol);

  /*!
   * \brief read from input files
   *
   * \param file_string   base name of input files
   * \param cycle         cycle counter
   * \param protocol      identifies I/O protocol
   */
  void read(const std::string& file_string,
            int cycle,
            const std::string& protocol);

  /*!
   * \brief read from a root file
   *
   * \param root_file     root file containing input data
   */
  void read(const std::string& root_file);

private:

  DISABLE_COPY_AND_ASSIGNMENT( IOManager );

  void createRootFile(const std::string& root_name,
                      const std::string& file_base,
                      int cycle);

  std::string getHDF5FileName(hid_t root_file_id, int rankgroup_id);

  int m_comm_size;  // num procs in the MPI communicator
  int m_my_rank;    // rank of this proc

  IOBaton m_baton;

  sidre::DataGroup * m_datagroup;
  int m_num_files;

  MPI_Comm m_mpi_comm;
};


} /* end namespace spio */
} /* end namespace asctoolkit */

#endif /* IOPARALLEL_HPP_ */

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
   * \param groups            Array of pointers to DataGroups.
   * \param num_datagroups    Size of the groups array
   * \param num_files         Number of files for I/O
   */
  IOManager(MPI_Comm com,
            sidre::DataGroup ** groups,
            int num_datagroups,
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
   * \brief read
   *
   * \param file_string   base name of input file
   * \param cycle         cycle counter
   * \param protocol      identifies I/O protocol
   */
  void read(const std::string& file_string,
            int cycle,
            const std::string& protocol);

  void read(const std::string& root_file, const std::string& protocol);

private:

  DISABLE_COPY_AND_ASSIGNMENT( IOManager );

  void createRootFile(const std::string& root_name,
                      const std::string& file_base,
                      int cycle);

  std::string getHDF5FileName(hid_t root_file_id, int rankgroup_id, int datagroup_id);

  int m_comm_size;  // num procs in the MPI communicator
  int m_my_rank;    // rank of this proc

  IOBaton m_baton;

  sidre::DataGroup ** m_datagroups;
  int m_num_datagroups;
  int m_num_files;

  MPI_Comm m_mpi_comm;
};


} /* end namespace spio */
} /* end namespace asctoolkit */

#endif /* IOPARALLEL_HPP_ */

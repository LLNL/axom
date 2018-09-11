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

/*!
 ******************************************************************************
 *
 * \file IOManager.hpp
 *
 * \brief   Header file containing definition of IOManager class.
 *
 ******************************************************************************
 */

#ifndef IOPARALLEL_HPP_
#define IOPARALLEL_HPP_

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/sidre/core/Group.hpp"

// Spio headers
#include "axom/sidre/spio/IOBaton.hpp"

#include "mpi.h"

namespace axom
{
namespace sidre
{


/*!
 * \class IOManager
 *
 * \brief IOManager manages and organizes the I/O operations.
 *
 * This class handles the bookkeeping and organizing tasks that must be done
 * before calling Group's I/O methods.  It uses IOBaton to control the
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
   * \param use_scr           Use SCR library for scalable I/O management.
   *                          If true, the calling code must have already
   *                          called SCR_Init() after MPI_Init().
   */
  IOManager(MPI_Comm com, bool use_scr = false);

  /*!
   * \brief Destructor
   */
  ~IOManager();

  /*!
   * \brief write a Group to output files
   *
   * The Group, including all of its child groups and views, is written
   * to files according to the given protocol.
   *
   * valid protocols:
   *
   *    sidre_hdf5
   *    sidre_conduit_json
   *    sidre_json
   *
   *    conduit_hdf5
   *    conduit_bin
   *    conduit_json
   *    json
   *
   * \note The sidre_hdf5 and conduit_hdf5 protocols are only available
   * when Axom is configured with hdf5.
   *
   * \param group         Group to write to output
   * \param num_files     number of output data files
   * \param file_string   base name for output files
   * \param protocol      identifies I/O protocol
   */
  void write(sidre::Group* group,
             int num_files,
             const std::string& file_string,
             const std::string& protocol);

  /*!
   * \brief write additional group to existing root file
   *
   * Should be called after write().  The native layout of the group will
   * be added to the root file.
   *
   * This may be called more than once to write multiple groups to the file.
   *
   * This is not an MPI collective call.  It writes one group from one rank to
   * one root file.
   *
   * This currently only works if the root file was created for protocol
   * sidre_hdf5.
   *
   * \param group         Group to add to root file
   * \param file_name     name of existing root file
   */
  void writeGroupToRootFile(sidre::Group* group,
                            const std::string& file_name);

  /*!
   * \brief write additional group to a path inside an existing root file
   *
   * This should be called on an existing root file that already contains a
   * Group, most likely from a previous call to writeGroupToRootFile.
   *
   * group_path is a path that specifies a location somewhere in the native
   * layout of a Group already in the file.  The group added in this
   * method will be stored as a child group of the group specified by the
   * path.
   *
   * This is not an MPI collective call.  It writes one group from one rank to
   * one root file.
   *
   * This currently only works if the root file was created for protocol
   * sidre_hdf5.
   *
   * \param group         Group to add to root file
   * \param file_name     name of existing root file
   * \param group_path    path to a location within the root file
   */
  void writeGroupToRootFileAtPath(sidre::Group* group,
                                  const std::string& file_name,
                                  const std::string& group_path);

  /*!
   * \brief write additional group to a path inside an existing root file
   *
   * This should be called on an existing root file that already contains a
   * Group, most likely from a previous call to writeGroupToRootFile.
   *
   * group_path is a path that specifies a location somewhere in the native
   * layout of a Group already in the file.  The view added in this
   * method will be stored as a child view of the group specified by the
   * path.
   *
   * This is not an MPI collective call.  It writes one view from one rank to
   * one root file.
   *
   * This currently only works if the root file was created for protocol
   * sidre_hdf5.
   *
   * \param view          View to add to root file
   * \param file_name     name of existing root file
   * \param group_path    path to a location within the root file
   */
  void writeViewToRootFileAtPath(sidre::View* view,
                                 const std::string& file_name,
                                 const std::string& group_path);

  void writeBlueprintIndexToRootFile(DataStore* datastore,
                                     const std::string& domain_path,
                                     const std::string& file_name,
                                     const std::string& index_path);

  /*!
   * \brief read from input files
   *
   * \param group         Group to fill with input data
   * \param file_string   base name of input files
   * \param protocol      identifies I/O protocol
   * \param preserve_contents   Preserves group's existing contents if true
   */
  void read(sidre::Group* group,
            const std::string& file_string,
            const std::string& protocol,
            bool preserve_contents = false);

  /*!
   * \brief read from a root file
   *
   * \param group      Group to fill with input data
   * \param root_file  root file containing input data
   * \param preserve_contents   Preserves group's existing contents if true
   * \param use_scr    Use SCR to find and read the files.  This should be
   *                   set to true only if the files were written with SCR.
   */
  void read(sidre::Group* group,
            const std::string& root_file,
            bool preserve_contents = false,
            bool use_scr = false);

  /**
   * \brief Finds conduit relay protocol corresponding to a sidre protocol
   *
   * \param sidre_protocol String representing the sidre protocol
   * \return The conduit relay protocol corresponding to \a sidre_protocol
   * Options are: "hdf5", "json" and "conduit_json"
   * \see Group::save() for a list of valid sidre protocols
   */
  static std::string correspondingRelayProtocol(
    const std::string& sidre_protocol);

  /*!
   * \brief load external data into a group
   *
   * This currently only works if the root file was created for protocol
   * sidre_hdf5.
   *
   * \param group         Group to fill with external data from input
   * \param root_file     root file containing input data
   */
  void loadExternalData(sidre::Group* group,
                        const std::string& root_file);

  /*!
   * \brief gets the number of files in the dataset from the specified root file
   */
  int getNumFilesFromRoot(const std::string& root_file);

private:

  DISABLE_COPY_AND_ASSIGNMENT( IOManager );

  void createRootFile(const std::string& file_base,
                      int num_files,
                      const std::string& protocol);

  std::string getProtocol(const std::string& root_name);

  /*!
   * Collective operation to get the file pattern from the root file.
   * The string is read on rank 0 and broadcast to the other ranks.
   * \note Works for all sidre protocols.
   */
  std::string getFilePatternFromRoot(const std::string& root_name,
                                     const std::string& protocol);

#ifdef AXOM_USE_HDF5
  std::string getHDF5FilePattern(const std::string& root_name);

  void readSidreHDF5(sidre::Group* group, const std::string& root_file,
                     bool preserve_contents = false);
#endif /* AXOM_USE_HDF5 */

  std::string getFileNameForRank(const std::string& file_pattern,
                                 const std::string& root_name,
                                 int rankgroup_id);

#ifdef AXOM_USE_SCR
  void readWithSCR(sidre::Group* group,
                   const std::string& root_file,
                   bool preserve_contents = false);
#endif

  int m_comm_size;  // num procs in the MPI communicator
  int m_my_rank;    // rank of this proc

  IOBaton* m_baton;

  MPI_Comm m_mpi_comm;

  bool m_use_scr;
  std::string m_scr_checkpoint_dir;

};


} /* end namespace sidre */
} /* end namespace axom */

#endif /* IOPARALLEL_HPP_ */

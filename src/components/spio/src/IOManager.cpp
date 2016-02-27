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

#include "hdf5.h"

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
 * Create manager
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
 * Destroy all contents
 *
 *************************************************************************
 */
IOManager::~IOManager()
{
}


/*
 *************************************************************************
 *
 * Write to file.
 *
 *************************************************************************
 */
void IOManager::write(const std::string& file_string, int cycle, const std::string& protocol)
{
  int group_id = m_baton.wait();
  std::ostringstream namestream;
  namestream << file_string << "_" << group_id << "_" << cycle;
  std::string file_name = namestream.str();
  if (protocol == "conduit_hdf5") {
    for (int i = 0; i < m_num_datagroups; ++i) {
      std::ostringstream savestream;
      savestream << file_name << ".group" << i << ".hdf5";
      std::string hdf5_name = savestream.str();

      hid_t h5_file_id;
      if (m_baton.isFirstInGroup()) {
        h5_file_id = H5Fcreate(hdf5_name.c_str(),
                               H5F_ACC_TRUNC,
                               H5P_DEFAULT,
                               H5P_DEFAULT);   
      } else {
        h5_file_id = H5Fopen(hdf5_name.c_str(),
                             H5F_ACC_RDWR,
                             H5P_DEFAULT); 
      }

      savestream << ":datagroup" << i << "_" << m_my_rank << "/";
      m_datagroups[i]->save(savestream.str(), protocol, h5_file_id);

      H5Fclose(h5_file_id);
    }
  } else if (protocol == "conduit") {
    for (int i = 0; i < m_num_datagroups; ++i) {
      std::ostringstream savestream;
      savestream << file_name << ".group" << i;
      std::string obase = savestream.str();
      m_datagroups[i]->save(obase, protocol);
    } 
  }
  (void)m_baton.pass();
}

/*
 *************************************************************************
 *
 * Read from file
 *
 *************************************************************************
 */
void IOManager::read(const std::string& file_string, int cycle, const std::string& protocol)
{
  int group_id = m_baton.wait();
  std::ostringstream namestream;
  namestream << file_string << "_" <<  group_id << "_" << cycle;
  std::string file_name = namestream.str();
  if (protocol == "conduit_hdf5") {
    for (int i = 0; i < m_num_datagroups; ++i) {
      std::ostringstream loadstream;
      loadstream << file_name << ".group" << i << ".hdf5";
      std::string hdf5_name = loadstream.str();
      hid_t h5_file_id = H5Fopen(hdf5_name.c_str(),
                                 H5F_ACC_RDONLY,
                                 H5P_DEFAULT);

      loadstream << ":datagroup" << i << "_" << m_my_rank << "/";

      m_datagroups[i]->load(loadstream.str(), protocol, h5_file_id);

      H5Fclose(h5_file_id);
    }
  } else {
    for (int i = 0; i < m_num_datagroups; ++i) {
      std::ostringstream loadstream;
      loadstream << file_name << ".group" << i;
      std::string obase = loadstream.str();
      m_datagroups[i]->load(obase, protocol);
    }
  }

  (void)m_baton.pass();
}




} /* end namespace spio */
} /* end namespace asctoolkit */

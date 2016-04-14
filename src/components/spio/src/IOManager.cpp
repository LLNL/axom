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
#include "sidre/DataStore.hpp"
#include "sidre/SidreTypes.hpp"

//This does not appear to be needed.  TODO - Ask Noah if there are future plans to
//use it.
//#include "relay_mpi.hpp"

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

      hid_t h5_file_id, h5_group_id;
      herr_t status;
(void)status;//be quiet compiler, i will use this later when checks are added.
      if (m_baton.isFirstInGroup()) {
        h5_file_id = H5Fcreate(hdf5_name.c_str(),
                               H5F_ACC_TRUNC,
                               H5P_DEFAULT,
                               H5P_DEFAULT);   
        // TODO - Ask Noah about adding a SLIC check here, make sure file create succeeded.
      } else {
        h5_file_id = H5Fopen(hdf5_name.c_str(),
                             H5F_ACC_RDWR,
                             H5P_DEFAULT); 
        // TODO - Ask Noah about adding a SLIC check here, make sure file create succeeded.
      }

      std::ostringstream group_stream;
      group_stream << "datagroup" << i << "_" << m_my_rank;
      std::string group_name = group_stream.str();
      h5_group_id = H5Gcreate(h5_file_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      // TODO - Ask Noah about adding a SLIC check here, make sure group create succeeded.
      m_datagroups[i]->getDataStore()->save(h5_group_id, m_datagroups[i] );
      status = H5Gclose(h5_group_id);
      // TODO - Ask Noah about adding a SLIC check here, make sure group close succeeded.
      status = H5Fclose(h5_file_id);
      // TODO - Ask Noah about adding a SLIC check here, make sure file close succeeded.
    }
  } else if (protocol == "conduit") {
    for (int i = 0; i < m_num_datagroups; ++i) {
      std::ostringstream savestream;
      savestream << file_name << ".group" << i;
      std::string obase = savestream.str();
      m_datagroups[i]->getDataStore()->save(obase, protocol, m_datagroups[i]);
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


      // TODO Add HDF5 call to change hdf5 internal directory to loadstream name.
      std::ostringstream groupstream;
      groupstream << "datagroup" << i << "_" << m_my_rank;
      std::string group_name = groupstream.str();
      hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
      m_datagroups[i]->getDataStore()->load(h5_group_id, m_datagroups[i]);

      H5Fclose(h5_file_id);
    }
  } else {
    for (int i = 0; i < m_num_datagroups; ++i) {
      std::ostringstream loadstream;
      loadstream << file_name << ".group" << i;
      std::string obase = loadstream.str();
      m_datagroups[i]->getDataStore()->load(obase, protocol, m_datagroups[i]);
    }
  }

  (void)m_baton.pass();
}




} /* end namespace spio */
} /* end namespace asctoolkit */

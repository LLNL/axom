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
  int num_files)
: m_comm_size(1),
  m_my_rank(0),
  m_baton(comm, num_files),
  m_num_files(num_files),
  m_mpi_comm(comm)
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
void IOManager::write(sidre::DataGroup * datagroup, const std::string& file_string, int cycle, const std::string& protocol)
{

  std::ostringstream rootstream;
  rootstream << file_string << cycle << ".root";
  std::string root_name = rootstream.str();

  if (m_my_rank == 0 && protocol == "conduit_hdf5") {
     createRootFile(root_name, file_string, cycle);
  }
  MPI_Barrier(m_mpi_comm);
  int group_id = m_baton.wait();

  std::ostringstream namestream;
  namestream << file_string << "_" << group_id << "_" << cycle;
  std::string file_name = namestream.str();

  if (protocol == "conduit_hdf5") {

    hid_t root_file_id = H5Fopen(root_name.c_str(),
                                 H5F_ACC_RDWR,
                                 H5P_DEFAULT);

    SLIC_ASSERT(root_file_id >= 0);

    herr_t status;

    std::string hdf5_name = getHDF5FileName(root_file_id, group_id);

    hid_t h5_file_id, h5_group_id;
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
    SLIC_ASSERT(h5_file_id >= 0);

    std::ostringstream group_stream;
    group_stream << "datagroup_" << m_my_rank;
    std::string group_name = group_stream.str();
    h5_group_id = H5Gcreate(h5_file_id,
                            group_name.c_str(),
                            H5P_DEFAULT,
                            H5P_DEFAULT,
                            H5P_DEFAULT);
    SLIC_ASSERT(h5_group_id >= 0);

    datagroup->getDataStore()->save(h5_group_id, datagroup);

    status = H5Gclose(h5_group_id);
    SLIC_ASSERT(status >= 0);
    status = H5Fclose(h5_file_id);
    SLIC_ASSERT(status >= 0);

    status = H5Fclose(root_file_id);
    SLIC_ASSERT(status >= 0);


  } else if (protocol == "conduit") {
    std::ostringstream savestream;
    savestream << file_name << ".group";
    std::string obase = savestream.str();
    datagroup->getDataStore()->save(obase, protocol, datagroup);
  }
  (void)m_baton.pass();
}

/*
 *************************************************************************
 *
 * Read from files
 *
 *************************************************************************
 */
void IOManager::read(
  sidre::DataGroup * datagroup,
  const std::string& file_string,
  int cycle,
  const std::string& protocol)
{
  int group_id = m_baton.wait();
  std::ostringstream namestream;
  namestream << file_string << "_" <<  group_id << "_" << cycle;
  std::string file_name = namestream.str();
  if (protocol == "conduit_hdf5") {

    std::ostringstream rootstream;
    rootstream << file_string << cycle << ".root";
    std::string root_name = rootstream.str();

    hid_t root_file_id = H5Fopen(root_name.c_str(),
                                 H5F_ACC_RDWR,
                                 H5P_DEFAULT);
    SLIC_ASSERT(root_file_id >= 0);

    herr_t errv;

    std::string hdf5_name = getHDF5FileName(root_file_id, group_id);

    hid_t h5_file_id = H5Fopen(hdf5_name.c_str(),
                               H5F_ACC_RDONLY,
                               H5P_DEFAULT);
    SLIC_ASSERT(h5_file_id >= 0);

    // TODO Add HDF5 call to change hdf5 internal directory to loadstream name.
    std::ostringstream groupstream;
    groupstream << "datagroup_" << m_my_rank;
    std::string group_name = groupstream.str();
    hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
    SLIC_ASSERT(h5_file_id >= 0);
    datagroup->getDataStore()->load(h5_group_id, datagroup);

    errv = H5Fclose(h5_file_id);
    SLIC_ASSERT(errv >= 0);

    errv = H5Fclose(root_file_id);
    SLIC_ASSERT(errv >= 0);

  } else {
    std::ostringstream loadstream;
    loadstream << file_name << ".group";
    std::string obase = loadstream.str();
    datagroup->getDataStore()->load(obase, protocol, datagroup);
  }

  (void)m_baton.pass();
}

/*
 *************************************************************************
 *
 * Read based on HDF5 root file.
 *
 *************************************************************************
 */
void IOManager::read(sidre::DataGroup * datagroup, const std::string& root_file)
{
  int group_id = m_baton.wait();

  hid_t root_file_id = H5Fopen(root_file.c_str(),
                               H5F_ACC_RDWR,
                               H5P_DEFAULT);

  SLIC_ASSERT(root_file_id >= 0);

  herr_t errv;

  std::string hdf5_name = getHDF5FileName(root_file_id, group_id);

  hid_t h5_file_id = H5Fopen(hdf5_name.c_str(),
                             H5F_ACC_RDONLY,
                             H5P_DEFAULT);
  SLIC_ASSERT(h5_file_id >= 0);

  std::ostringstream groupstream;
  groupstream << "datagroup_" << m_my_rank;
  std::string group_name = groupstream.str();
  hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
  SLIC_ASSERT(h5_group_id >= 0);
  datagroup->getDataStore()->load(h5_group_id, datagroup);

  errv = H5Fclose(h5_file_id);
  SLIC_ASSERT(errv >= 0);

  errv = H5Fclose(root_file_id);
  SLIC_ASSERT(errv >= 0);

  (void)m_baton.pass();
}

/*
 *************************************************************************
 *
 * Create a root file. 
 *
 *************************************************************************
 */
void IOManager::createRootFile(const std::string& root_name,
                               const std::string& file_base,
                               int cycle)
{
  hid_t root_file_id;

  root_file_id = H5Fcreate(root_name.c_str(),
                           H5F_ACC_TRUNC,
                           H5P_DEFAULT,
                           H5P_DEFAULT);
  SLIC_ASSERT(root_file_id >= 0);

  hsize_t dim[] = { 1 };
  hid_t int_space = H5Screate_simple(1, dim, 0);
  SLIC_ASSERT(int_space >= 0);

  hid_t filesset = H5Dcreate(root_file_id, "num_files", H5T_NATIVE_INT,
    int_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  SLIC_ASSERT(filesset >= 0);

  herr_t errv = H5Dwrite(filesset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, &m_num_files);
  SLIC_ASSERT(errv >= 0);

  hid_t ranksset = H5Dcreate(root_file_id, "num_ranks", H5T_NATIVE_INT,
    int_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  SLIC_ASSERT(ranksset >= 0);

  errv = H5Dwrite(ranksset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, &m_comm_size);
  SLIC_ASSERT(errv >= 0);

  hid_t cycleset = H5Dcreate(root_file_id, "cycle", H5T_NATIVE_INT,
    int_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  SLIC_ASSERT(cycleset >= 0);

  errv = H5Dwrite(cycleset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, &cycle);
  SLIC_ASSERT(errv >= 0);

  hid_t files_group_id = H5Gcreate2(root_file_id,
                                    "files",
                                    H5P_DEFAULT,
                                    H5P_DEFAULT,
                                    H5P_DEFAULT);
  SLIC_ASSERT(files_group_id >= 0);

  for (int i = 0; i < m_num_files; ++i) {
    std::ostringstream basestream;
    basestream << file_base << "_" << i << "_" << cycle;
    std::string base_name = basestream.str();

    std::ostringstream savestream;
    savestream << "file_" << i;
    std::string file_label = savestream.str();

    hid_t file_id = H5Gcreate2(files_group_id,
                               file_label.c_str(),
                               H5P_DEFAULT,
                               H5P_DEFAULT,
                               H5P_DEFAULT);
    SLIC_ASSERT(file_id >= 0);

    std::ostringstream groupstream;
    groupstream << "group";
    std::string group_name = groupstream.str();

    std::ostringstream h5namestream;
    h5namestream << base_name << ".hdf5";
    std::string hdf5_name = h5namestream.str();

    hid_t atype = H5Tcopy(H5T_C_S1);
    SLIC_ASSERT(atype >= 0);

    errv = H5Tset_size(atype, hdf5_name.size()+1);
    SLIC_ASSERT(errv >= 0);

    errv = H5Tset_strpad(atype, H5T_STR_NULLTERM);
    SLIC_ASSERT(errv >= 0);

    hid_t space = H5Screate_simple(1, dim, 0);
    SLIC_ASSERT(space >= 0);

    hid_t dataset = H5Dcreate(file_id, group_name.c_str(),
      atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    SLIC_ASSERT(dataset >= 0);

    errv = H5Dwrite(dataset, atype, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, hdf5_name.c_str());
    SLIC_ASSERT(errv >= 0);
  }
  errv = H5Fflush(root_file_id, H5F_SCOPE_LOCAL);
  SLIC_ASSERT(errv >= 0);
  errv = H5Fclose(root_file_id);
  SLIC_ASSERT(errv >= 0);

}

/*
 *************************************************************************
 *
 * Get file name for a file holding DataGroup data.
 *
 *************************************************************************
 */
std::string IOManager::getHDF5FileName(
  hid_t root_file_id,
  int rankgroup_id)
{
  std::ostringstream pathstream;
  pathstream << "/files/file_" << rankgroup_id << "/group";
  std::string path_name = pathstream.str();

  hid_t h5_name_id = H5Dopen(root_file_id, path_name.c_str(), H5P_DEFAULT);
  SLIC_ASSERT(h5_name_id >= 0);

  hid_t dtype = H5Dget_type(h5_name_id);
  SLIC_ASSERT(dtype >= 0);
  size_t dsize = H5Tget_size(dtype);

  char* h5_name_buf = new char[dsize];
  SLIC_ASSERT(h5_name_buf);
  herr_t errv = H5Dread(h5_name_id,
                        dtype,
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        h5_name_buf);
  SLIC_ASSERT(errv >= 0);

  std::string hdf5_name(h5_name_buf);
  delete[] h5_name_buf;

  return hdf5_name;
}



} /* end namespace spio */
} /* end namespace asctoolkit */

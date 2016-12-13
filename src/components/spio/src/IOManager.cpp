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
#include "common/FileUtilities.hpp"

// SiDRe project headers
#include "sidre/DataGroup.hpp"
#include "sidre/DataStore.hpp"
#include "sidre/SidreTypes.hpp"

#include "conduit_relay.hpp"
#include "conduit_relay_hdf5.hpp"
#include "fmt/fmt.hpp"

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
IOManager::IOManager(MPI_Comm comm)
: m_comm_size(1),
  m_my_rank(0),
  m_baton(ATK_NULLPTR),
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
  if (m_baton) {
    delete m_baton;
  }
}


/*
 *************************************************************************
 *
 * Write to file.
 *
 *************************************************************************
 */
void IOManager::write(sidre::DataGroup * datagroup, int num_files, const std::string& file_string, const std::string& protocol)
{
  if (m_baton) {
    if (m_baton->getNumFiles() != num_files) {
      delete m_baton;
      m_baton = ATK_NULLPTR;
    }
  }
 
  if (!m_baton) {
    m_baton = new IOBaton(m_mpi_comm, num_files);
  }

  std::string root_name = file_string + ".root";

  if (m_my_rank == 0 && protocol == "sidre_hdf5") {
     createRootFile(root_name, file_string, num_files);
  }
  MPI_Barrier(m_mpi_comm);

  if (protocol == "sidre_hdf5") {
    std::string file_pattern = getHDF5FilePattern(root_name);

    int group_id = m_baton->wait();

    std::string hdf5_name = getHDF5FileName(file_pattern, root_name, group_id);

    hid_t h5_file_id, h5_group_id;
    if (m_baton->isFirstInGroup()) {
      std::string dir_name; 
      utilities::filesystem::getDirName(dir_name, hdf5_name);
      if (!dir_name.empty()) {
        utilities::filesystem::makeDirsForPath(dir_name);
      }
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

    std::string group_name = fmt::sprintf("datagroup_%07d", m_my_rank);
    h5_group_id = H5Gcreate(h5_file_id,
                            group_name.c_str(),
                            H5P_DEFAULT,
                            H5P_DEFAULT,
                            H5P_DEFAULT);
    SLIC_ASSERT(h5_group_id >= 0);

    datagroup->save(h5_group_id);

    herr_t status = H5Gclose(h5_group_id);
    SLIC_ASSERT(status >= 0);
    status = H5Fflush(h5_file_id, H5F_SCOPE_LOCAL);
    SLIC_ASSERT(status >= 0);
    status = H5Fclose(h5_file_id);
    SLIC_ASSERT(status >= 0);

  } else {
    int group_id = m_baton->wait();
    std::string file_name = fmt::sprintf("%s_%07d", file_string, group_id);

    std::string obase = file_name + "." + protocol;
    datagroup->save(obase, protocol);
  }
  (void)m_baton->pass();
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
  const std::string& protocol)
{
  if (protocol == "sidre_hdf5") {
    std::string root_name = file_string + ".root";

    std::string file_pattern = getHDF5FilePattern(root_name);

    int group_id = m_baton->wait();

    herr_t errv;

    std::string hdf5_name = getHDF5FileName(file_pattern, root_name, group_id);

    hid_t h5_file_id = H5Fopen(hdf5_name.c_str(),
                               H5F_ACC_RDONLY,
                               H5P_DEFAULT);
    SLIC_ASSERT(h5_file_id >= 0);

    // TODO Add HDF5 call to change hdf5 internal directory to loadstream name.
    std::string group_name = fmt::sprintf("datagroup_%07d", m_my_rank);
    hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
    SLIC_ASSERT(h5_file_id >= 0);
    datagroup->load(h5_group_id);

    errv = H5Gclose(h5_group_id);
    SLIC_ASSERT(errv >= 0);

    errv = H5Fclose(h5_file_id);
    SLIC_ASSERT(errv >= 0);

  } else {
    int group_id = m_baton->wait();
    std::string file_name = fmt::sprintf("%s_%07d", file_string, group_id);

    std::string obase = file_name + ".group";
    datagroup->load(obase, protocol);
  }

  (void)m_baton->pass();
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
  int num_files = getNumFilesFromRoot(root_file);
  SLIC_ASSERT(num_files > 0);

  if (m_baton) {
    if (m_baton->getNumFiles() != num_files) {
      delete m_baton;
      m_baton = ATK_NULLPTR;
    }
  }
    
  if (!m_baton) {
    m_baton = new IOBaton(m_mpi_comm, num_files);
  }

  std::string file_pattern = getHDF5FilePattern(root_file);

  int group_id = m_baton->wait();

  herr_t errv;

  std::string hdf5_name = getHDF5FileName(file_pattern, root_file, group_id);

  hid_t h5_file_id = H5Fopen(hdf5_name.c_str(),
                             H5F_ACC_RDONLY,
                             H5P_DEFAULT);
  SLIC_ASSERT(h5_file_id >= 0);

  std::string group_name = fmt::sprintf("datagroup_%07d", m_my_rank);
  hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
  SLIC_ASSERT(h5_group_id >= 0);

  datagroup->load(h5_group_id);

  errv = H5Gclose(h5_group_id);
  SLIC_ASSERT(errv >= 0);

  errv = H5Fclose(h5_file_id);
  SLIC_ASSERT(errv >= 0);

  (void)m_baton->pass();
}


void IOManager::loadExternalData(sidre::DataGroup * datagroup, const std::string& root_file)
{
  int num_files = getNumFilesFromRoot(root_file);
  SLIC_ASSERT(num_files > 0);

  if (m_baton) {
    if (m_baton->getNumFiles() != num_files) {
      delete m_baton;
      m_baton = ATK_NULLPTR;
    }
  }

  if (!m_baton) {
    m_baton = new IOBaton(m_mpi_comm, num_files);
  }

  std::string file_pattern = getHDF5FilePattern(root_file);

  int group_id = m_baton->wait();

  herr_t errv;

  std::string hdf5_name = getHDF5FileName(file_pattern, root_file, group_id);

  hid_t h5_file_id = H5Fopen(hdf5_name.c_str(),
                             H5F_ACC_RDONLY,
                             H5P_DEFAULT);
  SLIC_ASSERT(h5_file_id >= 0);

  std::string group_name = fmt::sprintf("datagroup_%07d", m_my_rank);
  hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
  SLIC_ASSERT(h5_group_id >= 0);

  datagroup->loadExternalData(h5_group_id);

  errv = H5Gclose(h5_group_id);
  SLIC_ASSERT(errv >= 0);

  errv = H5Fclose(h5_file_id);
  SLIC_ASSERT(errv >= 0);

  (void)m_baton->pass();
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
                               int num_files)
{
  conduit::Node n;

  n["number_of_files"] = num_files;
  std::string local_file_base;
  std::string next;
  std::string slash = "/";
  conduit::utils::rsplit_string(file_base, slash, local_file_base, next);
  n["file_pattern"] = local_file_base + slash + local_file_base + "_" + "%07d.hdf5";
  n["number_of_trees"] = m_comm_size;
  
  n["tree_pattern"] = "datagroup_%07d";
  n["protocol/name"] = "sidre_hdf5";
  n["protocol/version"] = "0.0";
  
  conduit::relay::io::save(n,root_name,"hdf5");
}

/*
 *************************************************************************
 *
 * Get namescheme pattern for a file holding DataGroup data.
 *
 *************************************************************************
 */
std::string IOManager::getHDF5FilePattern(
  const std::string& root_name)
{
  std::string file_pattern;
  int buf_size = 0;
  if (m_my_rank == 0) {
    conduit::Node n;
    conduit::relay::io::load(root_name + ":file_pattern", "hdf5", n);

    file_pattern = n.as_string();

    buf_size = file_pattern.size() + 1;
  }

  MPI_Bcast(&buf_size, 1, MPI_INT, 0, m_mpi_comm);

  char name_buf[buf_size];
  if (m_my_rank == 0) {
    strcpy(name_buf, file_pattern.c_str());
  }

  MPI_Bcast(name_buf, buf_size, MPI_CHAR, 0, m_mpi_comm);

  if (m_my_rank != 0) {
    file_pattern = std::string(name_buf);
  }

  return file_pattern;
}

/*
 *************************************************************************
 *
 * Get the file name based on a pattern and a group id
 *
 *************************************************************************
 */


std::string IOManager::getHDF5FileName(
  const std::string& file_pattern,
  const std::string& root_name,
  int rankgroup_id)
{
  std::string hdf5_name = fmt::sprintf(file_pattern.c_str(), rankgroup_id);

  //If the root file was given as a path, find the directory and add it to
  //hdf5_name.
  std::string curr;
  std::string root_dir;
  std::string slash = "/";
  conduit::utils::rsplit_string(root_name, slash, curr, root_dir);

  if (!root_dir.empty()) {
    hdf5_name = root_dir + slash + hdf5_name;
  }

  return hdf5_name;
}

/*
 *************************************************************************
 *
 * Query the rootfile for the number of input files.
 *
 *************************************************************************
 */
int IOManager::getNumFilesFromRoot(const std::string& root_file)
{
  /*
   * Read num_files from rootfile on rank 0.
   */
  int read_num_files = 0;
  if (m_my_rank == 0) {
    conduit::Node n;
    conduit::relay::io::load(root_file + ":number_of_files","hdf5",n);
    read_num_files = n.to_int();
    SLIC_ASSERT(read_num_files > 0);
  }

  /*
   * Reduction sets num_files on all ranks.
   */
  int num_files;
  MPI_Allreduce(&read_num_files, &num_files, 1, MPI_INT, MPI_SUM, m_mpi_comm);
  SLIC_ASSERT(num_files > 0);

  return num_files;
}

/*
 *************************************************************************
 *
 * Write a group to an existing root file
 *
 *************************************************************************
 */

void IOManager::writeGroupToRootFile(sidre::DataGroup * group,
                                     const std::string& file_name)
{
  hid_t root_file_id = H5Fopen(file_name.c_str(),
                               H5F_ACC_RDWR,
                               H5P_DEFAULT);

  SLIC_ASSERT(root_file_id >= 0); 

  hid_t group_id = H5Gcreate2(root_file_id,
                              group->getName().c_str(),
                              H5P_DEFAULT,
                              H5P_DEFAULT,
                              H5P_DEFAULT);
  SLIC_ASSERT(group_id >= 0);

  conduit::Node data_holder;
  group->createNativeLayout(data_holder);

  conduit::relay::io::hdf5_write(data_holder, group_id);

  herr_t errv = H5Gclose(group_id);
  SLIC_ASSERT(errv >= 0);

  errv = H5Fflush(root_file_id, H5F_SCOPE_LOCAL);
  SLIC_ASSERT(errv >= 0);

  errv =  H5Fclose(root_file_id);
  SLIC_ASSERT(errv >= 0); 
}

/*
 *************************************************************************
 *
 * Write a group to an existing path inside a root file
 *
 *************************************************************************
 */

void IOManager::writeGroupToRootFileAtPath(sidre::DataGroup * group,
                                           const std::string& file_name,
                                           const std::string& group_path)
{
  hid_t root_file_id = H5Fopen(file_name.c_str(),
                               H5F_ACC_RDWR,
                               H5P_DEFAULT);

  SLIC_ASSERT(root_file_id >= 0);

  hid_t path_id = H5Gopen(root_file_id, group_path.c_str(), 0);

  SLIC_ASSERT(path_id >= 0);

  hid_t group_id = H5Gcreate2(path_id,
                              group->getName().c_str(),
                              H5P_DEFAULT,
                              H5P_DEFAULT,
                              H5P_DEFAULT);

  SLIC_ASSERT(group_id >= 0);

  conduit::Node data_holder;
  group->createNativeLayout(data_holder);

  conduit::relay::io::hdf5_write(data_holder, group_id);

  herr_t errv = H5Gclose(group_id);
  SLIC_ASSERT(errv >= 0);

  errv = H5Fflush(root_file_id, H5F_SCOPE_LOCAL);
  SLIC_ASSERT(errv >= 0);

  errv =  H5Fclose(root_file_id);
  SLIC_ASSERT(errv >= 0);

}

/*
 *************************************************************************
 *
 * Write a view to an existing path inside a root file
 *
 *************************************************************************
 */

void IOManager::writeViewToRootFileAtPath(sidre::DataView * view,
                                          const std::string& file_name,
                                          const std::string& group_path)
{
  hid_t root_file_id = H5Fopen(file_name.c_str(),
                               H5F_ACC_RDWR,
                               H5P_DEFAULT);

  SLIC_ASSERT(root_file_id >= 0);

  hid_t path_id = H5Gopen(root_file_id, group_path.c_str(), 0);

  SLIC_ASSERT(path_id >= 0);

  conduit::Node data_holder;
  view->createNativeLayout(data_holder[view->getName()]);

  conduit::relay::io::hdf5_write(data_holder, path_id);

  herr_t errv = H5Fflush(root_file_id, H5F_SCOPE_LOCAL);
  SLIC_ASSERT(errv >= 0);

  errv =  H5Fclose(root_file_id);
  SLIC_ASSERT(errv >= 0);

}



} /* end namespace spio */
} /* end namespace asctoolkit */

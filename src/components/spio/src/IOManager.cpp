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
 * \file IOManager.cpp
 *
 * \brief   Implementation file for IOManager class.
 *
 ******************************************************************************
 */

// Associated header file
#include "IOManager.hpp"

// Other axom headers
#include "axom_utils/FileUtilities.hpp"

// SiDRe project headers
#include "sidre/Group.hpp"
#include "sidre/DataStore.hpp"
#include "sidre/SidreTypes.hpp"
#include "fmt/fmt.hpp"

// Conduit headers
#include "conduit_relay.hpp"
#include "conduit_relay_hdf5.hpp"

namespace axom
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
                     bool use_scr)
: m_comm_size(1),
  m_my_rank(0),
  m_baton(AXOM_NULLPTR),
  m_mpi_comm(comm),
  m_scr_initialized(false)
{
  MPI_Comm_size(comm, &m_comm_size);
  MPI_Comm_rank(comm, &m_my_rank);
#ifdef AXOM_USE_SCR
  if (use_scr) {
//    SCR_Init();
    m_scr_initialized = true;
  }
#endif
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
#ifdef AXOM_USE_SCR
  if (m_scr_initialized) {
//    SCR_Finalize();
  }
#endif
}


/*
 *************************************************************************
 *
 * Write to file.
 *
 *************************************************************************
 */
void IOManager::write(sidre::Group * datagroup, int num_files, const std::string& file_string, const std::string& protocol)
{
  if (m_baton) {
    if (m_baton->getNumFiles() != num_files) {
      delete m_baton;
      m_baton = AXOM_NULLPTR;
    }
  }

  if (!m_baton) {
    m_baton = new IOBaton(m_mpi_comm, num_files);
  }

  std::string root_string = file_string;
#ifdef AXOM_USE_SCR
  if (m_scr_initialized) { 
    SCR_Start_checkpoint();

//    char checkpoint_file[256];
//    sprintf(checkpoint_file, "%s_%6d", file_string.c_str(), m_my_rank);
//    char scr_file[SCR_MAX_FILENAME];
//    SCR_Route_file(checkpoint_file, scr_file);
//    SCR_Route_file(file_string.c_str(), scr_file);
//    root_string = std::string(scr_file);
  }
#endif
  if (m_my_rank == 0) {
    createRootFile(root_string, num_files, protocol);
  }
  MPI_Barrier(m_mpi_comm);
////**** BROADCAST m_checkpoint_dir to all procs!!!!!!

  std::string root_name = root_string + ".root";
  if (m_scr_initialized) {
    root_name = m_scr_checkpoint_dir + "/" + root_name;
  }

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

    herr_t status;
    AXOM_DEBUG_VAR(status);

    status = H5Gclose(h5_group_id);
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

#ifdef AXOM_USE_SCR
  if (m_scr_initialized) {
    MPI_Barrier(m_mpi_comm);
    int valid = 1;
    SCR_Complete_checkpoint(valid);
  }
#endif
}

/*
 *************************************************************************
 *
 * Read from files
 *
 *************************************************************************
 */
void IOManager::read(
  sidre::Group * datagroup,
  const std::string& root_file,
  const std::string& protocol,
  bool preserve_contents)
{
  if (protocol == "sidre_hdf5") {
    readSidreHDF5(datagroup, root_file, preserve_contents);
  } else {
    if (m_baton) {
      if (m_baton->getNumFiles() != 1) {
        delete m_baton;
        m_baton = AXOM_NULLPTR;
      }
    }

    if (!m_baton) {
      m_baton = new IOBaton(m_mpi_comm, 1);
    }

    int group_id = m_baton->wait();

    std::string file_name = getRankGroupFileName(root_file, group_id, protocol);

    datagroup->load(file_name, protocol, preserve_contents);

    (void)m_baton->pass();

  }
}

/*
 *************************************************************************
 *
 * Read based on HDF5 root file.
 *
 *************************************************************************
 */
void IOManager::read(sidre::Group * datagroup, const std::string& root_file, bool preserve_contents)
{
  std::string protocol = getProtocol(root_file);

  read(datagroup, root_file, protocol, preserve_contents);
}

/*
 *************************************************************************
 *
 * Read based on HDF5 root file.
 *
 *************************************************************************
 */
void IOManager::readSidreHDF5(sidre::Group * datagroup,
                              const std::string& root_file,
                              bool preserve_contents)
{
  int num_files = getNumFilesFromRoot(root_file);
  SLIC_ASSERT(num_files > 0);

  if (m_baton) {
    if (m_baton->getNumFiles() != num_files) {
      delete m_baton;
      m_baton = AXOM_NULLPTR;
    }
  }

  if (!m_baton) {
    m_baton = new IOBaton(m_mpi_comm, num_files);
  }

  std::string file_pattern = getHDF5FilePattern(root_file);

  int group_id = m_baton->wait();

  herr_t errv;
  AXOM_DEBUG_VAR(errv);

  std::string hdf5_name = getHDF5FileName(file_pattern, root_file, group_id);

  hid_t h5_file_id = H5Fopen(hdf5_name.c_str(),
                             H5F_ACC_RDONLY,
                             H5P_DEFAULT);
  SLIC_ASSERT(h5_file_id >= 0);

  std::string group_name = fmt::sprintf("datagroup_%07d", m_my_rank);
  hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
  SLIC_ASSERT(h5_group_id >= 0);

  datagroup->load(h5_group_id, "sidre_hdf5", preserve_contents);

  errv = H5Gclose(h5_group_id);
  SLIC_ASSERT(errv >= 0);

  errv = H5Fclose(h5_file_id);
  SLIC_ASSERT(errv >= 0);

  (void)m_baton->pass();
}


void IOManager::loadExternalData(sidre::Group * datagroup, const std::string& root_file)
{
  int num_files = getNumFilesFromRoot(root_file);
  SLIC_ASSERT(num_files > 0);

  if (m_baton) {
    if (m_baton->getNumFiles() != num_files) {
      delete m_baton;
      m_baton = AXOM_NULLPTR;
    }
  }

  if (!m_baton) {
    m_baton = new IOBaton(m_mpi_comm, num_files);
  }

  std::string file_pattern = getHDF5FilePattern(root_file);

  int group_id = m_baton->wait();

  herr_t errv;
  AXOM_DEBUG_VAR(errv);

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
void IOManager::createRootFile(const std::string& file_base,
                               int num_files,
                               const std::string& protocol)
{

  conduit::Node n;
  std::string root_file_name;
  std::string conduit_protocol;
  std::string local_file_base;

  if (protocol == "sidre_hdf5" || protocol == "conduit_hdf5") {

    n["number_of_files"] = num_files;
    if (protocol == "sidre_hdf5") {
      std::string next;
      std::string slash = "/";
      conduit::utils::rsplit_string(file_base, slash, local_file_base, next);
      n["file_pattern"] = local_file_base + slash + local_file_base + "_" + "%07d.hdf5";
    } else {
      n["file_pattern"] = file_base + "_" + "%07d.conduit_hdf5";
    }
    n["number_of_trees"] = m_comm_size;

    n["tree_pattern"] = "datagroup_%07d";
    n["protocol/name"] = protocol;
    n["protocol/version"] = "0.0";

    root_file_name = file_base + ".root";
    conduit_protocol = "hdf5";
//    conduit::relay::io::save(n, file_base + ".root", "hdf5");
  } else {

    n["number_of_files"] = num_files;
    n["file_pattern"] = file_base + "_" + "%07d." + protocol;
    n["number_of_trees"] = m_comm_size;

    n["tree_pattern"] = "datagroup_%07d";
    n["protocol/name"] = protocol;
    n["protocol/version"] = "0.0";

    root_file_name = file_base + ".json.root";
    if (protocol == "sidre_conduit_json") {
      conduit_protocol = "conduit_json";
    } else if (protocol == "sidre_json" || protocol == "conduit_bin") {
      conduit_protocol = "json";
    } else {
      conduit_protocol = protocol;
    }
  }

  if (!m_scr_initialized) {
    conduit::relay::io::save(n, root_file_name, conduit_protocol);
  } else {
#ifdef AXOM_USE_SCR
//    SCR_Start_checkpoint();

    if (protocol == "sidre_hdf5") {
      n["file_pattern"] = local_file_base + "_" + "%07d.hdf5"; 
    }

    std::string root_name = root_file_name;
    char checkpoint_file[256];
    sprintf(checkpoint_file, "%s", root_name.c_str());
    char scr_file[SCR_MAX_FILENAME];
    SCR_Route_file(checkpoint_file, scr_file);
    root_file_name = scr_file;

    std::string dir_name; 
    utilities::filesystem::getDirName(dir_name, root_file_name);
    if (!dir_name.empty()) {
      utilities::filesystem::makeDirsForPath(dir_name);
    }
    m_scr_checkpoint_dir = dir_name;
    conduit::relay::io::save(n, root_file_name, conduit_protocol);
//    int valid = 1;
//    SCR_Complete_checkpoint(valid);
#endif
  }


}

/*
 *************************************************************************
 *
 * Get namescheme pattern for a file holding Group data.
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
 * Get protocol from root file
 *
 *************************************************************************
 */
std::string IOManager::getProtocol(
  const std::string& root_name)
{
  std::string extension;
  std::string base;
  std::string dot = ".";
  conduit::utils::rsplit_string(root_name, dot, extension, base);

  // sidre_hdf5 protocol is ".root" others are "*.protocol.root"
  if( extension.find(dot) == std::string::npos )
  {
    extension = "hdf5";
  }
  else
  {
    std::string new_base = base;
    conduit::utils::rsplit_string(new_base, dot, extension, base);
  }


  std::string protocol;
  int buf_size = 0;
  if (m_my_rank == 0) {
    conduit::Node n;
    conduit::relay::io::load(root_name, extension, n);

    protocol = n["protocol/name"].as_string();

    buf_size = protocol.size() + 1;
  }

  MPI_Bcast(&buf_size, 1, MPI_INT, 0, m_mpi_comm);

  char protocol_buf[buf_size];
  if (m_my_rank == 0) {
    strcpy(protocol_buf, protocol.c_str());
  }

  MPI_Bcast(protocol_buf, buf_size, MPI_CHAR, 0, m_mpi_comm);

  if (m_my_rank != 0) {
    protocol = std::string(protocol_buf);
  }

  return protocol;
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


std::string IOManager::getRankGroupFileName(
  const std::string& root_name,
  int rankgroup_id,
  const std::string& protocol)
{

  std::string file_name = "file";
  if (protocol == "sidre_hdf5" || protocol == "conduit_hdf5") {
    std::string file_pattern = getHDF5FilePattern(root_name );
    file_name = getHDF5FileName(file_pattern, root_name , rankgroup_id);
  } else {
    conduit::Node n;
    std::string relay_protocol = protocol;
    if (protocol == "sidre_json" || protocol == "conduit_bin") {
      relay_protocol = "json";
    } else if (protocol == "sidre_conduit_json") {
      relay_protocol = "conduit_json";
    }
    conduit::relay::io::load(root_name, relay_protocol, n);
    std::string file_pattern = n["file_pattern"].as_string();
    file_name = getHDF5FileName(file_pattern, root_name, rankgroup_id);
  }

  return file_name;
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

void IOManager::writeGroupToRootFile(sidre::Group * group,
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

  herr_t errv;
  AXOM_DEBUG_VAR(errv);

  errv = H5Gclose(group_id);
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

void IOManager::writeGroupToRootFileAtPath(sidre::Group * group,
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

  herr_t errv;
  AXOM_DEBUG_VAR(errv);

  errv = H5Gclose(group_id);
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

void IOManager::writeViewToRootFileAtPath(sidre::View * view,
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

  herr_t errv;
  AXOM_DEBUG_VAR(errv);

  errv = H5Fflush(root_file_id, H5F_SCOPE_LOCAL);
  SLIC_ASSERT(errv >= 0);

  errv =  H5Fclose(root_file_id);
  SLIC_ASSERT(errv >= 0);

}



} /* end namespace spio */
} /* end namespace axom */

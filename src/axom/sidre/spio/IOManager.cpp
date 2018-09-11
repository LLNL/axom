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

// Associated header file
#include "IOManager.hpp"

// Other axom headers
#include "axom/core/utilities/FileUtilities.hpp"

// SiDRe project headers
#include "axom/sidre/core/Group.hpp"
#include "axom/sidre/core/DataStore.hpp"
#include "axom/sidre/core/SidreTypes.hpp"
#include "fmt/fmt.hpp"

// Conduit headers
#include "conduit_relay.hpp"

#ifdef AXOM_USE_HDF5
#include "conduit_relay_hdf5.hpp"
#include "hdf5.h"
#endif

#ifdef AXOM_USE_SCR
#include "scr.h"
#endif

namespace
{

/*!
 *  Utility function to broadcast a string from rank 0 to all other ranks
 */
std::string broadcastString(const std::string& str, MPI_Comm comm, int rank)
{
  // Find the string size and send to all ranks
  int buf_size = 0;
  if (rank == 0)
  {
    buf_size = str.size() + 1;
  }
  MPI_Bcast(&buf_size, 1, MPI_INT, 0, comm);

  // Allocate a buffer, copy and broadcast
  char* buf = new char[buf_size];
  if (rank == 0)
  {
    strcpy(buf, str.c_str());
  }
  MPI_Bcast(buf, buf_size, MPI_CHAR, 0, comm);

  // Create a string, delete the buffer and return
  std::string res(buf);
  delete[] buf;
  buf = nullptr;

  return res;
}

} // end anonymous namespace

namespace axom
{
namespace sidre
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
  m_baton(nullptr),
  m_mpi_comm(comm),
  m_use_scr(use_scr)
{
  MPI_Comm_size(comm, &m_comm_size);
  MPI_Comm_rank(comm, &m_my_rank);
#ifndef AXOM_USE_SCR
  if (m_use_scr)
  {
    SLIC_WARNING(
      "IOManager constructor called with use_scr = true, but Axom was not compiled with SCR. IOManager will operate without SCR.");
  }
  m_use_scr = false;
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
  if (m_baton)
  {
    delete m_baton;
  }
}


/*
 ******************************************************************************
 *
 * Returns the conduit relay protocol corresponding to the given sidre protocol
 *
 ******************************************************************************
 */
std::string IOManager::correspondingRelayProtocol(
  const std::string& sidre_protocol)
{
#ifdef AXOM_USE_HDF5
  const std::string DEFAULT_PROTOCOL = "hdf5";
#else
  const std::string DEFAULT_PROTOCOL = "json";
#endif

  if(sidre_protocol == "sidre_hdf5"
     || sidre_protocol == "conduit_hdf5")
  {
    return "hdf5";
  }
  else if(sidre_protocol == "sidre_json"
          || sidre_protocol == "conduit_bin"
          || sidre_protocol == "json")
  {
    return "json";
  }
  else if(sidre_protocol == "sidre_conduit_json"
          || sidre_protocol == "conduit_json")
  {
    return "conduit_json";
  }

  SLIC_WARNING("'" << sidre_protocol << "' is not a valid sidre protocol.");
  return DEFAULT_PROTOCOL;
}

/*
 *************************************************************************
 *
 * Write to file.
 *
 *************************************************************************
 */
void IOManager::write(sidre::Group* datagroup, int num_files,
                      const std::string& file_string,
                      const std::string& protocol)
{
  if (m_baton)
  {
    if (m_baton->getNumFiles() != num_files)
    {
      delete m_baton;
      m_baton = nullptr;
    }
  }

  if (!m_baton)
  {
    m_baton = new IOBaton(m_mpi_comm, num_files);
  }

  std::string root_string = file_string;
#ifdef AXOM_USE_SCR
  if (m_use_scr)
  {
    SCR_Start_checkpoint();
  }
#endif
  if (m_my_rank == 0)
  {
    createRootFile(root_string, num_files, protocol);
  }
  MPI_Barrier(m_mpi_comm);

  std::string root_name = root_string + ".root";
  if (m_use_scr)
  {
    m_scr_checkpoint_dir =
      broadcastString(m_scr_checkpoint_dir, m_mpi_comm, m_my_rank);
    root_name = m_scr_checkpoint_dir + "/" + root_name;
  }

  MPI_Barrier(m_mpi_comm);

  if (protocol == "sidre_hdf5")
  {
#ifdef AXOM_USE_HDF5
    std::string file_pattern = getHDF5FilePattern(root_name);

    int group_id = m_baton->wait();

    std::string hdf5_name =
      getFileNameForRank(file_pattern, root_name, group_id);

    hid_t h5_file_id, h5_group_id;
    if (m_baton->isFirstInGroup())
    {
      std::string dir_name;
      utilities::filesystem::getDirName(dir_name, hdf5_name);
      if (!dir_name.empty())
      {
        utilities::filesystem::makeDirsForPath(dir_name);
      }
      h5_file_id = conduit::relay::io::hdf5_create_file(hdf5_name);
    }
    else
    {
      h5_file_id =
        conduit::relay::io::hdf5_open_file_for_read_write(hdf5_name);
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
#else
    SLIC_WARNING("'sidre_hdf5' protocol only available "
                 << "when axom is configured with hdf5");
#endif /* AXOM_USE_HDF5 */
  }
  else
  {
    int group_id = m_baton->wait();
    std::string file_name = fmt::sprintf("%s_%07d", file_string, group_id);

    std::string obase = file_name + "." + protocol;
    datagroup->save(obase, protocol);
  }
  (void)m_baton->pass();

#ifdef AXOM_USE_SCR
  if (m_use_scr)
  {
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
  sidre::Group* datagroup,
  const std::string& root_file,
  const std::string& protocol,
  bool preserve_contents)
{
  if (protocol == "sidre_hdf5")
  {
#ifdef AXOM_USE_HDF5
    readSidreHDF5(datagroup, root_file, preserve_contents);
#else
    SLIC_WARNING("'sidre_hdf5' protocol only available "
                 << "when axom is configured with hdf5");
#endif /* AXOM_USE_HDF5 */
  }
  else
  {
    if (m_baton)
    {
      if (m_baton->getNumFiles() != 1)
      {
        delete m_baton;
        m_baton = nullptr;
      }
    }

    if (!m_baton)
    {
      m_baton = new IOBaton(m_mpi_comm, m_comm_size);
    }

    std::string file_pattern = getFilePatternFromRoot(root_file, protocol);

    int group_id = m_baton->wait();

    std::string file_name =
      getFileNameForRank(file_pattern, root_file, group_id);

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
void IOManager::read(sidre::Group* datagroup,
                     const std::string& root_file,
                     bool preserve_contents,
                     bool read_with_scr)
{

  if (!read_with_scr)
  {
    std::string protocol = getProtocol(root_file);
    read(datagroup, root_file, protocol, preserve_contents);
  }
  else
  {
#ifdef AXOM_USE_SCR

    if (m_use_scr)
    {
      int valid = -1;
      SCR_Have_restart(&valid, 0);
      if (valid == 1)
      {
        SCR_Start_restart(0);
        readWithSCR(datagroup, root_file, preserve_contents);
        SCR_Complete_restart(valid);
      }
      else
      {
        SLIC_WARNING(
          "IOManager::read() requested to read files using SCR, but no SCR restart was found. Will attempt to restart without SCR.");
        read(datagroup, root_file, preserve_contents, false);
      }
    }
    else
    {
      SLIC_WARNING(
        "IOManager::read() requested to read files using SCR, but IOManager was constructed to not use SCR. SCR will not be used");
      read(datagroup, root_file, preserve_contents, false);
    }

#else

    SLIC_WARNING(
      "IOManager::read() requested to read files using SCR, but Axom was not compiled with SCR. SCR will not be used");
    read(datagroup, root_file, preserve_contents, false);

#endif
  }
}

/*
 *************************************************************************
 *
 * Read based on root file that was dumped in an SCR checkpoint.
 *
 *************************************************************************
 */
#ifdef AXOM_USE_SCR
void IOManager::readWithSCR(
  sidre::Group* datagroup,
  const std::string& root_file,
  bool preserve_contents)
{
  SLIC_ASSERT(m_use_scr);
  char file[SCR_MAX_FILENAME];
  if (SCR_Route_file(root_file.c_str(), file) == SCR_SUCCESS)
  {
    std::string scr_root(file);
    std::string protocol = getProtocol(scr_root);
    read(datagroup, scr_root, protocol, preserve_contents);
  }
  else
  {
    SLIC_WARNING(
      "Root file: "
      << root_file << " not found by SCR. Attempting to read without SCR.");

    read(datagroup, root_file, preserve_contents, false);
  }
}
#endif

void IOManager::loadExternalData(sidre::Group* datagroup,
                                 const std::string& root_file)
{
  int num_files = getNumFilesFromRoot(root_file);
  SLIC_ASSERT(num_files > 0);

  if (m_baton)
  {
    if (m_baton->getNumFiles() != num_files)
    {
      delete m_baton;
      m_baton = nullptr;
    }
  }

  if (!m_baton)
  {
    m_baton = new IOBaton(m_mpi_comm, num_files);
  }

#ifdef AXOM_USE_HDF5
  std::string file_pattern = getHDF5FilePattern(root_file);

  int group_id = m_baton->wait();

  herr_t errv;
  AXOM_DEBUG_VAR(errv);

  std::string hdf5_name = getFileNameForRank(file_pattern, root_file, group_id);

  hid_t h5_file_id = conduit::relay::io::hdf5_open_file_for_read(hdf5_name);
  SLIC_ASSERT(h5_file_id >= 0);

  std::string group_name = fmt::sprintf("datagroup_%07d", m_my_rank);
  hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
  SLIC_ASSERT(h5_group_id >= 0);

  datagroup->loadExternalData(h5_group_id);

  errv = H5Gclose(h5_group_id);
  SLIC_ASSERT(errv >= 0);

  errv = H5Fclose(h5_file_id);
  SLIC_ASSERT(errv >= 0);
#else
  AXOM_DEBUG_VAR(datagroup);
  SLIC_WARNING("Loading external data only only available "
               << "when Axom is configured with hdf5");
#endif /* AXOM_USE_HDF5 */


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
  std::string local_file_base;

  std::string relay_protocol = correspondingRelayProtocol(protocol);

  if (protocol == "sidre_hdf5" || protocol == "conduit_hdf5")
  {
#ifdef AXOM_USE_HDF5
    n["number_of_files"] = num_files;
    if (protocol == "sidre_hdf5")
    {
      std::string next;
      std::string slash = "/";
      conduit::utils::rsplit_string(file_base, slash, local_file_base, next);
      n["file_pattern"] = local_file_base + slash + local_file_base + "_" +
                          "%07d.hdf5";
    }
    else
    {
      n["file_pattern"] = file_base + "_" + "%07d.conduit_hdf5";
    }
    n["number_of_trees"] = m_comm_size;

    n["tree_pattern"] = "datagroup_%07d";
    n["protocol/name"] = protocol;
    n["protocol/version"] = "0.0";

    root_file_name = file_base + ".root";
#else
    SLIC_WARNING("IOManager::createRootFile() -- '"
                 << protocol <<"' protocol only available "
                 << "when Axom is configured with hdf5");
#endif /* AXOM_USE_HDF5 */
  }
  else
  {

    n["number_of_files"] = num_files;
    n["file_pattern"] = file_base + "_" + "%07d." + protocol;
    n["number_of_trees"] = m_comm_size;

    n["tree_pattern"] = "datagroup_%07d";
    n["protocol/name"] = protocol;
    n["protocol/version"] = "0.0";

    root_file_name = file_base + ".root";
  }

#ifdef AXOM_USE_SCR
  if (m_use_scr)
  {
    if (protocol == "sidre_hdf5")
    {
      n["file_pattern"] = local_file_base + "_" + "%07d.hdf5";
    }

    std::string root_name = root_file_name;
    char checkpoint_file[256];
    sprintf(checkpoint_file, "%s", root_name.c_str());
    char scr_file[SCR_MAX_FILENAME];
    if (SCR_Route_file(checkpoint_file, scr_file) == SCR_SUCCESS)
    {
      root_file_name = scr_file;
    }
    else
    {
      SLIC_WARNING("Attempt to create SCR route for file: " << root_name <<
                   " failed. Writing root file without SCR.");
    }

    std::string dir_name;
    utilities::filesystem::getDirName(dir_name, root_file_name);
    if (!dir_name.empty())
    {
      utilities::filesystem::makeDirsForPath(dir_name);
    }
    m_scr_checkpoint_dir = dir_name;
  }
#endif

  conduit::relay::io::save(n, root_file_name, relay_protocol);


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

  // Separate the first extension from the root file name
  // It should always be "root"
  conduit::utils::rsplit_string(root_name, dot, extension, base);

  SLIC_CHECK_MSG(extension == "root",
                 "The root file name should always end in 'root'."
                 << " File name was '"<< root_name <<"'");

  std::string relay_protocol = "json";
#ifdef AXOM_USE_HDF5
  // Attempt to open the root file using HDF5.  If it succeeds, set
  // relay_protocol to "hdf5", otherwise we assume a json protocol.
  //
  // Suppress error output for H5Fopen, since failure is acceptable here.
  H5E_auto2_t herr_func;
  void* old_client_data;
  H5Eget_auto(H5E_DEFAULT, &herr_func, &old_client_data);
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  hid_t file_id = H5Fopen(root_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id > 0)
  {
    relay_protocol = "hdf5";
    herr_t errv = H5Fclose(file_id);
    SLIC_ASSERT(errv >= 0);
  }

  // Restore error output
  H5Eset_auto(H5E_DEFAULT, herr_func, old_client_data);
#endif

  std::string protocol;
  if (m_my_rank == 0)
  {
    conduit::Node n;
    conduit::relay::io::load(root_name, relay_protocol, n);

    conduit::Node& pnode = n["protocol/name"];
    if (pnode.dtype().is_string())
    {
      protocol = n["protocol/name"].as_string();
    }
    else if (relay_protocol == "json")
    {
      //If relay_protocol "json" didn't work, try "conduit_json"
      n.reset();
      conduit::relay::io::load(root_name, "conduit_json", n);
      protocol = n["protocol/name"].as_string();
    }

    if (protocol.empty())
    {
      // Did not find protocol name, issue warning and assign a default guess.

      SLIC_WARNING(
        "'" << root_name
            << "/protocol/name' does not contain a valid Sidre protocol name.  "
            << "Will attempt to use a default protocol.");

      if (relay_protocol == "hdf5")
      {
        protocol = "sidre_hdf5";
      }
      else
      {
        protocol = "sidre_json";
      }
    }
  }

  protocol = broadcastString(protocol, m_mpi_comm, m_my_rank);
  return protocol;
}



std::string IOManager::getFilePatternFromRoot(const std::string& root_name,
                                              const std::string& protocol)
{
  std::string file_pattern;
  if (m_my_rank == 0)
  {
    conduit::Node n;
    std::string relay_protocol = correspondingRelayProtocol(protocol);
    conduit::relay::io::load(root_name, relay_protocol, n);
    file_pattern = n["file_pattern"].as_string();
  }

  file_pattern = broadcastString(file_pattern, m_mpi_comm, m_my_rank);
  return file_pattern;
}

// Several private functions are only relevant when HDF5 is enabled
#ifdef AXOM_USE_HDF5

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
  if (m_my_rank == 0)
  {
    conduit::Node n;
    conduit::relay::io::load(root_name + ":file_pattern", "hdf5", n);

    file_pattern = n.as_string();
  }
  file_pattern = broadcastString(file_pattern, m_mpi_comm, m_my_rank);

  return file_pattern;
}

/*
 *************************************************************************
 *
 * Read based on HDF5 root file.
 *
 *************************************************************************
 */
void IOManager::readSidreHDF5(sidre::Group* datagroup,
                              const std::string& root_file,
                              bool preserve_contents)
{
  int num_files = getNumFilesFromRoot(root_file);
  SLIC_ASSERT(num_files > 0);

  if (m_baton)
  {
    if (m_baton->getNumFiles() != num_files)
    {
      delete m_baton;
      m_baton = nullptr;
    }
  }

  if (!m_baton)
  {
    m_baton = new IOBaton(m_mpi_comm, num_files);
  }

  std::string file_pattern = getHDF5FilePattern(root_file);

  int group_id = m_baton->wait();

  herr_t errv;
  AXOM_DEBUG_VAR(errv);

  std::string hdf5_name = getFileNameForRank(file_pattern, root_file, group_id);

  hid_t h5_file_id = conduit::relay::io::hdf5_open_file_for_read(hdf5_name);
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
#endif /* AXOM_USE_HDF5 */

/*
 *************************************************************************
 *
 * Get the file name based on a pattern and a group id
 *
 *************************************************************************
 */
std::string IOManager::getFileNameForRank(
  const std::string& file_pattern,
  const std::string& root_name,
  int rankgroup_id)
{
  std::string file_name = fmt::sprintf(file_pattern.c_str(), rankgroup_id);

  //If the root file was given as a path,
  //find the directory and add it to file_name
  std::string curr;
  std::string root_dir;
  std::string slash = "/";
  conduit::utils::rsplit_string(root_name, slash, curr, root_dir);

  if (!root_dir.empty())
  {
    file_name = root_dir + slash + file_name;
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
  if (m_my_rank == 0)
  {
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

void IOManager::writeGroupToRootFile(sidre::Group* group,
                                     const std::string& file_name)
{
#ifdef AXOM_USE_HDF5
  hid_t root_file_id =
    conduit::relay::io::hdf5_open_file_for_read_write(file_name);

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
#else
  AXOM_DEBUG_VAR(group);
  AXOM_DEBUG_VAR(file_name);

  SLIC_WARNING(
    "Axom configured without hdf5. "
    <<"IOManager::writeGroupToRootFile() only currently implemented "
    <<"for 'sidre_hdf5' protocol. ");
#endif /* AXOM_USE_HDF5 */
}

/*
 *************************************************************************
 *
 * Write a group to an existing path inside a root file
 *
 *************************************************************************
 */

void IOManager::writeGroupToRootFileAtPath(sidre::Group* group,
                                           const std::string& file_name,
                                           const std::string& group_path)
{
#ifdef AXOM_USE_HDF5
  hid_t root_file_id =
    conduit::relay::io::hdf5_open_file_for_read_write(file_name);

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
#else
  AXOM_DEBUG_VAR(group);
  AXOM_DEBUG_VAR(file_name);
  AXOM_DEBUG_VAR(group_path);

  SLIC_WARNING("Axom configured without hdf5. "
               <<"IOManager::writeGroupToRootFileAtPath() only currently "
               <<"implemented for 'sidre_hdf5' protocol. ");
#endif /* AXOM_USE_HDF5 */
}

/*
 *************************************************************************
 *
 * Write a view to an existing path inside a root file
 *
 *************************************************************************
 */

void IOManager::writeViewToRootFileAtPath(sidre::View* view,
                                          const std::string& file_name,
                                          const std::string& group_path)
{
#ifdef AXOM_USE_HDF5
  hid_t root_file_id =
    conduit::relay::io::hdf5_open_file_for_read_write(file_name);

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
#else
  AXOM_DEBUG_VAR(view);
  AXOM_DEBUG_VAR(file_name);
  AXOM_DEBUG_VAR(group_path);

  SLIC_WARNING("Axom configured without hdf5. "
               <<"writeViewToRootFileAtPath() only currently implemented "
               <<"for 'sidre_hdf5' protocol. ");
#endif /* AXOM_USE_HDF5 */
}

void IOManager::writeBlueprintIndexToRootFile(DataStore* datastore,
                                              const std::string& domain_path,
                                              const std::string& file_name,
                                              const std::string& mesh_name)
{
#ifdef AXOM_USE_HDF5
  hid_t root_file_id =
    conduit::relay::io::hdf5_open_file_for_read_write(file_name);

  SLIC_ASSERT(root_file_id >= 0);

  std::string bp_index("blueprint_index/" + mesh_name); 
  datastore->generateBlueprintIndex(domain_path,
     mesh_name, bp_index, 1);

  Group* ind_group = datastore->getRoot()->getGroup("blueprint_index");
  writeGroupToRootFile(ind_group, file_name);

#else
  AXOM_DEBUG_VAR(datastore);
  AXOM_DEBUG_VAR(domain_path);
  AXOM_DEBUG_VAR(file_name);
  AXOM_DEBUG_VAR(index_path);
#endif /* AXOM_USE_HDF5 */
}



} /* end namespace sidre */
} /* end namespace axom */

// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

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
#include "conduit_relay_io_hdf5.hpp"
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
    m_baton = new IOBaton(m_mpi_comm, num_files, m_comm_size);
  }

  SLIC_ERROR_IF(m_use_scr && num_files != m_comm_size,
                "SCR requires a file per process");

  std::string root_string = file_string;
  if (m_my_rank == 0)
  {
    createRootFile(root_string, num_files, protocol);
  }
  MPI_Barrier(m_mpi_comm);

  std::string root_name = root_string + ".root";

  if (protocol == "sidre_hdf5")
  {
#ifdef AXOM_USE_HDF5
    std::string file_pattern = getHDF5FilePattern(root_name);

    int set_id = m_baton->wait();

    std::string hdf5_name =
      getFileNameForRank(file_pattern, root_name, set_id);

    if (m_use_scr)
    {
      hdf5_name = getSCRPath(hdf5_name);
    }

    hid_t h5_file_id, h5_group_id;
    if (m_baton->isFirstInGroup())
    {
      // no need to create directories in SCR
      if (!m_use_scr)
      {
        std::string dir_name;
        utilities::filesystem::getDirName(dir_name, hdf5_name);
        if (!dir_name.empty())
        {
          utilities::filesystem::makeDirsForPath(dir_name);
        }
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
    int set_id = m_baton->wait();
    std::string file_name = fmt::sprintf("%s_%07d", file_string, set_id);

    std::string obase = file_name + "." + protocol;
    datagroup->save(obase, protocol);
  }
  (void)m_baton->pass();

  MPI_Barrier(m_mpi_comm);
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
  MPI_Barrier(m_mpi_comm);

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
      m_baton = new IOBaton(m_mpi_comm, m_comm_size, m_comm_size);
    }

    std::string file_pattern = getFilePatternFromRoot(root_file, protocol);

    int set_id = m_baton->wait();

    std::string file_name =
      getFileNameForRank(file_pattern, root_file, set_id);

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
                     bool preserve_contents)
{
  MPI_Barrier(m_mpi_comm);

#ifdef AXOM_USE_SCR
  if (m_use_scr)
  {
    readWithSCR(datagroup, root_file, preserve_contents);
  }
  else
#endif
  {
    std::string protocol = getProtocol(root_file);
    read(datagroup, root_file, protocol, preserve_contents);
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

  std::string scr_root = root_file;
  if (m_my_rank == 0)
  {
    scr_root = getSCRPath(scr_root);
  }
  std::string protocol = getProtocol(scr_root);

  read(datagroup, root_file, protocol, preserve_contents);
}
#endif

std::string IOManager::getSCRPath(const std::string & path)
{
#ifdef AXOM_USE_SCR
  SLIC_ASSERT(m_use_scr);
  char scr_name[SCR_MAX_FILENAME];
  std::string scr_path;
  if (SCR_Route_file(path.c_str(), scr_name) == SCR_SUCCESS)
  {
    scr_path = scr_name;
  }
  else
  {
    SLIC_WARNING(
      "SCR routing of path "
      << path << " was unsuccessful. Attempting to proceed without routing.");
    scr_path = path;
  }
  return scr_path;
#else
  return path;
#endif
}

void IOManager::loadExternalData(sidre::Group* datagroup,
                                 const std::string& root_file)
{
  int num_files = getNumFilesFromRoot(root_file);
  int num_groups = getNumGroupsFromRoot(root_file);
  SLIC_ASSERT(num_files > 0);
  SLIC_ASSERT(num_groups > 0);

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
    m_baton = new IOBaton(m_mpi_comm, num_files, num_groups);
  }

#ifdef AXOM_USE_HDF5
  std::string file_pattern = getHDF5FilePattern(root_file);

  int set_id = m_baton->wait();

  if (num_groups <= m_comm_size)
  {
    if (m_my_rank < num_groups)
    {
      herr_t errv;
      AXOM_DEBUG_VAR(errv);

      std::string hdf5_name =
        getFileNameForRank(file_pattern, root_file, set_id);

      if (m_use_scr)
      {
        hdf5_name = getSCRPath(hdf5_name);
      }

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
    }
  }
  else
  {
    for (int input_rank = m_my_rank ; input_rank < num_groups ;
         input_rank += m_comm_size)
    {

      herr_t errv;
      AXOM_DEBUG_VAR(errv);

      std::string hdf5_name =
        getFileNameForRank(file_pattern, root_file, input_rank);

      if (m_use_scr)
      {
        hdf5_name = getSCRPath(hdf5_name);
      }

      hid_t h5_file_id = conduit::relay::io::hdf5_open_file_for_read(hdf5_name);
      SLIC_ASSERT(h5_file_id >= 0);

      std::string group_name = fmt::sprintf("datagroup_%07d", input_rank);
      hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
      SLIC_ASSERT(h5_group_id >= 0);

      std::string input_name =
        fmt::sprintf("rank_%07d/sidre_input", input_rank);
      Group* one_rank_input = datagroup->getGroup(input_name);

      one_rank_input->loadExternalData(h5_group_id);

      errv = H5Gclose(h5_group_id);
      SLIC_ASSERT(errv >= 0);

      errv = H5Fclose(h5_file_id);
      SLIC_ASSERT(errv >= 0);

    }
  }

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
  SLIC_ASSERT(m_my_rank == 0);

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
    root_file_name = getSCRPath(root_file_name);
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
  std::string protocol;

  if (m_my_rank == 0)
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
    H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);

    hid_t file_id = H5Fopen(root_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id > 0)
    {
      relay_protocol = "hdf5";
      herr_t errv = H5Fclose(file_id);
      AXOM_DEBUG_VAR(errv);
      SLIC_ASSERT(errv >= 0);
    }

    // Restore error output
    H5Eset_auto(H5E_DEFAULT, herr_func, old_client_data);

#endif

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
    std::string root_path = root_name;
    if (m_use_scr)
    {
      root_path = getSCRPath(root_path);
    }

    conduit::Node n;
    std::string relay_protocol = correspondingRelayProtocol(protocol);
    conduit::relay::io::load(root_path, relay_protocol, n);
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
    std::string root_path = root_name;
    if (m_use_scr)
    {
      root_path = getSCRPath(root_path);
    }

    conduit::Node n;
    conduit::relay::io::load(root_path + ":file_pattern", "hdf5", n);

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
  int num_groups = getNumGroupsFromRoot(root_file);
  SLIC_ASSERT(num_files > 0);
  SLIC_ERROR_IF(num_groups > m_comm_size && num_groups != num_files,
                "IOManager attempted to read using a smaller number of processors "
                << "than were used to produce the I/O files.  This only can work if "
                << "those files were created in file-per-processor mode.");

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
    m_baton = new IOBaton(m_mpi_comm, num_files, num_groups);
  }

  std::string file_pattern = getHDF5FilePattern(root_file);

  int set_id = m_baton->wait();
  if (num_groups <= m_comm_size)
  {
    if (m_my_rank < num_groups)
    {

      herr_t errv;
      AXOM_DEBUG_VAR(errv);

      std::string hdf5_name =
        getFileNameForRank(file_pattern, root_file, set_id);

      if (m_use_scr)
      {
        hdf5_name = getSCRPath(hdf5_name);
      }

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
    }
  }
  else
  {

    datagroup->createViewScalar("reduced_input_ranks", num_groups);

    for (int input_rank = m_my_rank ; input_rank < num_groups ;
         input_rank += m_comm_size)
    {

      herr_t errv;
      AXOM_DEBUG_VAR(errv);

      std::string hdf5_name =
        getFileNameForRank(file_pattern, root_file, input_rank);

      if (m_use_scr)
      {
        hdf5_name = getSCRPath(hdf5_name);
      }

      hid_t h5_file_id = conduit::relay::io::hdf5_open_file_for_read(hdf5_name);
      SLIC_ASSERT(h5_file_id >= 0);

      std::string group_name = fmt::sprintf("datagroup_%07d", input_rank);
      hid_t h5_group_id = H5Gopen(h5_file_id, group_name.c_str(), 0);
      SLIC_ASSERT(h5_group_id >= 0);

      std::string input_name =
        fmt::sprintf("rank_%07d/sidre_input", input_rank);
      Group* one_rank_input = datagroup->createGroup(input_name);

      one_rank_input->load(h5_group_id, "sidre_hdf5", preserve_contents);

      errv = H5Gclose(h5_group_id);
      SLIC_ASSERT(errv >= 0);

      errv = H5Fclose(h5_file_id);
      SLIC_ASSERT(errv >= 0);
    }
  }
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
    std::string root_path = root_file;
    if (m_use_scr)
    {
      root_path = getSCRPath(root_path);
    }

    conduit::Node n;
    conduit::relay::io::load(root_path + ":number_of_files","hdf5",n);
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

int IOManager::getNumGroupsFromRoot(const std::string& root_file)
{
  /*
   * Read number_of_trees from rootfile on rank 0.
   */
  int read_num_trees = 0;
  if (m_my_rank == 0)
  {
    std::string root_path = root_file;
    if (m_use_scr)
    {
      root_path = getSCRPath(root_path);
    }

    conduit::Node n;
    conduit::relay::io::load(root_path + ":number_of_trees","hdf5",n);
    read_num_trees = n.to_int();
    SLIC_ASSERT(read_num_trees > 0);
  }

  /*
   * Reduction sets num_groups on all ranks.
   */
  int num_groups;
  MPI_Allreduce(&read_num_trees, &num_groups, 1, MPI_INT, MPI_SUM, m_mpi_comm);
  SLIC_ASSERT(num_groups > 0);

  return num_groups;
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
  std::string tmp_name = file_name;
  if (m_use_scr)
  {
    tmp_name = getSCRPath(tmp_name);
  }

  hid_t root_file_id =
    conduit::relay::io::hdf5_open_file_for_read_write(tmp_name);

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
  std::string tmp_name = file_name;
  if (m_use_scr)
  {
    tmp_name = getSCRPath(tmp_name);
  }

  hid_t root_file_id =
    conduit::relay::io::hdf5_open_file_for_read_write(tmp_name);

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
  std::string tmp_name = file_name;
  if (m_use_scr)
  {
    tmp_name = getSCRPath(tmp_name);
  }

  hid_t root_file_id =
    conduit::relay::io::hdf5_open_file_for_read_write(tmp_name);

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
                                              const std::string& mesh_path)
{
#ifdef AXOM_USE_HDF5
  std::string tmp_name = file_name;
  if (m_use_scr)
  {
    tmp_name = getSCRPath(tmp_name);
  }

  hid_t root_file_id =
    conduit::relay::io::hdf5_open_file_for_read_write(tmp_name);

  AXOM_DEBUG_VAR(root_file_id);
  SLIC_ASSERT(root_file_id >= 0);

  std::string blueprint_name;
  std::string path_to_mesh;
  std::string delimiter(1, datastore->getRoot()->getPathDelimiter());

  //The final name in mesh_path will be used as the name of the
  //blueprint index.
  conduit::utils::rsplit_string(mesh_path, delimiter,
                                blueprint_name, path_to_mesh);

  std::string bp_index("blueprint_index/" + blueprint_name);

  bool success = datastore->generateBlueprintIndex(domain_path,
                                                   mesh_path, bp_index,
                                                   m_comm_size);

  if (success)
  {
    Group* ind_group = datastore->getRoot()->getGroup("blueprint_index");
    writeGroupToRootFile(ind_group, file_name);
  }
  else
  {
    SLIC_WARNING("DataStore failed to generate Blueprint Index "
                 <<"based on group at path "
                 << domain_path);
  }
#else
  AXOM_DEBUG_VAR(datastore);
  AXOM_DEBUG_VAR(domain_path);
  AXOM_DEBUG_VAR(file_name);
  AXOM_DEBUG_VAR(mesh_path);

  SLIC_WARNING("Axom configured without hdf5. "
               <<"writeBlueprintIndexToRootFile() only currently implemented "
               <<"for 'sidre_hdf5' protocol. ");
#endif /* AXOM_USE_HDF5 */
}



} /* end namespace sidre */
} /* end namespace axom */

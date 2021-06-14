// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"

// Slic includes
#include "axom/slic/interface/slic.hpp"  // for SLIC macros
#include "axom/slic/streams/GenericOutputStream.hpp"
#if defined(AXOM_USE_MPI) && defined(AXOM_USE_LUMBERJACK)
  #include "axom/slic/streams/LumberjackStream.hpp"
#elif defined(AXOM_USE_MPI) && !defined(AXOM_USE_LUMBERJACK)
  #include "axom/slic/streams/SynchronizedStream.hpp"
#endif

// Mint includes
#include "axom/mint/mesh/UnstructuredMesh.hpp"  // for mint::UnstructuredMesh

// Quest includes
#ifdef AXOM_USE_MPI
  #include "axom/quest/stl/PSTLReader.hpp"
#endif

#include <limits>

namespace axom
{
namespace quest
{
namespace internal
{
/// MPI Helper/Wrapper Methods

#ifdef AXOM_USE_MPI

/*
 * Deallocates the specified MPI window object.
 */
void mpi_win_free(MPI_Win* window)
{
  if(*window != MPI_WIN_NULL)
  {
    MPI_Win_free(window);
  }
}

/*
 * Deallocates the specified MPI communicator object.
 */
void mpi_comm_free(MPI_Comm* comm)
{
  if(*comm != MPI_COMM_NULL)
  {
    MPI_Comm_free(comm);
  }
}

/*
 * Reads the mesh on rank 0 and exchanges the mesh metadata, i.e., the
 * number of nodes and faces with all other ranks.
 */
int read_and_exchange_mesh_metadata(int global_rank_id,
                                    MPI_Comm global_comm,
                                    quest::STLReader& reader,
                                    axom::IndexType mesh_metadata[2])
{
  constexpr int NUM_NODES = 0;
  constexpr int NUM_FACES = 1;
  constexpr int ROOT_RANK = 0;

  switch(global_rank_id)
  {
  case 0:
    if(reader.read() == READ_SUCCESS)
    {
      mesh_metadata[NUM_NODES] = reader.getNumNodes();
      mesh_metadata[NUM_FACES] = reader.getNumFaces();
    }
    else
    {
      SLIC_WARNING("reading STL file failed, setting mesh to NULL");
      mesh_metadata[NUM_NODES] = READ_FAILED;
      mesh_metadata[NUM_FACES] = READ_FAILED;
    }
    MPI_Bcast(mesh_metadata,
              2,
              axom::mpi_traits<axom::IndexType>::type,
              ROOT_RANK,
              global_comm);
    break;
  default:
    MPI_Bcast(mesh_metadata,
              2,
              axom::mpi_traits<axom::IndexType>::type,
              ROOT_RANK,
              global_comm);
  }

  int rc = (mesh_metadata[NUM_NODES] == READ_FAILED) ? READ_FAILED : READ_SUCCESS;
  return rc;
}

#endif /* AXOM_USE_MPI */

#ifdef AXOM_USE_MPI3
/*
 * Creates inter-node and intra-node communicators from the given global
 * MPI communicator handle.
 */
void create_communicators(MPI_Comm global_comm,
                          MPI_Comm& intra_node_comm,
                          MPI_Comm& inter_node_comm,
                          int& global_rank_id,
                          int& local_rank_id,
                          int& intercom_rank_id)
{
  // Sanity checks
  SLIC_ASSERT(global_comm != MPI_COMM_NULL);
  SLIC_ASSERT(intra_node_comm == MPI_COMM_NULL);
  SLIC_ASSERT(inter_node_comm == MPI_COMM_NULL);

  constexpr int IGNORE_KEY = 0;

  // STEP 0: get global rank, used to order ranks in the inter-node comm.
  MPI_Comm_rank(global_comm, &global_rank_id);

  // STEP 1: create the intra-node communicator
  MPI_Comm_split_type(global_comm,
                      MPI_COMM_TYPE_SHARED,
                      IGNORE_KEY,
                      MPI_INFO_NULL,
                      &intra_node_comm);
  MPI_Comm_rank(intra_node_comm, &local_rank_id);
  SLIC_ASSERT(local_rank_id >= 0);

  // STEP 2: create inter-node communicator
  const int color = (local_rank_id == 0) ? 1 : MPI_UNDEFINED;
  MPI_Comm_split(global_comm, color, global_rank_id, &inter_node_comm);

  if(color == 1)
  {
    MPI_Comm_rank(inter_node_comm, &intercom_rank_id);
  }

  SLIC_ASSERT(intra_node_comm != MPI_COMM_NULL);
}
#endif

#if defined(AXOM_USE_MPI) && defined(AXOM_USE_MPI3)
/*
 * Allocates a shared memory buffer for the mesh that is shared among
 * all the ranks within the same compute node.
 */
MPI_Aint allocate_shared_buffer(int local_rank_id,
                                MPI_Comm intra_node_comm,
                                const axom::IndexType mesh_metadata[2],
                                double*& x,
                                double*& y,
                                double*& z,
                                axom::IndexType*& conn,
                                unsigned char*& mesh_buffer,
                                MPI_Win& shared_window)
{
  constexpr int ROOT_RANK = 0;

  const int nnodes = mesh_metadata[0];
  const int nfaces = mesh_metadata[1];

  int disp = sizeof(unsigned char);
  MPI_Aint bytesize =
    nnodes * 3 * sizeof(double) + nfaces * 3 * sizeof(axom::IndexType);
  MPI_Aint window_size = (local_rank_id != ROOT_RANK) ? 0 : bytesize;

  MPI_Win_allocate_shared(window_size,
                          disp,
                          MPI_INFO_NULL,
                          intra_node_comm,
                          &mesh_buffer,
                          &shared_window);
  MPI_Win_shared_query(shared_window, ROOT_RANK, &bytesize, &disp, &mesh_buffer);

  // calculate offset to the coordinates & cell connectivity in the buffer
  int baseOffset = nnodes * sizeof(double);
  int x_offset = 0;
  int y_offset = baseOffset;
  int z_offset = y_offset + baseOffset;
  int conn_offset = z_offset + baseOffset;

  x = reinterpret_cast<double*>(&mesh_buffer[x_offset]);
  y = reinterpret_cast<double*>(&mesh_buffer[y_offset]);
  z = reinterpret_cast<double*>(&mesh_buffer[z_offset]);
  conn = reinterpret_cast<axom::IndexType*>(&mesh_buffer[conn_offset]);

  return (bytesize);
}
#endif

/// Mesh I/O methods

#if defined(AXOM_USE_MPI) && defined(AXOM_USE_MPI3)
/*
 * Reads in the surface mesh from the specified file into a shared
 * memory buffer that is attached to the given MPI shared window.
 */
int read_mesh_shared(const std::string& file,
                     MPI_Comm global_comm,
                     unsigned char*& mesh_buffer,
                     mint::Mesh*& m,
                     MPI_Comm& intra_node_comm,
                     MPI_Win& shared_window)
{
  SLIC_ASSERT(global_comm != MPI_COMM_NULL);
  SLIC_ASSERT(intra_node_comm == MPI_COMM_NULL);
  SLIC_ASSERT(shared_window == MPI_WIN_NULL);

  // NOTE: STL meshes are always 3D mesh consisting of triangles.
  using TriangleMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

  // STEP 0: check input mesh pointer
  if(m != nullptr)
  {
    SLIC_WARNING("supplied mesh pointer is not null!");
    return READ_FAILED;
  }

  if(mesh_buffer != nullptr)
  {
    SLIC_WARNING("supplied mesh buffer should be null!");
    return READ_FAILED;
  }

  // STEP 1: create intra-node and inter-node MPI communicators
  int global_rank_id = -1;
  int local_rank_id = -1;
  int intercom_rank_id = -1;
  MPI_Comm inter_node_comm = MPI_COMM_NULL;
  create_communicators(global_comm,
                       intra_node_comm,
                       inter_node_comm,
                       global_rank_id,
                       local_rank_id,
                       intercom_rank_id);

  // STEP 2: Exchange mesh metadata
  constexpr int NUM_NODES = 0;
  constexpr int NUM_FACES = 1;
  axom::IndexType mesh_metadata[2] = {0, 0};

  quest::STLReader reader;
  reader.setFileName(file);
  int rc = read_and_exchange_mesh_metadata(global_rank_id,
                                           global_comm,
                                           reader,
                                           mesh_metadata);
  if(rc != READ_SUCCESS)
  {
    return READ_FAILED;
  }

  // STEP 3: allocate shared buffer and wire pointers
  double* x = nullptr;
  double* y = nullptr;
  double* z = nullptr;
  axom::IndexType* conn = nullptr;
  MPI_Aint numBytes = allocate_shared_buffer(local_rank_id,
                                             intra_node_comm,
                                             mesh_metadata,
                                             x,
                                             y,
                                             z,
                                             conn,
                                             mesh_buffer,
                                             shared_window);
  SLIC_ASSERT(x != nullptr);
  SLIC_ASSERT(y != nullptr);
  SLIC_ASSERT(z != nullptr);
  SLIC_ASSERT(conn != nullptr);

  // STEP 5: allocate corresponding mesh object with external pointers.
  m = new TriangleMesh(mint::TRIANGLE,
                       mesh_metadata[NUM_FACES],
                       conn,
                       mesh_metadata[NUM_NODES],
                       x,
                       y,
                       z);

  // STEP 4: read in data to shared buffer
  if(global_rank_id == 0)
  {
    reader.getMesh(static_cast<TriangleMesh*>(m));
  }

  // STEP 5: inter-node communication
  if(intercom_rank_id >= 0)
  {
    MPI_Bcast(mesh_buffer, numBytes, MPI_UNSIGNED_CHAR, 0, inter_node_comm);
  }

  // STEP 6 free communicators

  MPI_Barrier(global_comm);
  mpi_comm_free(&inter_node_comm);
  return READ_SUCCESS;
}
#endif

/*
 * Reads in the surface mesh from the specified file.
 */
int read_mesh(const std::string& file, mint::Mesh*& m, MPI_Comm comm)
{
  // NOTE: STL meshes are always 3D
  constexpr int DIMENSION = 3;
  using TriangleMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

  // STEP 0: check input mesh pointer
  if(m != nullptr)
  {
    SLIC_WARNING("supplied mesh pointer is not null!");
    return READ_FAILED;
  }

  // STEP 1: allocate output mesh object
  m = new TriangleMesh(DIMENSION, mint::TRIANGLE);

  // STEP 2: allocate reader
  quest::STLReader* reader = nullptr;
#ifdef AXOM_USE_MPI
  reader = new quest::PSTLReader(comm);
#else
  AXOM_UNUSED_VAR(comm);
  reader = new quest::STLReader();
#endif

  // STEP 3: read the mesh from the STL file
  reader->setFileName(file);
  int rc = reader->read();
  if(rc == READ_SUCCESS)
  {
    reader->getMesh(static_cast<TriangleMesh*>(m));
  }
  else
  {
    SLIC_WARNING("reading STL file failed, setting mesh to NULL");
    delete m;
    m = nullptr;
  }

  // STEP 4: delete the reader
  delete reader;
  reader = nullptr;

  return rc;
}

/// Mesh Helper Methods

/*
 * Computes the bounds of the given mesh.
 */
void compute_mesh_bounds(const mint::Mesh* mesh, double* lo, double* hi)
{
  SLIC_ASSERT(mesh != nullptr);
  SLIC_ASSERT(lo != nullptr);
  SLIC_ASSERT(hi != nullptr);

  const int ndims = mesh->getDimension();

  // STEP 0: initialize lo,hi
  for(int i = 0; i < ndims; ++i)
  {
    lo[i] = std::numeric_limits<double>::max();
    hi[i] = std::numeric_limits<double>::lowest();
  }  // END for all dimensions

  // STEP 1: compute lo,hi
  double pt[3];
  const axom::IndexType numNodes = mesh->getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < numNodes; ++inode)
  {
    mesh->getNode(inode, pt);
    for(int i = 0; i < ndims; ++i)
    {
      lo[i] = (pt[i] < lo[i]) ? pt[i] : lo[i];
      hi[i] = (pt[i] > hi[i]) ? pt[i] : hi[i];
    }  // END for all dimensions

  }  // END for all nodes
}

/// Logger Initialize/Finalize Methods

/*
 * Helper method to initialize the Slic logger if needed.
 */
void logger_init(bool& isInitialized, bool& mustFinalize, bool verbose, MPI_Comm comm)
{
  if(isInitialized)
  {
    // Query has already initialized the logger
    return;
  }

  if(slic::isInitialized())
  {
    // logger is initialized by an application, the application will finalize
    isInitialized = true;
    mustFinalize = false;
    return;
  }

  // The SignedDistance Query must initialize the Slic logger and is then
  // also responsible for finalizing it when done
  isInitialized = true;
  mustFinalize = true;
  slic::initialize();

  slic::LogStream* ls = nullptr;
  std::string msgfmt = "[<LEVEL>]: <MESSAGE>\n";

#if defined(AXOM_USE_MPI) && defined(AXOM_USE_LUMBERJACK)
  constexpr int RLIMIT = 8;
  ls = new slic::LumberjackStream(&std::cout, comm, RLIMIT, msgfmt);
#elif defined(AXOM_USE_MPI) && !defined(AXOM_USE_LUMBERJACK)
  msgfmt.insert(0, "[<RANK>]", 8);
  ls = new slic::SynchronizedStream(&std::cout, comm, msgfmt);
#else
  AXOM_UNUSED_VAR(comm);
  ls = new slic::GenericOutputStream(&std::cout, msgfmt);
#endif

  slic::addStreamToAllMsgLevels(ls);
  slic::setLoggingMsgLevel((verbose) ? slic::message::Info
                                     : slic::message::Error);
}

/*
 * Finalizes the Slic logger (if needed)
 */
void logger_finalize(bool mustFinalize)
{
  if(mustFinalize)
  {
    slic::finalize();
  }
}

} /* end namespace internal */
} /* end namespace quest    */
} /* end namespace axom     */

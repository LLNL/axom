// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/quest/interface/signed_distance.hpp"

// Axom includes
#include "axom/quest/interface/internal/QuestHelpers.hpp"
#include "axom/quest/SignedDistance.hpp"

// Mint includes
#include "axom/mint/mesh/Mesh.hpp"

#include "axom/slic/interface/slic.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

namespace axom
{
namespace quest
{
//------------------------------------------------------------------------------
// INTERNAL DATA STRUCTURES
//------------------------------------------------------------------------------
namespace
{
constexpr int INIT_FAILED = -1;
constexpr int INIT_SUCCESS = 0;

using SignedDistance3D = SignedDistance<3>;
using SignedDistance2D = SignedDistance<2>;

/*!
 * \brief Holds the options for the SignedDistance query.
 */
static struct parameters_t
{
  int dimension; /*!< the dimension, 2 or 3 */

  int max_levels;         /*!< max levels of subdivision for the BVH */
  int max_occupancy;      /*!< max occupancy per BVH bin */
  bool verbose;           /*!< logger verbosity */
  bool is_closed_surface; /*!< indicates if the input is a closed surface */
  bool use_shared_memory; /*!< use MPI-3 shared memory for the surface mesh */
  bool compute_sign;      /*!< indicates if sign should be computed */

  /*!
   * \brief Default Constructor. Sets default values for the parameters.
   */
  parameters_t()
    : dimension(3)
    , max_levels(12)
    , max_occupancy(5)
    , verbose(false)
    , is_closed_surface(true)
    , use_shared_memory(false)
    , compute_sign(true)
  { }

} Parameters;

// TODO: note the SignedDistance query is currently only supported in 3-D
static SignedDistance3D* s_query = nullptr;
static mint::Mesh* s_surface_mesh = nullptr;
static bool s_must_delete_mesh = false;
static bool s_must_finalize_logger = false;
static bool s_logger_is_initialized = false;

#if defined(AXOM_USE_MPI) && defined(AXOM_USE_MPI3)
static unsigned char* s_shared_mesh_buffer = nullptr;
MPI_Comm s_intra_node_comm = MPI_COMM_NULL;
MPI_Win s_window = MPI_WIN_NULL;
#endif

}  // end anonymous namespace

//------------------------------------------------------------------------------
// SIGNED DISTANCE QUERY INTERFACE IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
int signed_distance_init(const std::string& file, MPI_Comm comm)
{
  // STEP 0: initialize logger
  internal::logger_init(s_logger_is_initialized,
                        s_must_finalize_logger,
                        Parameters.verbose,
                        comm);

  SLIC_ASSERT(s_query == nullptr);

  if(Parameters.dimension != 3)
  {
    SLIC_WARNING("the SignedDistance Query is currently only supported in 3D");
    return INIT_FAILED;
  }

  // STEP 0: read the STL mesh
  int rc = INIT_FAILED;

#if defined(AXOM_USE_MPI) && defined(AXOM_USE_MPI3)

  if(Parameters.use_shared_memory)
  {
    rc = internal::read_mesh_shared(file,
                                    comm,
                                    s_shared_mesh_buffer,
                                    s_surface_mesh,
                                    s_intra_node_comm,
                                    s_window);
  }
  else
  {
    rc = internal::read_mesh(file, s_surface_mesh, comm);
  }

#else

  SLIC_WARNING_IF(Parameters.use_shared_memory,
                  "Shared memory requires MPI3 and building Axom with "
                  "AXOM_USE_MPI3 set to ON");

  rc = internal::read_mesh(file, s_surface_mesh, comm);
#endif

  if(rc != 0)
  {
    SLIC_WARNING("reading mesh from [" << file << "] failed!");
    return INIT_FAILED;
  }

  // STEP 1: initialized the signed distance query
  s_must_delete_mesh = true;
  rc = signed_distance_init(s_surface_mesh, comm);
  return rc;
}

//------------------------------------------------------------------------------
int signed_distance_init(const mint::Mesh* m, MPI_Comm comm)
{
  internal::logger_init(s_logger_is_initialized,
                        s_must_finalize_logger,
                        Parameters.verbose,
                        comm);

  SLIC_ERROR_IF(signed_distance_initialized(),
                "signed distance query has already been initialized!");
  SLIC_ERROR_IF(m->getDimension() != 3,
                "signed distance query currently only support 3-D meshes");
  SLIC_ERROR_IF(
    m->getMeshType() != mint::UNSTRUCTURED_MESH,
    "signed distance query currently only supports unstructured meshes");
  SLIC_ERROR_IF(
    m->hasMixedCellTypes() == true,
    "signed distance query does not support meshes with mixed shape topology");
  SLIC_ERROR_IF(
    m->getCellType() != mint::TRIANGLE,
    "signed distance currently only support 3D triangular surface meshes");

  if(s_surface_mesh != m)
  {
    SLIC_ASSERT(s_surface_mesh == nullptr);
    s_surface_mesh = const_cast<mint::Mesh*>(m);
    s_must_delete_mesh = false;
  }

  s_query = new SignedDistance3D(s_surface_mesh,
                                 Parameters.is_closed_surface,
                                 Parameters.max_occupancy,
                                 Parameters.max_levels,
                                 Parameters.compute_sign);

  return INIT_SUCCESS;
}

//------------------------------------------------------------------------------
bool signed_distance_initialized() { return (s_query != nullptr); }

//------------------------------------------------------------------------------
void signed_distance_get_mesh_bounds(double* lo, double* hi)
{
  SLIC_ERROR_IF(!signed_distance_initialized(),
                "signed distance query must be initialized prior to"
                  << "calling get_mesh_bounds()");
  SLIC_ERROR_IF(lo == nullptr, "supplied buffer is null");
  SLIC_ERROR_IF(hi == nullptr, "supplied buffer is null");

  internal::compute_mesh_bounds(s_surface_mesh, lo, hi);
}

//------------------------------------------------------------------------------
void signed_distance_set_dimension(int dim)
{
  SLIC_ERROR_IF(dim != 3, "The signed distance query only support 3D");
  SLIC_ERROR_IF(
    signed_distance_initialized(),
    "signed distance query already initialized; setting option has no effect!");

  Parameters.dimension = dim;
}

//------------------------------------------------------------------------------
void signed_distance_set_closed_surface(bool status)
{
  SLIC_ERROR_IF(
    signed_distance_initialized(),
    "signed distance query already initialized; setting option has no effect!");

  Parameters.is_closed_surface = status;
}

//------------------------------------------------------------------------------
void signed_distance_set_compute_signs(bool computeSign)
{
  SLIC_ERROR_IF(
    signed_distance_initialized(),
    "signed distance query already initialized; setting option has no effect!");

  Parameters.compute_sign = computeSign;
}

//------------------------------------------------------------------------------
void signed_distance_set_max_levels(int maxLevels)
{
  SLIC_ERROR_IF(
    signed_distance_initialized(),
    "signed distance query already initialized; setting option has no effect!");

  Parameters.max_levels = maxLevels;
}

//------------------------------------------------------------------------------
void signed_distance_set_max_occupancy(int threshold)
{
  SLIC_ERROR_IF(
    signed_distance_initialized(),
    "signed distance query already initialized; setting option has no effect!");

  Parameters.max_occupancy = threshold;
}

//------------------------------------------------------------------------------
void signed_distance_set_verbose(bool status)
{
  SLIC_ERROR_IF(
    signed_distance_initialized(),
    "signed distance query already initialized; setting option has no effect!");

  Parameters.verbose = status;
}

//------------------------------------------------------------------------------
void signed_distance_use_shared_memory(bool status)
{
  SLIC_ERROR_IF(
    signed_distance_initialized(),
    "signed distance query already initialized; setting option has no effect!");

  Parameters.use_shared_memory = status;

#ifndef AXOM_USE_MPI3
  SLIC_WARNING_IF(Parameters.use_shared_memory,
                  "Enabling shared memory requires MPI-3. Option is ignored!");
#endif
}

//------------------------------------------------------------------------------
double signed_distance_evaluate(double x, double y, double z)
{
  SLIC_ERROR_IF(
    !signed_distance_initialized(),
    "signed distance query must be initialized prior to calling evaluate()!");

  double phi = s_query->computeDistance(x, y, z);
  return (phi);
}

//------------------------------------------------------------------------------
void signed_distance_evaluate(const double* x,
                              const double* y,
                              const double* z,
                              int npoints,
                              double* phi)
{
  SLIC_ERROR_IF(
    !signed_distance_initialized(),
    "signed distance query must be initialized prior to calling evaluate()!");
  SLIC_ERROR_IF(x == nullptr, "x-coords array is null");
  SLIC_ERROR_IF(y == nullptr, "y-coords array is null");
  SLIC_ERROR_IF(z == nullptr, "z-coords array is null");
  SLIC_ERROR_IF(phi == nullptr, "output phi array is null");

#ifdef AXOM_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for(int i = 0; i < npoints; ++i)
  {
    phi[i] = s_query->computeDistance(x[i], y[i], z[i]);
  }
}

//------------------------------------------------------------------------------
void signed_distance_finalize()
{
  if(s_query != nullptr)
  {
    delete s_query;
    s_query = nullptr;
  }

  if(s_surface_mesh != nullptr && s_must_delete_mesh)
  {
    delete s_surface_mesh;
  }

  s_surface_mesh = nullptr;

  SLIC_ASSERT(!signed_distance_initialized());
  internal::logger_finalize(s_must_finalize_logger);

#if defined(AXOM_USE_MPI) && defined(AXOM_USE_MPI3)
  internal::mpi_comm_free(&s_intra_node_comm);
  internal::mpi_win_free(&s_window);
  s_shared_mesh_buffer = nullptr;
#endif
}

}  // end namespace quest
}  // end namespace axom

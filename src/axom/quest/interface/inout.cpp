// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/interface/inout.hpp"

// Axom includes
#include "axom/slic.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"

#include "axom/quest/interface/internal/QuestHelpers.hpp"
#include "axom/quest/InOutOctree.hpp"

namespace axom
{
namespace quest
{
namespace internal
{
/// Parameter variables for the InOutQuery
struct InOutParameters
{
  bool m_verbose {false};
  int m_dimension {3};
  int m_segmentsPerPiece {100};  /// Used when linearizing curves
  double m_vertexWeldThreshold {1E-9};

  void setDefault()
  {
    m_verbose = false;
    m_dimension = 3;
    m_segmentsPerPiece = 100;
    m_vertexWeldThreshold = 1E-9;
  }
};

/*!
 * \struct InOutHelper
 * \brief A simple struct to hold the state of an InOut query
 *
 * This helper class is called by the quest::inout API functions.
 * \warning This class has minimal error checking since we assume the calling
 * functions have done the proper checks.
 */
template <int DIM>
struct InOutHelper
{
  using GeometricBoundingBox = primal::BoundingBox<double, DIM>;
  using SpacePt = primal::Point<double, DIM>;
  using SpaceVec = primal::Vector<double, DIM>;

  static_assert(DIM == 2 || DIM == 3, "InOutHelper only supports 2D and 3D");

  /// State variables for the InOutQuery
  struct State
  {
    bool m_initialized;
    bool m_logger_is_initialized;
    bool m_should_delete_logger;
    bool m_should_delete_mesh;
    slic::message::Level m_previousLevel;

    void setDefault()
    {
      m_initialized = false;
      m_logger_is_initialized = false;
      m_should_delete_logger = false;
      m_should_delete_mesh = false;
      m_previousLevel = slic::message::Num_Levels;
    }
  };

  InOutHelper() : m_surfaceMesh(nullptr), m_inoutTree(nullptr)
  {
    m_params.setDefault();
    m_state.setDefault();
  }

  ~InOutHelper() { finalize(); }

  /// Predicate to check if InOut query has been initialized
  bool isInitialized() const { return m_state.m_initialized; }

  /// Returns the spatial dimension of the mesh
  int getDimension() const { return DIM; }

  /// Sets the verbosity parameter
  void setVerbose(bool verbose) { m_params.m_verbose = verbose; }

  void setVertexWeldThreshold(double thresh)
  {
    m_params.m_vertexWeldThreshold = thresh;
  }

  void setSegmentsPerPiece(int numSegments)
  {
    m_params.m_segmentsPerPiece = numSegments;
  }

  /// Saves the current slic logging level
  void saveLoggingLevel()
  {
    if(slic::isInitialized())
    {
      m_state.m_previousLevel = slic::getLoggingMsgLevel();
    }
  }

  /// Restores the saved slic logging level
  void restoreLoggingLevel()
  {
    if(slic::isInitialized())
    {
      slic::setLoggingMsgLevel(m_state.m_previousLevel);
      slic::flushStreams();
    }
  }

  /*!
   * Initializes the InOut query from an stl file
   *
   * \sa inout_init
   */
  int initialize(const std::string& file, MPI_Comm comm)
  {
    mint::Mesh* mesh = nullptr;

    // load the mesh
    int rc = QUEST_INOUT_FAILED;

    switch(DIM)
    {
    case 2:
      rc = internal::read_c2c_mesh(file, m_params.m_segmentsPerPiece, mesh, comm);
      break;
    case 3:
      rc = internal::read_stl_mesh(file, mesh, comm);
      break;
    default:  // no-op
      break;
    }

    if(rc != QUEST_INOUT_SUCCESS)
    {
      SLIC_WARNING("reading mesh from [" << file << "] failed!");
      return QUEST_INOUT_FAILED;
    }
    m_state.m_should_delete_mesh = true;

    // the rest of the initialization is handled by the other function
    return initialize(mesh, comm);
  }

  /*!
   * Initializes the InOut query from a preloaded mesh
   *
   * \sa inout_init
   */
  int initialize(mint::Mesh*& mesh, MPI_Comm comm)
  {
    // initialize logger, if necessary
    internal::logger_init(m_state.m_logger_is_initialized,
                          m_state.m_should_delete_logger,
                          m_params.m_verbose,
                          comm);

    // Update log level based on verbosity
    this->saveLoggingLevel();
    slic::setLoggingMsgLevel(m_params.m_verbose ? slic::message::Debug
                                                : slic::message::Warning);

    // handle mesh pointer, with some error checking
    if(mesh == nullptr)
    {
      SLIC_WARNING("Cannot initialize: mesh was NULL");
      return QUEST_INOUT_FAILED;
    }
    m_surfaceMesh = mesh;

    if(m_surfaceMesh->getDimension() != m_params.m_dimension)
    {
      SLIC_WARNING("Incorrect dimensionality for mesh."
                   << "Expected " << m_params.m_dimension << ", "
                   << "but got " << m_surfaceMesh->getDimension());
      return QUEST_INOUT_FAILED;
    }

    // compute the mesh bounding box and center of mass
    m_meshBoundingBox.clear();
    m_meshCenterOfMass = SpacePt::zero();
    SpacePt pt;

    const int numMeshNodes = m_surfaceMesh->getNumberOfNodes();
    if(numMeshNodes > 0)
    {
      for(int i = 0; i < numMeshNodes; ++i)
      {
        m_surfaceMesh->getNode(i, pt.data());

        m_meshBoundingBox.addPoint(pt);
        m_meshCenterOfMass.array() += pt.array();
      }

      m_meshCenterOfMass.array() /= numMeshNodes;
      SLIC_ASSERT(m_meshBoundingBox.isValid());
    }

    // initialize InOutOctree
    m_inoutTree = new InOutOctree<DIM>(m_meshBoundingBox, m_surfaceMesh);

    // set params
    m_inoutTree->setVertexWeldThreshold(m_params.m_vertexWeldThreshold);

    // initialize the spatial index
    m_inoutTree->generateIndex();

    // Update the mesh parameter since the InOutOctree modifies the mesh
    mesh = m_surfaceMesh;

    this->restoreLoggingLevel();

    // set the initialized flag to true
    m_state.m_initialized = true;

    return QUEST_INOUT_SUCCESS;
  }

  /*!
   * Finalizes the InOut query
   *
   * \sa inout_finalize
   */
  int finalize()
  {
    // deal with spatial index
    if(m_inoutTree != nullptr)
    {
      delete(m_inoutTree);
      m_inoutTree = nullptr;
    }

    // deal with mesh
    if(m_state.m_should_delete_mesh)
    {
      delete(m_surfaceMesh);
    }
    m_surfaceMesh = nullptr;

    // deal with logging
    internal::logger_finalize(m_state.m_should_delete_logger);

    m_state.m_initialized = false;

    // Clear state and input parameters
    m_state.setDefault();
    m_params.setDefault();

    return QUEST_INOUT_SUCCESS;
  }

  /*!
   * Returns the precomputed mesh bounding box
   */
  const GeometricBoundingBox& getBoundingBox() const
  {
    return m_meshBoundingBox;
  }

  /*!
   * Returns the precomputed mesh center of mass
   */
  const SpacePt& getCenterOfMass() const { return m_meshCenterOfMass; }

  /// Predicate to determine if a point is inside the surface
  bool within(double x, double y, double z = 0.) const
  {
    return m_inoutTree->within(SpacePt {x, y, z});
  }

  /*!
   * Batched containment query over an array of points
   *
   * \sa inout_evaluate
   */
  int within(const double* x,
             const double* y,
             const double* z,
             int npoints,
             int* res) const
  {
    if(z == nullptr)
    {
#ifdef AXOM_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
      for(int i = 0; i < npoints; ++i)
      {
        const bool ins = m_inoutTree->within(SpacePt {x[i], y[i]});
        res[i] = ins ? 1 : 0;
      }
    }
    else
    {
#ifdef AXOM_USE_OPENMP
  #pragma omp parallel for schedule(static)
#endif
      for(int i = 0; i < npoints; ++i)
      {
        const bool ins = m_inoutTree->within(SpacePt {x[i], y[i], z[i]});
        res[i] = ins ? 1 : 0;
      }
    }

    return QUEST_INOUT_SUCCESS;
  }

private:
  mint::Mesh* m_surfaceMesh;
  InOutOctree<DIM>* m_inoutTree;
  GeometricBoundingBox m_meshBoundingBox;
  SpacePt m_meshCenterOfMass;

  InOutParameters m_params;
  State m_state;
};
}  // end namespace internal

/// Static instance of the InOutHelper in 3D.
/// Used by the quest::inout API functions
static internal::InOutHelper<3> s_inoutHelper3D;

/// Static instance of the InOutHelper in 2D
/// Used by the quest::inout API functions
static internal::InOutHelper<2> s_inoutHelper2D;

static internal::InOutParameters s_inoutParams;

//------------------------------------------------------------------------------

bool inout_initialized()
{
  const int dim = inout_get_dimension();
  return (dim == 2) ? s_inoutHelper2D.isInitialized()
                    : s_inoutHelper3D.isInitialized();
}

int inout_init(const std::string& file, MPI_Comm comm)
{
  const int dim = inout_get_dimension();
  int rc = QUEST_INOUT_FAILED;

  if(inout_initialized())
  {
    SLIC_WARNING("quest inout query already initialized ");
    return rc;
  }

  switch(dim)
  {
  case 2:
    s_inoutHelper2D.setVerbose(s_inoutParams.m_verbose);
    s_inoutHelper2D.setSegmentsPerPiece(s_inoutParams.m_segmentsPerPiece);
    s_inoutHelper2D.setVertexWeldThreshold(s_inoutParams.m_vertexWeldThreshold);

    rc = s_inoutHelper2D.initialize(file, comm);
    if(rc == QUEST_INOUT_FAILED)
    {
      s_inoutHelper2D.restoreLoggingLevel();
    }
    break;

  case 3:
    s_inoutHelper3D.setVerbose(s_inoutParams.m_verbose);
    s_inoutHelper3D.setVertexWeldThreshold(s_inoutParams.m_vertexWeldThreshold);

    rc = s_inoutHelper3D.initialize(file, comm);
    if(rc == QUEST_INOUT_FAILED)
    {
      s_inoutHelper3D.restoreLoggingLevel();
    }
    break;

  default:
    rc = QUEST_INOUT_FAILED;
    break;
  }

  return rc;
}

int inout_init(mint::Mesh*& mesh, MPI_Comm comm)
{
  const int dim = inout_get_dimension();
  int rc = QUEST_INOUT_FAILED;

  if(inout_initialized())
  {
    SLIC_WARNING("quest inout query already initialized ");
    return QUEST_INOUT_FAILED;
  }

  switch(dim)
  {
  case 2:
    s_inoutHelper2D.setVerbose(s_inoutParams.m_verbose);
    s_inoutHelper2D.setSegmentsPerPiece(s_inoutParams.m_segmentsPerPiece);
    s_inoutHelper2D.setVertexWeldThreshold(s_inoutParams.m_vertexWeldThreshold);

    rc = s_inoutHelper2D.initialize(mesh, comm);
    if(rc == QUEST_INOUT_FAILED)
    {
      s_inoutHelper2D.restoreLoggingLevel();
    }

  case 3:
    s_inoutHelper3D.setVerbose(s_inoutParams.m_verbose);
    s_inoutHelper3D.setSegmentsPerPiece(s_inoutParams.m_segmentsPerPiece);
    s_inoutHelper3D.setVertexWeldThreshold(s_inoutParams.m_vertexWeldThreshold);

    rc = s_inoutHelper3D.initialize(mesh, comm);
    if(rc == QUEST_INOUT_FAILED)
    {
      s_inoutHelper3D.restoreLoggingLevel();
    }
    break;

  default:
    rc = QUEST_INOUT_FAILED;
    break;
  }

  return rc;
}

int inout_finalize()
{
  const int dim = inout_get_dimension();

  // Finalize the 2D and 3D structures and reset the params
  int rc2 = s_inoutHelper2D.finalize();
  int rc3 = s_inoutHelper3D.finalize();
  s_inoutParams.setDefault();

  return (dim == 2) ? rc2 : rc3;
}

//------------------------------------------------------------------------------

int inout_mesh_min_bounds(double* coords)
{
  if(!inout_initialized())
  {
    SLIC_WARNING("quest inout query must be initialized "
                 << "prior to calling quest inout interface functions");

    return QUEST_INOUT_FAILED;
  }

  SLIC_ERROR_IF(coords == nullptr, "supplied buffer 'coords' is null");

  const int dim = inout_get_dimension();
  int rc = QUEST_INOUT_SUCCESS;
  switch(dim)
  {
  case 2:
    s_inoutHelper2D.getBoundingBox().getMin().array().to_array(coords);
    break;
  case 3:
    s_inoutHelper3D.getBoundingBox().getMin().array().to_array(coords);
    break;
  default:
    rc = QUEST_INOUT_FAILED;
    break;
  }

  return rc;
}

int inout_mesh_max_bounds(double* coords)
{
  if(!inout_initialized())
  {
    SLIC_WARNING("quest inout query must be initialized "
                 << "prior to calling quest inout interface functions");

    return QUEST_INOUT_FAILED;
  }

  SLIC_ERROR_IF(coords == nullptr, "supplied buffer 'coords' is null");

  const int dim = inout_get_dimension();
  int rc = QUEST_INOUT_SUCCESS;
  switch(dim)
  {
  case 2:
    s_inoutHelper2D.getBoundingBox().getMax().array().to_array(coords);
    break;
  case 3:
    s_inoutHelper3D.getBoundingBox().getMax().array().to_array(coords);
    break;
  default:
    rc = QUEST_INOUT_FAILED;
    break;
  }

  return rc;
}

int inout_mesh_center_of_mass(double* coords)
{
  if(!inout_initialized())
  {
    SLIC_WARNING("quest inout query must be initialized "
                 << "prior to calling quest inout interface functions");

    return QUEST_INOUT_FAILED;
  }

  SLIC_ERROR_IF(coords == nullptr, "supplied buffer 'coords' is null");

  const int dim = inout_get_dimension();
  int rc = QUEST_INOUT_SUCCESS;
  switch(dim)
  {
  case 2:
    s_inoutHelper2D.getCenterOfMass().array().to_array(coords);
    break;
  case 3:
    s_inoutHelper3D.getCenterOfMass().array().to_array(coords);
    break;
  default:
    rc = QUEST_INOUT_FAILED;
    break;
  }

  return rc;
}

//------------------------------------------------------------------------------

bool inout_evaluate(double x, double y)
{
  if(!inout_initialized())
  {
    SLIC_WARNING("quest inout query must be initialized "
                 << "prior to calling quest inout interface functions");

    return false;
  }

  const int dim = inout_get_dimension();
  return (dim == 2) ? s_inoutHelper2D.within(x, y) : s_inoutHelper3D.within(x, y);
}

bool inout_evaluate(double x, double y, double z)
{
  if(!inout_initialized())
  {
    SLIC_WARNING("quest inout query must be initialized "
                 << "prior to calling quest inout interface functions");

    return false;
  }

  const int dim = inout_get_dimension();
  return (dim == 2) ? s_inoutHelper2D.within(x, y, z)
                    : s_inoutHelper3D.within(x, y, z);
}

int inout_evaluate(const double* x, const double* y, int npoints, int* res)
{
  if(!inout_initialized())
  {
    SLIC_WARNING("quest inout query must be initialized "
                 << "prior to calling quest inout interface functions");

    return QUEST_INOUT_FAILED;
  }

  const int dim = inout_get_dimension();
  if(x == nullptr || y == nullptr || res == nullptr)
  {
    SLIC_WARNING("supplied buffers must NOT be null");
    return QUEST_INOUT_FAILED;
  }

  return (dim == 2) ? s_inoutHelper2D.within(x, y, nullptr, npoints, res)
                    : s_inoutHelper3D.within(x, y, nullptr, npoints, res);
}

int inout_evaluate(const double* x,
                   const double* y,
                   const double* z,
                   int npoints,
                   int* res)
{
  if(!inout_initialized())
  {
    SLIC_WARNING("quest inout query must be initialized "
                 << "prior to calling quest inout interface functions");

    return QUEST_INOUT_FAILED;
  }

  const int dim = inout_get_dimension();
  if(x == nullptr || y == nullptr || (dim == 3 && z == nullptr) || res == nullptr)
  {
    SLIC_WARNING("supplied buffers must NOT be null");
    return QUEST_INOUT_FAILED;
  }

  return (dim == 2) ? s_inoutHelper2D.within(x, y, z, npoints, res)
                    : s_inoutHelper3D.within(x, y, z, npoints, res);
}

//------------------------------------------------------------------------------

int inout_get_dimension() { return s_inoutParams.m_dimension; }

//------------------------------------------------------------------------------

int inout_set_verbose(bool verbosity)
{
  if(inout_initialized())
  {
    SLIC_WARNING("quest inout query must NOT be initialized "
                 << "prior to calling 'inout_set_verbose'");

    return QUEST_INOUT_FAILED;
  }

  s_inoutParams.m_verbose = verbosity;

  return QUEST_INOUT_SUCCESS;
}

int inout_set_dimension(int dim)
{
  if(inout_initialized())
  {
    SLIC_WARNING("quest inout query must NOT be initialized "
                 << "prior to calling 'inout_set_dimension'");

    return QUEST_INOUT_FAILED;
  }

  if(!(dim == 2 || dim == 3))
  {
    SLIC_WARNING("quest inout query only supports 2D or 3D queries."
                 << " Supplied dimension was " << dim);

    return QUEST_INOUT_FAILED;
  }

  s_inoutParams.m_dimension = dim;

  return QUEST_INOUT_SUCCESS;
}

int inout_set_segments_per_piece(int num_segments)
{
  if(inout_initialized())
  {
    SLIC_WARNING("quest inout query must NOT be initialized "
                 << "prior to calling 'inout_set_segments_per_piece'");

    return QUEST_INOUT_FAILED;
  }

  if(num_segments < 1)
  {
    SLIC_WARNING("quest inout query: segments per piece must be at least 1."
                 << " Supplied value was " << num_segments);

    return QUEST_INOUT_FAILED;
  }

  s_inoutParams.m_segmentsPerPiece = num_segments;

  return QUEST_INOUT_SUCCESS;
}

int inout_set_vertex_weld_threshold(double thresh)
{
  if(inout_initialized())
  {
    SLIC_WARNING("quest inout query must NOT be initialized "
                 << "prior to calling 'inout_set_vertex_weld_threshold'");

    return QUEST_INOUT_FAILED;
  }

  if(thresh < 0.)
  {
    SLIC_WARNING("quest inout query: vertex weld threshold must be postive."
                 << " Supplied value was " << thresh);

    return QUEST_INOUT_FAILED;
  }

  s_inoutParams.m_vertexWeldThreshold = thresh;

  return QUEST_INOUT_SUCCESS;
}

}  // end namespace quest
}  // end namespace axom

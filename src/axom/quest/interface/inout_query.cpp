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
#include "axom/quest/interface/inout_query.hpp"

// Axom includes
#include "axom/slic.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"

#include "axom/quest/interface/internal/QuestHelpers.hpp"
#include "axom/quest/geom/InOutOctree.hpp"


#include <algorithm>  // for std::copy

namespace axom
{
namespace quest
{

namespace internal
{
/*!
 * \struct
 * \brief A simple struct to hold the state of an InOut query
 */
template<int DIM>
struct InOutHelper
{
  typedef primal::BoundingBox< double, DIM> GeometricBoundingBox;
  typedef primal::Point< double, DIM> SpacePt;
  typedef primal::Vector< double, DIM> SpaceVec;

  struct Parameters
  {
    bool m_verbose;
    int m_dimension;

    void setDefault()
    {
      m_verbose = false;
      m_dimension = 3;
    }
  };

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

  InOutHelper() :
    m_surfaceMesh(nullptr),
    m_inoutTree(nullptr)
  {
    m_params.setDefault();
    m_state.setDefault();
  }

  ~InOutHelper()
  {
    finalize();
  }

  bool isInitialized() const { return m_state.m_initialized; }

  void setVerbose(bool verbose)
  {
    m_params.m_verbose=verbose;
  }


  void saveLoggingLevel()
  {
    if(slic::isInitialized())
    {
      m_state.m_previousLevel = slic::getLoggingMsgLevel();
    }
  }

  void restoreLoggingLevel()
  {
    if(slic::isInitialized())
    {
      slic::setLoggingMsgLevel(m_state.m_previousLevel);
    }
  }

  int initialize(const std::string& file, MPI_Comm comm)
  {
    mint::Mesh* mesh = nullptr;

    // load the mesh
    int rc = QUEST_INOUT_FAILED;
    rc = internal::read_mesh(file, mesh, comm);

    if (rc != QUEST_INOUT_SUCCESS)
    {
      SLIC_WARNING("reading mesh from [" << file << "] failed!");
      return QUEST_INOUT_FAILED;
    }
    m_state.m_should_delete_mesh = true;

    // the rest of the initialization is handled by the other function
    return initialize(mesh, comm);
  }

  int initialize(mint::Mesh*& mesh, MPI_Comm comm)
  {
    // initialize logger, if necessary
    internal::logger_init(
      m_state.m_logger_is_initialized,
      m_state.m_should_delete_logger,
      m_params.m_verbose,
      comm
      );

    // Update log level based on verbosity
    this->saveLoggingLevel();
    slic::setLoggingMsgLevel(m_params.m_verbose
                             ? slic::message::Debug
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
                   << "Expected " << m_params.m_dimension <<", "
                   << "but got " << m_surfaceMesh->getDimension() );
      return QUEST_INOUT_FAILED;
    }

    // compute the mesh bounding box and center of mass
    m_meshBoundingBox.clear();
    m_meshCenterOfMass = SpacePt::zero();
    SpacePt pt;

    const int numMeshNodes = m_surfaceMesh->getNumberOfNodes();
    if(numMeshNodes > 0)
    {
      for ( int i=0 ; i < numMeshNodes ; ++i )
      {
        m_surfaceMesh->getNode( i, pt.data() );

        m_meshBoundingBox.addPoint( pt );
        m_meshCenterOfMass.array() += pt.array();
      }

      m_meshCenterOfMass.array() /= numMeshNodes;
      SLIC_ASSERT( m_meshBoundingBox.isValid() );
    }

    // initialize InOutOctree
    m_inoutTree = new InOutOctree<DIM>(m_meshBoundingBox, m_surfaceMesh);
    m_inoutTree->generateIndex();

    // Update the mesh parameter since the InOutOctree modifies the mesh
    mesh = m_surfaceMesh;

    this->restoreLoggingLevel();

    // set the initialized flag to true
    m_state.m_initialized = true;

    return QUEST_INOUT_SUCCESS;
  }

  int finalize()
  {
    // deal with spatial index
    if (m_inoutTree != nullptr)
    {
      delete(m_inoutTree);
      m_inoutTree = nullptr;
    }

    // deal with mesh
    if (m_state.m_should_delete_mesh)
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

  const GeometricBoundingBox& getBoundingBox() const
  {
    return m_meshBoundingBox;
  }

  const SpacePt& getCenterOfMass() const
  {
    return m_meshCenterOfMass;
  }

  bool within(double x, double y, double z) const
  {
    return m_inoutTree->within(SpacePt::make_point(x,y,z));
  }

private:
  mint::Mesh* m_surfaceMesh;
  InOutOctree<DIM>* m_inoutTree;
  GeometricBoundingBox m_meshBoundingBox;
  SpacePt m_meshCenterOfMass;
  Parameters m_params;
  State m_state;

};
} // end namespace internal

// Static instance
static internal::InOutHelper<3> s_inoutHelper;

bool inout_initialized()
{
  return s_inoutHelper.isInitialized();
}


int inout_init(const std::string& file, MPI_Comm comm)
{
  if(inout_initialized())
  {
    SLIC_WARNING( "quest inout query already initialized ");
    return QUEST_INOUT_FAILED;
  }

  int rc = s_inoutHelper.initialize(file, comm);
  if(rc == QUEST_INOUT_FAILED)
  {
    s_inoutHelper.restoreLoggingLevel();
  }
  return rc;
}

int inout_init(mint::Mesh*& mesh, MPI_Comm comm)
{
  if(inout_initialized())
  {
    SLIC_WARNING( "quest inout query already initialized ");
    return QUEST_INOUT_FAILED;
  }

  int rc = s_inoutHelper.initialize(mesh, comm);
  if(rc == QUEST_INOUT_FAILED)
  {
    s_inoutHelper.restoreLoggingLevel();
  }
  return rc;
}

int inout_finalize()
{
  return s_inoutHelper.finalize();
}


int inout_set_verbose(bool verbosity)
{
  if(inout_initialized())
  {
    SLIC_WARNING( "quest inout query must NOT be initialized "
                  << "prior to calling 'inout_set_verbose'");

    return QUEST_INOUT_FAILED;
  }

  s_inoutHelper.setVerbose(verbosity);

  return QUEST_INOUT_SUCCESS;
}

int inout_mesh_min_bounds(double* coords)
{
  if(!inout_initialized())
  {
    SLIC_WARNING( "quest inout query must be initialized "
                  << "prior to calling quest inout interface functions");

    return QUEST_INOUT_FAILED;
  }

  SLIC_ERROR_IF(coords==nullptr, "supplied buffer 'coords' is null");

  auto& bbox = s_inoutHelper.getBoundingBox();
  std::copy_n(bbox.getMin().data(), bbox.dimension(), coords);

  return QUEST_INOUT_SUCCESS;
}

int inout_mesh_max_bounds(double* coords)
{
  if(!inout_initialized())
  {
    SLIC_WARNING( "quest inout query must be initialized "
                  << "prior to calling quest inout interface functions");

    return QUEST_INOUT_FAILED;
  }

  SLIC_ERROR_IF(coords==nullptr, "supplied buffer 'coords' is null");

  auto& bbox = s_inoutHelper.getBoundingBox();
  std::copy_n(bbox.getMax().data(), bbox.dimension(), coords);

  return QUEST_INOUT_SUCCESS;
}

int inout_mesh_center_of_mass(double* coords)
{
  if(!inout_initialized())
  {
    SLIC_WARNING( "quest inout query must be initialized "
                  << "prior to calling quest inout interface functions");

    return QUEST_INOUT_FAILED;
  }

  SLIC_ERROR_IF(coords==nullptr, "supplied buffer 'coords' is null");

  auto& cmass = s_inoutHelper.getCenterOfMass();
  std::copy_n(cmass.data(), cmass.dimension(), coords);

  return QUEST_INOUT_SUCCESS;
}


bool inout_inside(double x, double y, double z)
{
  if(!inout_initialized())
  {
    SLIC_WARNING( "quest inout query must be initialized "
                  << "prior to calling 'inout_get_mesh_bounds'");

    return false;
  }

  return s_inoutHelper.within(x,y,z);
}

} // end namespace quest
} // end namespace axom

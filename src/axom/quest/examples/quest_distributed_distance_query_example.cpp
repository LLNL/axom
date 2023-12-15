// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_distributed_distance_query_example.cpp
 * \brief Driver for a distributed distance query
 */

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/sidre.hpp"
#include "axom/quest.hpp"
#include "axom/slam.hpp"
#include "axom/core/Types.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_blueprint_mpi.hpp"
#include "conduit_relay_io_blueprint.hpp"
#include "conduit_relay_mpi_io_blueprint.hpp"

#include "axom/quest/DistributedClosestPoint.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#ifndef AXOM_USE_MPI
  #error This example requires Axom to be configured with MPI
#endif
#include "mpi.h"

// C/C++ includes
#include <string>
#include <limits>
#include <map>
#include <vector>
#include <cmath>

namespace quest = axom::quest;
namespace slic = axom::slic;
namespace sidre = axom::sidre;
namespace slam = axom::slam;
namespace spin = axom::spin;
namespace primal = axom::primal;
namespace mint = axom::mint;
namespace numerics = axom::numerics;

using RuntimePolicy = axom::runtime_policy::Policy;

// converts the input string into an 80 character string
// padded on both sides with '=' symbols
std::string banner(const std::string& str)
{
  return axom::fmt::format("{:=^80}", str);
}

/// Struct to parse and store the input parameters
struct Input
{
public:
  std::string meshFile;
  std::string distanceFile {"cp_coords"};
  std::string objectFile {"object_mesh"};

  double circleRadius {1.0};
  std::vector<double> circleCenter {0.0, 0.0};
  // TODO: Ensure that circleCenter size matches dimensionality.
  int circlePoints {100};
  RuntimePolicy policy {RuntimePolicy::seq};

  double distThreshold {std::numeric_limits<double>::max()};

  bool checkResults {false};

  bool randomSpacing {true};

  std::vector<unsigned int> objDomainCountRange {1, 1};

private:
  bool m_verboseOutput {false};
  double m_emptyRankProbability {0.};

public:
  bool isVerbose() const { return m_verboseOutput; }
  double percentEmptyRanks() const { return m_emptyRankProbability; }

  std::string getDCMeshName() const
  {
    using axom::utilities::string::removeSuffix;

    // Remove the parent directories and file suffix
    std::string name = axom::Path(meshFile).baseName();
    name = removeSuffix(name, ".root");

    return name;
  }

  std::string getMdMeshName() const
  {
    using axom::utilities::string::removeSuffix;

    // Remove the parent directories and file suffix
    std::string name = axom::Path(meshFile).baseName();
    name = removeSuffix(name, ".root");

    return name;
  }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-m,--mesh-file", meshFile)
      ->description(
        "Path to multidomain computational mesh following conduit blueprint "
        "convention.")
      ->check(axom::CLI::ExistingFile);

    app.add_option("-s,--distance-file", distanceFile)
      ->description("Name of output mesh file containing closest distance.")
      ->capture_default_str();

    app.add_option("-o,--object-file", objectFile)
      ->description("Name of output file containing object mesh.")
      ->capture_default_str();

    app.add_flag("-v,--verbose,!--no-verbose", m_verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_option("--empty-rank-probability", m_emptyRankProbability)
      ->description(
        "Probability that a rank's data is empty "
        "(tests code's ability to handle empty ranks)")
      ->check(axom::CLI::Range(0., 1.))
      ->capture_default_str();

    app.add_option("-r,--radius", circleRadius)
      ->description("Radius for circle")
      ->capture_default_str();

    auto* circle_options =
      app.add_option_group("circle",
                           "Options for setting up the circle of points");
    circle_options->add_option("--center", circleCenter)
      ->description("Center for object (x,y[,z])")
      ->expected(2, 3);

    circle_options->add_option("--obj-domain-count-range", objDomainCountRange)
      ->description("Range of object domain counts/rank (min, max)")
      ->expected(2);

    app.add_flag("--random-spacing,!--no-random-spacing", randomSpacing)
      ->description("Enable/disable random spacing of circle points")
      ->capture_default_str();

    app.add_option("-d,--dist-threshold", distThreshold)
      ->check(axom::CLI::NonNegativeNumber)
      ->description("Distance threshold to search")
      ->capture_default_str();

    app.add_option("-n,--num-samples", circlePoints)
      ->description("Number of points for circle")
      ->capture_default_str();

    app.add_option("-p, --policy", policy)
      ->description("Set runtime policy for point query method")
      ->capture_default_str()
      ->transform(
        axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

    app.add_flag("-c,--check-results,!--no-check-results", checkResults)
      ->description(
        "Enable/disable checking results against analytical solution")
      ->capture_default_str();

    app.get_formatter()->column_width(60);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(m_verboseOutput ? slic::message::Debug
                                             : slic::message::Info);
  }
};

/**
 *  \brief Simple wrapper to a blueprint particle mesh
 *
 *  Given a sidre Group, creates the stubs for a mesh blueptint particle mesh
 */
struct BlueprintParticleMesh
{
public:
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;
  using PointArray2D = axom::Array<Point2D>;
  using PointArray3D = axom::Array<Point3D>;

  explicit BlueprintParticleMesh(sidre::Group* group,
                                 const std::string& topology,
                                 const std::string& coordset)
    : m_topologyName(topology)
    , m_coordsetName(coordset)
    , m_group(group)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_nranks);
  }

  explicit BlueprintParticleMesh(sidre::Group* group)
    : m_topologyName()
    , m_coordsetName()
    , m_group(group)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_nranks);
  }

  /// Gets the root group for this mesh blueprint
  sidre::Group* root_group() const { return m_group; }

  /// Gets number of domains in the multidomain particle mesh
  axom::IndexType domain_count() const { return m_group->getNumGroups(); }

  /// Gets a domain group.
  sidre::Group* domain_group(axom::IndexType groupIdx) const
  {
    SLIC_ASSERT(groupIdx < m_group->getNumGroups());
    return m_group->getGroup(groupIdx);
  }
  /// Gets the parent group for the blueprint coordinate set
  sidre::Group* coords_group(axom::IndexType groupIdx) const
  {
    return domain_group(groupIdx)->getGroup("coordsets")->getGroup(m_coordsetName);
  }
  /// Gets the parent group for the blueprint mesh topology
  sidre::Group* topo_group(axom::IndexType groupIdx) const
  {
    return domain_group(groupIdx)->getGroup("topologies")->getGroup(m_topologyName);
  }
  /// Gets the parent group for the blueprint fields
  sidre::Group* fields_group(axom::IndexType groupIdx) const
  {
    return domain_group(groupIdx)->getGroup("fields");
  }

  const std::string& getTopologyName() const { return m_topologyName; }
  const std::string& getCoordsetName() const { return m_coordsetName; }

  /// Gets the MPI rank for this mesh
  int getRank() const { return m_rank; }
  /// Gets the number of ranks in the problem
  int getNumRanks() const { return m_nranks; }

  /*!
    @brief Returns the number of points in a particle mesh domain
    including ghost points.
  */
  int numPoints(axom::IndexType dIdx) const
  {
    int rval = 0;
    auto* cg = coords_group(dIdx);
    SLIC_ASSERT(cg != nullptr && cg->hasView("values/x"));
    rval = cg->getView("values/x")->getNumElements();
    return rval;
  }
  /// Returns the number of points in the particle mesh
  int numPoints() const
  {
    int rval = 0;
    const axom::IndexType domCount = domain_count();
    for(axom::IndexType dIdx = 0; dIdx < domCount; ++dIdx)
    {
      rval += numPoints(dIdx);
    }
    return rval;
  }

  int dimension() const { return m_dimension; }

  /*!
    @brief Read a blueprint mesh and store it internally in m_group.

    If the topology wasn't specified in the constructor, the first
    topology from the file is used.  The coordset name will be
    replaced with the one corresponding to the topology.
  */
  void read_blueprint_mesh(const std::string& meshFilename)
  {
    SLIC_ASSERT(!meshFilename.empty());

    conduit::Node mdMesh;
    conduit::relay::mpi::io::blueprint::load_mesh(meshFilename,
                                                  mdMesh,
                                                  MPI_COMM_WORLD);
    assert(conduit::blueprint::mesh::is_multi_domain(mdMesh));
    conduit::index_t domCount =
      conduit::blueprint::mesh::number_of_domains(mdMesh);

    if(domCount > 0)
    {
      m_coordsAreStrided = mdMesh[0]
                             .fetch_existing("topologies/mesh/elements/dims")
                             .has_child("strides");
      if(m_coordsAreStrided)
      {
        SLIC_WARNING(axom::fmt::format(
          "Mesh '{}' is strided.  Stride support is under development.",
          meshFilename));
      }
    }

    if(domCount > 0)
    {
      if(m_topologyName.empty())
      {
        // No topology given.  Pick the first one.
        m_topologyName = mdMesh[0].fetch_existing("topologies")[0].name();
      }
      auto topologyPath = axom::fmt::format("topologies/{}", m_topologyName);

      m_coordsetName =
        mdMesh[0].fetch_existing(topologyPath + "/coordset").as_string();
      const conduit::Node coordsetNode =
        mdMesh[0].fetch_existing("coordsets").fetch_existing(m_coordsetName);
      m_dimension = conduit::blueprint::mesh::coordset::dims(coordsetNode);
    }

    MPI_Allreduce(MPI_IN_PLACE, &m_dimension, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    SLIC_ASSERT(m_dimension > 0);

    if(domCount > 0)
    {
      // Put mdMesh into sidre Group.
      const bool goodImport = m_group->importConduitTree(mdMesh, false);
      SLIC_ASSERT(goodImport);
      SLIC_ASSERT(m_group->getNumGroups() == domCount);
      AXOM_UNUSED_VAR(goodImport);
    }

    bool valid = isValid();
    SLIC_ASSERT(valid);
    AXOM_UNUSED_VAR(valid);
  }

  /*!  @brief Set the coordinate data from an array of primal Points

    The points are assigned to a new domain (in the multidomain context).

    This method is for manually creating the mesh.  Don't use it if
    the mesh is read in.
  */
  template <int NDIMS>
  void setPoints(axom::IndexType domainId,
                 const axom::Array<primal::Point<double, NDIMS>>& pts)
  {
    axom::IndexType localIdx = createBlueprintStubs();
    SLIC_ASSERT(domain_group(localIdx) != nullptr);
    domain_group(localIdx)->createViewScalar<std::int64_t>("state/domain_id",
                                                           domainId);

    const int SZ = pts.size();

    if(m_dimension == -1)
    {
      m_dimension = NDIMS;
    }
    else
    {
      SLIC_ASSERT(NDIMS == m_dimension);
    }

    // lambda to create a strided view into the buffer
    // uses workaround for empty meshes since apply() requires size > 0
    auto createAndApplyView = [=](sidre::Group* grp,
                                  const std::string& path,
                                  sidre::Buffer* buf,
                                  int dim,
                                  int sz) {
      if(sz > 0)
      {
        grp->createView(path)->attachBuffer(buf)->apply(sz, dim, NDIMS);
      }
      else
      {
        grp->createViewAndAllocate(path, sidre::DOUBLE_ID, 0);
      }
    };

    // create views into a shared buffer for the coordinates, with stride NDIMS
    {
      auto* buf = domain_group(localIdx)
                    ->getDataStore()
                    ->createBuffer(sidre::DOUBLE_ID, NDIMS * SZ)
                    ->allocate();

      createAndApplyView(coords_group(localIdx), "values/x", buf, 0, SZ);
      if(NDIMS > 1)
      {
        createAndApplyView(coords_group(localIdx), "values/y", buf, 1, SZ);
      }
      if(NDIMS > 2)
      {
        createAndApplyView(coords_group(localIdx), "values/z", buf, 2, SZ);
      }

      // copy coordinate data into the buffer
      const std::size_t nbytes = sizeof(double) * SZ * NDIMS;
      axom::copy(buf->getVoidPtr(), pts.data(), nbytes);
    }

    // set the default connectivity
    // May be required by an old version of visit.  May not be needed by newer versions of visit.
    sidre::Array<int> arr(
      topo_group(localIdx)->createView("elements/connectivity"),
      SZ,
      SZ);
    for(int i = 0; i < SZ; ++i)
    {
      arr[i] = i;
    }
  }

  template <int NDIMS>
  axom::Array<primal::Point<double, NDIMS>> getPoints(int domainIdx)
  {
    auto* cGroup = coords_group(domainIdx);
    auto* xView = cGroup->getView("values/x");
    auto* yView = cGroup->getView("values/y");
    auto* zView = NDIMS >= 3 ? cGroup->getView("values/z") : nullptr;
    const auto ptCount = xView->getNumElements();
    assert(xView->getStride() == 1);
    assert(yView->getStride() == 1);
    assert(zView == nullptr || zView->getStride() == 1);
    double* xs = xView->getArray();
    double* ys = yView->getArray();
    double* zs = zView ? (double*)(zView->getArray()) : nullptr;

    using PointType = primal::Point<double, NDIMS>;
    axom::Array<PointType> pts;
    pts.resize(ptCount);
    for(int i = 0; i < ptCount; ++i)
    {
      pts[i][0] = xs[i];
    }
    for(int i = 0; i < ptCount; ++i)
    {
      pts[i][1] = ys[i];
    }
    if(NDIMS == 3)
    {
      for(int i = 0; i < ptCount; ++i)
      {
        pts[i][0] = zs[i];
      }
    }
    return pts;
  }

  template <typename T>
  void registerNodalScalarField(const std::string& fieldName)
  {
    for(axom::IndexType dIdx = 0; dIdx < domain_count(); ++dIdx)
    {
      auto* fld = fields_group(dIdx)->createGroup(fieldName);
      fld->createViewString("association", "vertex");
      fld->createViewString("topology", topo_group(dIdx)->getName());
      if(m_coordsAreStrided)
      {
        auto* offsets = topo_group(dIdx)->getView("elements/dims/offsets");
        auto* strides = topo_group(dIdx)->getView("elements/dims/strides");
        fld->copyView(offsets);
        fld->copyView(strides);
      }
      fld->createViewAndAllocate("values",
                                 sidre::detail::SidreTT<T>::id,
                                 numPoints(dIdx));
    }
  }

  template <typename T>
  void registerNodalVectorField(const std::string& fieldName)
  {
    const int DIM = dimension();
    for(axom::IndexType dIdx = 0; dIdx < domain_count(); ++dIdx)
    {
      const int SZ = numPoints(dIdx);

      auto* fld = fields_group(dIdx)->createGroup(fieldName);
      fld->createViewString("association", "vertex");
      fld->createViewString("topology", topo_group(dIdx)->getName());
      if(m_coordsAreStrided)
      {
        auto* offsets = topo_group(dIdx)->getView("elements/dims/offsets");
        auto* strides = topo_group(dIdx)->getView("elements/dims/strides");
        fld->copyView(offsets);
        fld->copyView(strides);
      }

      // create views into a shared buffer for the coordinates, with stride DIM
      auto* buf = domain_group(dIdx)
                    ->getDataStore()
                    ->createBuffer(sidre::detail::SidreTT<T>::id, DIM * SZ)
                    ->allocate();
      switch(DIM)
      {
      case 3:
        fld->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, DIM);
        fld->createView("values/y")->attachBuffer(buf)->apply(SZ, 1, DIM);
        fld->createView("values/z")->attachBuffer(buf)->apply(SZ, 2, DIM);
        break;
      case 2:
        fld->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, DIM);
        fld->createView("values/y")->attachBuffer(buf)->apply(SZ, 1, DIM);
        break;
      default:
        fld->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, DIM);
        break;
      }
    }
  }

  bool hasScalarField(const std::string& fieldName, int domainIdx = 0) const
  {
    auto* domain = m_group->getGroup(domainIdx);
    auto* fields = domain->getGroup("fields");
    auto has = fields->hasGroup(fieldName);
    return has;
  }

  bool hasVectorField(const std::string& fieldName, int domainIdx = 0) const
  {
    return m_group->getGroup(domainIdx)->getGroup("fields")->hasGroup(fieldName);
  }

  template <typename T>
  axom::ArrayView<T> getNodalScalarField(const std::string& fieldName,
                                         int domainIdx)
  {
    SLIC_ASSERT_MSG(
      domainIdx >= 0 && axom::IndexType(domainIdx) < domain_count(),
      axom::fmt::format("Rank {} has no domain {}, only {} domains",
                        m_rank,
                        domainIdx,
                        domain_count()));

    auto* domain = m_group->getGroup(domainIdx);
    auto* fields = domain->getGroup("fields");
    auto* field = fields->getGroup(fieldName);
    T* data =
      field ? static_cast<T*>(field->getView("values")->getVoidPtr()) : nullptr;

    return field ? axom::ArrayView<T>(data, numPoints(domainIdx))
                 : axom::ArrayView<T>();
  }

  template <typename T>
  axom::ArrayView<T> getNodalVectorField(const std::string& fieldName,
                                         int domainIdx)
  {
    SLIC_ASSERT_MSG(
      domainIdx >= 0 && axom::IndexType(domainIdx) < domain_count(),
      axom::fmt::format("Rank {} has only {} domains, no domain index {}",
                        m_rank,
                        domain_count(),
                        domainIdx));

    // Note: the implementation currently assumes that the field data is
    // interleaved, so it is safe to get a pointer to the beginning of the
    // x-coordinate's data. This will be relaxed in the future, and we will
    // need to modify this implementation accordingly.
    T* data = nullptr;
    axom::IndexType npts = 0;
    bool has =
      m_group->getGroup(domainIdx)->getGroup("fields")->hasGroup(fieldName);
    if(has)
    {
      auto xView =
        m_group->getGroup(domainIdx)->getGroup("fields")->getGroup(fieldName)->getView(
          "values/x");
      data = static_cast<T*>(xView->getVoidPtr());
      npts = xView->getNumElements();
    }
    return axom::ArrayView<T>(data, npts);
  }

  sidre::Group* getDomain(axom::IndexType domain)
  {
    return m_group->getGroup(domain);
  }
  sidre::Group* getFields(axom::IndexType domainIdx)
  {
    auto* fields = m_group->getGroup(domainIdx)->getGroup("fields");
    return fields;
  }

  /// Checks whether the blueprint is valid and prints diagnostics
  bool isValid() const
  {
    {
      conduit::Node meshNode;
      m_group->createNativeLayout(meshNode);
      conduit::Node info;
      if(!conduit::blueprint::mpi::verify("mesh", meshNode, info, MPI_COMM_WORLD))
      {
        SLIC_INFO("Invalid blueprint for particle mesh: \n" << info.to_yaml());
        slic::flushStreams();
        return false;
      }
    }
    return true;
  }

  /// Outputs the particle mesh to disk
  void saveMesh(const std::string& filename)
  {
    conduit::Node meshNode;
    m_group->createNativeLayout(meshNode);
    conduit::relay::mpi::io::blueprint::save_mesh(meshNode,
                                                  filename,
                                                  "hdf5",
                                                  MPI_COMM_WORLD);
  }

  void print_mesh_info() const
  {
    // Copy to conduit::Node.  It's output is easier to read, especially in parallel.
    conduit::Node meshNode;
    m_group->createNativeLayout(meshNode);
    meshNode.print();
  }

private:
  /// Creates blueprint stubs for this mesh
  // for the "coordset", "topologies", "fields" and "state"
  // Return the domain index created.
  axom::IndexType createBlueprintStubs()
  {
    SLIC_ASSERT(m_group != nullptr);

    auto* domainGroup = m_group->createUnnamedGroup();

    auto* coordsGroup =
      domainGroup->createGroup("coordsets")->createGroup(m_coordsetName);
    coordsGroup->createViewString("type", "explicit");
    coordsGroup->createGroup("values");

    auto* topoGroup =
      domainGroup->createGroup("topologies")->createGroup(m_topologyName);
    topoGroup->createViewString("coordset", m_coordsetName);
    topoGroup->createViewString("type", "unstructured");
    topoGroup->createViewString("elements/shape", "point");

    domainGroup->createGroup("fields");
    domainGroup->createGroup("state");

    return m_group->getNumGroups() - 1;
  }

private:
  //!@brief Whether stride/offsets are given for blueprint mesh coordinates data.
  bool m_coordsAreStrided = false;
  std::string m_topologyName;
  std::string m_coordsetName;
  /// Parent group for the entire mesh
  sidre::Group* m_group;

  int m_rank;
  int m_nranks;
  int m_dimension {-1};
};  // BlueprintParticleMesh

/**
 * Helper class to generate a mesh blueprint-conforming particle mesh for the input object.
 * The mesh is represented using a Sidre hierarchy
 */
class ObjectMeshWrapper
{
public:
  using Circle = primal::Sphere<double, 2>;

  ObjectMeshWrapper(sidre::Group* group) : m_objectMesh(group, "mesh", "coords")
  {
    SLIC_ASSERT(group != nullptr);
  }

  BlueprintParticleMesh& getParticleMesh() { return m_objectMesh; }

  /// Get a pointer to the root group for this mesh
  sidre::Group* getBlueprintGroup() const { return m_objectMesh.root_group(); }

  std::string getTopologyName() const { return m_objectMesh.getTopologyName(); }
  std::string getCoordsetName() const { return m_objectMesh.getCoordsetName(); }

  void setVerbosity(bool verbose) { m_verbose = verbose; }

  /**
   * Generates a collection of \a numPoints points along a circle.
   * Point spacing can be random (default) or uniform.
   */
  void generateCircleMesh(const Circle& circle,
                          int totalNumPoints,
                          int localDomainCount,
                          bool randomSpacing = true)
  {
    using axom::utilities::random_real;

    constexpr int DIM = 2;
    using PointType = primal::Point<double, DIM>;
    using PointArray = axom::Array<PointType>;

    int rank = m_objectMesh.getRank();
    int nranks = m_objectMesh.getNumRanks();

    // perform scan on ranks to compute totalNumPoints, thetaStart and thetaEnd
    axom::Array<int> sums(nranks, nranks);
    {
      axom::Array<int> indivDomainCounts(nranks, nranks);
      indivDomainCounts.fill(-1);
      MPI_Allgather(&localDomainCount,
                    1,
                    MPI_INT,
                    indivDomainCounts.data(),
                    1,
                    MPI_INT,
                    MPI_COMM_WORLD);

      SLIC_DEBUG_IF(m_verbose,
                    axom::fmt::format("After all gather: [{}]",
                                      axom::fmt::join(indivDomainCounts, ",")));

      sums[0] = indivDomainCounts[0];
      for(int i = 1; i < nranks; ++i)
      {
        sums[i] = sums[i - 1] + indivDomainCounts[i];
      }
      // If no rank has any domains, force last one to 1 domain.
      if(sums[nranks - 1] == 0)
      {
        sums[nranks - 1] = 1;
        if(rank == nranks - 1)
        {
          localDomainCount = 1;
        }
      }
    }

    SLIC_DEBUG_IF(
      m_verbose,
      axom::fmt::format("After scan: [{}]", axom::fmt::join(sums, ",")));

    int globalDomainCount = sums[nranks - 1];
    totalNumPoints = std::max(totalNumPoints, globalDomainCount);
    int ptsPerDomain = totalNumPoints / globalDomainCount;
    int domainsWithExtraPt = totalNumPoints % globalDomainCount;

    int myDomainBegin = rank == 0 ? 0 : sums[rank - 1];
    int myDomainEnd = sums[rank];
    assert(myDomainEnd - myDomainBegin == localDomainCount);

    double radius = circle.getRadius();
    const auto& center = circle.getCenter();
    const double avgAng = 2. * M_PI / totalNumPoints;

    for(int di = myDomainBegin; di < myDomainEnd; ++di)
    {
      int pBegin = di * ptsPerDomain + std::min(di, domainsWithExtraPt);
      int pEnd = (di + 1) * ptsPerDomain + std::min((di + 1), domainsWithExtraPt);
      int domainPointCount = pEnd - pBegin;
      PointArray pts(0, domainPointCount);

      for(int pi = pBegin; pi < pEnd; ++pi)
      {
        const double ang = randomSpacing
          ? random_real(avgAng * pBegin, avgAng * pEnd)
          : pi * avgAng;
        const double rsinT = center[1] + radius * std::sin(ang);
        const double rcosT = center[0] + radius * std::cos(ang);
        pts.push_back(PointType {rcosT, rsinT});
      }
      m_objectMesh.setPoints(di, pts);
    }

    axom::slic::flushStreams();
    SLIC_ASSERT(m_objectMesh.isValid());
  }

  /// Outputs the object mesh to disk
  void saveMesh(const std::string& filename = "object_mesh")
  {
    SLIC_INFO(
      banner(axom::fmt::format("Saving particle mesh '{}' to disk", filename)));

    m_objectMesh.saveMesh(filename);
  }

private:
  BlueprintParticleMesh m_objectMesh;
  bool m_verbose {false};
};

class QueryMeshWrapper
{
public:
  using Circle = primal::Sphere<double, 2>;

  //!@brief Construct with blueprint mesh.
  QueryMeshWrapper(sidre::Group* group, const std::string& meshFilename)
    : m_queryMesh(group)
  {
    // Test reading in multidomain mesh.
    m_queryMesh.read_blueprint_mesh(meshFilename);
    // setupParticleMesh();
  }

  BlueprintParticleMesh& getParticleMesh() { return m_queryMesh; }

  sidre::Group* getBlueprintGroup() const { return m_queryMesh.root_group(); }

  std::string getTopologyName() const { return m_queryMesh.getTopologyName(); }
  std::string getCoordsetName() const { return m_queryMesh.getCoordsetName(); }

  /// Returns an array containing the positions of the mesh vertices
  template <typename PointArray>
  PointArray getVertexPositions(int domainIdx)
  {
    // SLIC_ERROR("TODO: get a PointArray from a sidre coordset/values group.");
    sidre::Group* cvg = m_queryMesh.domain_group(domainIdx)->getGroup(
      axom::fmt::format("coordsets/{}/values", m_queryMesh.getCoordsetName()));
    int ndim = cvg->getNumViews();
    sidre::View* xv = cvg->getView("x");
    sidre::View* yv = cvg->getView("y");
    sidre::View* zv = ndim == 3 ? cvg->getView("z") : nullptr;
    axom::IndexType npts = xv->getNumElements();
    double* xp = xv->getData();
    double* yp = yv->getData();
    double* zp = zv ? (double*)(zv->getData()) : nullptr;
    double* xyzs[3] {xp, yp, zp};
    PointArray rval(npts, npts);
    for(int d = 0; d < ndim; ++d)
    {
      double* vs = xyzs[d];
      for(int i = 0; i < npts; ++i)
      {
        rval[i][d] = vs[i];
      }
    }
    typename PointArray::value_type pt;
    axom::primal::Point<double, PointArray::value_type::DIMENSION> pt1;
    return rval;
  }

  /// Saves the mesh to disk
  void saveMesh(const std::string& filename)
  {
    SLIC_INFO(
      banner(axom::fmt::format("Saving query mesh '{}' to disk", filename)));

    m_queryMesh.saveMesh(filename);
  }

  void setupParticleMesh()
  {
    {
      m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_rank");
      m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_index");
      m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_domain_index");
      m_queryMesh.registerNodalScalarField<double>("cp_distance");
      m_queryMesh.registerNodalVectorField<double>("cp_coords");
    }

    SLIC_ASSERT(m_queryMesh.isValid());
  }

  /// Prints some info about the mesh
  void print_mesh_info() { m_queryMesh.print_mesh_info(); }

  /*!
    @brief Update results from closest point search.
  */
  void update_closest_points(const conduit::Node& node)
  {
    sidre::Group* dstDomains = m_queryMesh.root_group();
    bool isMultidomain = conduit::blueprint::mesh::is_multi_domain(node);
    if(!isMultidomain)
    {
      SLIC_ASSERT(!isMultidomain ||
                  dstDomains->getNumGroups() == node.number_of_children());
    }
    const int domainCount = dstDomains->getNumGroups();
    for(int d = 0; d < domainCount; ++d)
    {
      sidre::Group& domGroup = *dstDomains->getGroup(d);
      const conduit::Node& domNode = isMultidomain ? node.child(d) : node;

      sidre::Group& dstFieldsGroup = *domGroup.getGroup("fields");
      const conduit::Node& srcFieldsNode = domNode.fetch_existing("fields");
      {
        if(!m_queryMesh.hasScalarField("cp_rank"))
        {
          m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_rank");
        }
        auto dst = dstFieldsGroup.getGroup("cp_rank");
        auto src = srcFieldsNode.fetch_existing("cp_rank");
        bool goodImport = dst->importConduitTree(src);
        SLIC_ASSERT(goodImport);
        AXOM_UNUSED_VAR(goodImport);
      }
      {
        if(!m_queryMesh.hasScalarField("cp_index"))
        {
          m_queryMesh.registerNodalScalarField<axom::IndexType>("cp_index");
        }
        auto dst = dstFieldsGroup.getGroup("cp_index");
        auto src = srcFieldsNode.fetch_existing("cp_index");
        bool goodImport = dst->importConduitTree(src);
        SLIC_ASSERT(goodImport);
        AXOM_UNUSED_VAR(goodImport);
      }
      if(srcFieldsNode.has_child("cp_domain_index"))
      {
        if(!m_queryMesh.hasScalarField("cp_domain_index"))
        {
          m_queryMesh.registerNodalScalarField<axom::IndexType>(
            "cp_domain_index");
        }
        auto src = srcFieldsNode.fetch_existing("cp_domain_index");
        auto dst = dstFieldsGroup.getGroup("cp_domain_index");
        bool goodImport = dst->importConduitTree(src);
        SLIC_ASSERT(goodImport);
        AXOM_UNUSED_VAR(goodImport);
      }
      {
        if(!m_queryMesh.hasVectorField("cp_coords"))
        {
          m_queryMesh.registerNodalVectorField<double>("cp_coords");
        }
        auto dstGroup = dstFieldsGroup.getGroup("cp_coords");
        auto srcNode = srcFieldsNode.fetch_existing("cp_coords");
        int dim = srcNode.fetch_existing("values").number_of_children();
        for(int d = 0; d < dim; ++d)
        {
          conduit::float64_array dst =
            dstGroup->getGroup("values")->getView(d)->getArray();
          const conduit::float64_array src =
            srcNode.fetch_existing("values").child(d).value();
          SLIC_ASSERT(src.number_of_elements() == dst.number_of_elements());
          int nPts = src.number_of_elements();
          for(int i = 0; i < nPts; ++i)
          {
            dst[i] = src[i];
          }
        }
      }
    }
  }

  /**
   * Check for error in the search.
   * - check that points within threshold have a closest point
   *   on the object.
   * - check that found closest-point is near its corresponding
   *   closest point on the circle (within tolerance)
   *
   * Return number of errors found on the local mesh partition.
   * Populate "error_flag" field with the number of errors, for
   * visualization.
   *
   * Randomized circle points (--random-spacing switch) can cause
   * false positives, so when it's on, distance inaccuracy is a warning
   * (not an error) for the purpose of checking.
   */
  template <int NDIMS>
  int checkClosestPoints(const Circle& circle, const Input& params)
  {
    using PointType = Circle::PointType;
    using PointArray = axom::Array<PointType>;

    m_queryMesh.registerNodalScalarField<axom::IndexType>("error_flag");

    int sumErrCount = 0;
    int sumWarningCount = 0;
    for(axom::IndexType dIdx = 0; dIdx < m_queryMesh.domain_count(); ++dIdx)
    {
      PointArray queryPts = m_queryMesh.getPoints<NDIMS>(dIdx);

      axom::ArrayView<PointType> cpCoords =
        m_queryMesh.getNodalVectorField<PointType>("cp_coords", dIdx);
      SLIC_INFO(axom::fmt::format("Closest points ({}):", cpCoords.size()));

      axom::ArrayView<axom::IndexType> cpIndices =
        m_queryMesh.getNodalScalarField<axom::IndexType>("cp_index", dIdx);

      axom::ArrayView<axom::IndexType> errorFlag =
        m_queryMesh.getNodalScalarField<axom::IndexType>("error_flag", dIdx);

      SLIC_ASSERT(queryPts.size() == cpCoords.size());
      SLIC_ASSERT(queryPts.size() == cpIndices.size());

      if(params.isVerbose())
      {
        SLIC_INFO(axom::fmt::format("Closest points ({}):", cpCoords.size()));
      }

      /*
        Allowable slack is half the arclength between 2 adjacent circle
        points.  A query point on the circle can correctly have that
        closest-distance, even though the analytical distance is zero.
        If spacing is random, distance between adjacent points is not
        predictable, leading to false positives.  We don't claim errors
        for this in when using random.
      */
      const double avgObjectRes =
        2 * M_PI * params.circleRadius / params.circlePoints;
      const double allowableSlack = avgObjectRes / 2;

      using IndexSet = slam::PositionSet<>;
      for(auto i : IndexSet(queryPts.size()))
      {
        bool errf = false;

        const auto& qPt = queryPts[i];
        const auto& cpCoord = cpCoords[i];
        double analyticalDist = std::fabs(circle.computeSignedDistance(qPt));
        const bool closestPointFound = (cpIndices[i] == -1);
        if(closestPointFound)
        {
          if(analyticalDist < params.distThreshold - allowableSlack)
          {
            errf = true;
            SLIC_INFO(
              axom::fmt::format("***Error: Query point {} ({}) is within "
                                "threshold by {} but lacks closest point.",
                                i,
                                qPt,
                                params.distThreshold - analyticalDist));
          }
        }
        else
        {
          if(analyticalDist >= params.distThreshold + allowableSlack)
          {
            errf = true;
            SLIC_INFO(
              axom::fmt::format("***Error: Query point {} ({}) is outside "
                                "threshold by {} but has closest point at {}.",
                                i,
                                qPt,
                                analyticalDist - params.distThreshold,
                                cpCoord));
          }

          if(!axom::utilities::isNearlyEqual(circle.computeSignedDistance(cpCoord),
                                             0.0))
          {
            errf = true;
            SLIC_INFO(axom::fmt::format(
              "***Error: Closest point ({}) for index {} is not on the circle.",
              cpCoords[i],
              i));
          }

          double dist = sqrt(primal::squared_distance(qPt, cpCoord));
          if(!axom::utilities::isNearlyEqual(dist, analyticalDist, allowableSlack))
          {
            if(params.randomSpacing)
            {
              ++sumWarningCount;
              SLIC_INFO(axom::fmt::format(
                "***Warning: Closest distance for index {} is {}, off by {}.",
                i,
                dist,
                dist - analyticalDist));
            }
            else
            {
              errf = true;
              SLIC_INFO(
                axom::fmt::format("***Error: Closest distance for index {} is "
                                  "{}, off by {}.",
                                  i,
                                  dist,
                                  dist - analyticalDist));
            }
          }
        }
        errorFlag[i] = errf;
        sumErrCount += errf;
      }
    }

    SLIC_INFO(axom::fmt::format(
      "Local partition has {} errors, {} warnings in closest distance results.",
      sumErrCount,
      sumWarningCount));

    return sumErrCount;
  }

private:
  BlueprintParticleMesh m_queryMesh;
};

void make_coords_contiguous(conduit::Node& coordValues)
{
  bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
  if(isInterleaved)
  {
    conduit::Node oldValues = coordValues;
    conduit::blueprint::mcarray::to_contiguous(oldValues, coordValues);
  }
}

void make_coords_interleaved(conduit::Node& coordValues)
{
  bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
  if(!isInterleaved)
  {
    conduit::Node oldValues = coordValues;
    conduit::blueprint::mcarray::to_interleaved(oldValues, coordValues);
  }
}

/// Utility function to initialize the logger
void initializeLogger()
{
  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

  slic::LogStream* logStream;

#ifdef AXOM_USE_MPI
  std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  #ifdef AXOM_USE_LUMBERJACK
  const int RLIMIT = 8;
  logStream = new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, RLIMIT, fmt);
  #else
  logStream = new slic::SynchronizedStream(&std::cout, MPI_COMM_WORLD, fmt);
  #endif
#else
  std::string fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::GenericOutputStream(&std::cout, fmt);
#endif  // AXOM_USE_MPI

  slic::addStreamToAllMsgLevels(logStream);
}

/// Utility function to finalize the logger
void finalizeLogger()
{
  if(slic::isInitialized())
  {
    slic::flushStreams();
    slic::finalize();
  }
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  initializeLogger();
  //slic::setAbortOnWarning(true);

  //---------------------------------------------------------------------------
  // Set up and parse command line arguments
  //---------------------------------------------------------------------------
  Input params;
  axom::CLI::App app {"Driver for distributed distance query"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    if(my_rank == 0)
    {
      retval = app.exit(e);
    }

    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();

    exit(retval);
  }

  // Issue warning about result-checking requiring good resolution.
  if(params.checkResults && params.randomSpacing)
  {
    SLIC_INFO(axom::fmt::format(
      "***Warning: Result-checking may yield false positive (warnings) when "
      "circle points have random spacing.  High resolution helps limit this."
      "We recommend at least 500 points for each radius length unit."));
  }

  constexpr int DIM = 2;

  using PointType = primal::Point<double, DIM>;
  using PointArray = axom::Array<PointType>;
  using IndexSet = slam::PositionSet<>;
  using Circle = primal::Sphere<double, DIM>;

#if defined(AXOM_USE_UMPIRE)
  //---------------------------------------------------------------------------
  // Memory resource.  For testing, choose device memory if appropriate.
  //---------------------------------------------------------------------------
  const std::string umpireResourceName = params.policy == RuntimePolicy::seq
    ? "HOST"
    :
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
    params.policy == RuntimePolicy::omp ? "HOST" :
  #endif
  #if defined(UMPIRE_ENABLE_DEVICE)
                                        "DEVICE"
  #elif defined(UMPIRE_ENABLE_UM)
    "UM"
  #elif defined(UMPIRE_ENABLE_PINNED)
    "PINNED"
  #else
    "HOST"
  #endif
    ;
  auto& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator umpireAllocator = rm.getAllocator(umpireResourceName);
#endif

  //---------------------------------------------------------------------------
  // Load/generate object mesh
  //---------------------------------------------------------------------------
  const Circle circle(
    PointType(params.circleCenter.data(), params.circleCenter.size()),
    params.circleRadius);

  sidre::DataStore dataStore;

  ObjectMeshWrapper objectMeshWrapper(
    dataStore.getRoot()->createGroup("object_mesh", true));
  objectMeshWrapper.setVerbosity(params.isVerbose());

  {
    SLIC_ASSERT(params.objDomainCountRange[1] >= params.objDomainCountRange[0]);
    const unsigned int omin = params.objDomainCountRange[0];
    const unsigned int omax = params.objDomainCountRange[1];
    const double prob = axom::utilities::random_real(0., 1.);
    int localDomainCount = omin + int(0.5 + prob * (omax - omin));
    objectMeshWrapper.generateCircleMesh(circle,
                                         params.circlePoints,
                                         localDomainCount,
                                         params.randomSpacing);
  }

  SLIC_INFO_IF(
    params.isVerbose(),
    axom::fmt::format("Object mesh has {} points in {} domains locally",
                      objectMeshWrapper.getParticleMesh().numPoints(),
                      objectMeshWrapper.getParticleMesh().domain_count()));

  objectMeshWrapper.saveMesh(params.objectFile);
  slic::flushStreams();

  //---------------------------------------------------------------------------
  // Load computational mesh and generate a particle mesh over its nodes
  // These will be used to query the closest points on the object mesh(es)
  //---------------------------------------------------------------------------
  QueryMeshWrapper queryMeshWrapper(
    dataStore.getRoot()->createGroup("queryMesh", true),
    params.meshFile);
  // queryMeshWrapper.print_mesh_info();

  SLIC_INFO_IF(
    params.isVerbose(),
    axom::fmt::format("Query mesh has {} points in {} domains locally",
                      queryMeshWrapper.getParticleMesh().numPoints(),
                      queryMeshWrapper.getParticleMesh().domain_count()));
  slic::flushStreams();

  auto getIntMinMax = [](int inVal, int& minVal, int& maxVal, int& sumVal) {
    MPI_Allreduce(&inVal, &minVal, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&inVal, &maxVal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&inVal, &sumVal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  };

  // Output some global mesh size stats
  {
    int mn, mx, sum;
    getIntMinMax(objectMeshWrapper.getParticleMesh().numPoints(), mn, mx, sum);
    SLIC_INFO(axom::fmt::format(
      "Object mesh has {{min:{}, max:{}, sum:{}, avg:{}}} points",
      mn,
      mx,
      sum,
      (double)sum / num_ranks));
  }
  {
    int mn, mx, sum;
    getIntMinMax(objectMeshWrapper.getParticleMesh().domain_count(), mn, mx, sum);
    SLIC_INFO(axom::fmt::format(
      "Object mesh has {{min:{}, max:{}, sum:{}, avg:{}}} domains",
      mn,
      mx,
      sum,
      (double)sum / num_ranks));
  }
  {
    int mn, mx, sum;
    getIntMinMax(queryMeshWrapper.getParticleMesh().numPoints(), mn, mx, sum);
    SLIC_INFO(axom::fmt::format(
      "Query mesh has {{min:{}, max:{}, sum:{}, avg:{}}} points",
      mn,
      mx,
      sum,
      (double)sum / num_ranks));
  }
  {
    int mn, mx, sum;
    getIntMinMax(queryMeshWrapper.getParticleMesh().domain_count(), mn, mx, sum);
    SLIC_INFO(axom::fmt::format(
      "Query mesh has {{min:{}, max:{}, sum:{}, avg:{}}} domains",
      mn,
      mx,
      sum,
      (double)sum / num_ranks));
  }

  slic::flushStreams();

  //---------------------------------------------------------------------------
  // Initialize spatial index for querying points, and run query
  //---------------------------------------------------------------------------

  auto init_str =
    banner(axom::fmt::format("Initializing BVH tree over {} points",
                             params.circlePoints));

  axom::utilities::Timer initTimer(false);
  axom::utilities::Timer queryTimer(false);

  // Convert blueprint representation from sidre to conduit
  conduit::Node objectMeshNode;
  if(objectMeshWrapper.getParticleMesh().numPoints() > 0)
  {
    objectMeshWrapper.getBlueprintGroup()->createNativeLayout(objectMeshNode);
  }

  // Put sidre data into Conduit Node.
  conduit::Node queryMeshNode;
  queryMeshWrapper.getBlueprintGroup()->createNativeLayout(queryMeshNode);

  // To test with contiguous and interleaved coordinate storage,
  // make half them contiguous.
  for(int di = 0; di < objectMeshNode.number_of_children(); ++di)
  {
    auto& dom = objectMeshNode.child(di);
    if((my_rank + di) % 2 == 1)
    {
      make_coords_contiguous(dom.fetch_existing("coordsets/coords/values"));
    }
  }
  for(int di = 0; di < queryMeshNode.number_of_children(); ++di)
  {
    auto& dom = queryMeshNode.child(di);
    if((my_rank + di) % 2 == 1)
    {
      make_coords_contiguous(dom.fetch_existing("coordsets/coords/values"));
    }
  }

  // Create distributed closest point query object and set some parameters
  quest::DistributedClosestPoint query;
  query.setRuntimePolicy(params.policy);
#if defined(AXOM_USE_UMPIRE)
  query.setAllocatorID(umpireAllocator.getId());
#endif
  query.setMpiCommunicator(MPI_COMM_WORLD, true);
  query.setDimension(DIM);
  query.setVerbosity(params.isVerbose());
  query.setDistanceThreshold(params.distThreshold);
  // To test support for single-domain format, use single-domain when possible.
  query.setObjectMesh(
    objectMeshNode.number_of_children() == 1 ? objectMeshNode[0] : objectMeshNode,
    objectMeshWrapper.getTopologyName());

  // Build the spatial index over the object on each rank
  SLIC_INFO(init_str);
  slic::flushStreams();
  initTimer.start();
  query.generateBVHTree();
  initTimer.stop();

  // Run the distributed closest point query over the nodes of the computational mesh
  // To test support for single-domain format, use single-domain when possible.
  slic::flushStreams();
  queryTimer.start();
  query.computeClosestPoints(
    queryMeshNode.number_of_children() == 1 ? queryMeshNode[0] : queryMeshNode,
    queryMeshWrapper.getTopologyName());
  queryTimer.stop();

  auto getDoubleMinMax =
    [](double inVal, double& minVal, double& maxVal, double& sumVal) {
      MPI_Allreduce(&inVal, &minVal, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&inVal, &maxVal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&inVal, &sumVal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    };

  // Output some timing stats
  {
    double minInit, maxInit, sumInit;
    getDoubleMinMax(initTimer.elapsedTimeInSec(), minInit, maxInit, sumInit);

    double minQuery, maxQuery, sumQuery;
    getDoubleMinMax(queryTimer.elapsedTimeInSec(), minQuery, maxQuery, sumQuery);

    SLIC_INFO(axom::fmt::format(
      "Initialization with policy {} took {{avg:{}, min:{}, max:{}}} seconds",
      params.policy,
      sumInit / num_ranks,
      minInit,
      maxInit));
    SLIC_INFO(axom::fmt::format(
      "Query with policy {} took {{avg:{}, min:{}, max:{}}} seconds",
      params.policy,
      sumQuery / num_ranks,
      minQuery,
      maxQuery));
  }
  slic::flushStreams();
  queryMeshWrapper.update_closest_points(queryMeshNode);

  int errCount = 0;
  int localErrCount = 0;
  if(params.checkResults)
  {
    localErrCount = queryMeshWrapper.checkClosestPoints<DIM>(circle, params);
  }
  MPI_Allreduce(&localErrCount, &errCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  //---------------------------------------------------------------------------
  // Transform closest points to distances and directions
  //---------------------------------------------------------------------------
  using primal::squared_distance;

  auto& queryMesh = queryMeshWrapper.getParticleMesh();
  PointType nowhere(std::numeric_limits<double>::signaling_NaN());
  const double nodist = std::numeric_limits<double>::signaling_NaN();
  queryMesh.registerNodalScalarField<double>("distance");
  queryMesh.registerNodalVectorField<double>("direction");
  for(axom::IndexType di = 0; di < queryMesh.domain_count(); ++di)
  {
    auto cpCoords = queryMesh.getNodalVectorField<PointType>("cp_coords", di);

    auto cpIndices =
      queryMesh.getNodalScalarField<axom::IndexType>("cp_index", di);

    PointArray qPts = queryMeshWrapper.getVertexPositions<PointArray>(di);
    axom::ArrayView<double> distances =
      queryMesh.getNodalScalarField<double>("distance", di);
    axom::ArrayView<PointType> directions =
      queryMesh.getNodalVectorField<PointType>("direction", di);
    axom::IndexType ptCount = queryMeshWrapper.getParticleMesh().numPoints(di);
    for(auto ptIdx : IndexSet(ptCount))
    {
      const bool has_cp = cpIndices[ptIdx] >= 0;
      const PointType& cp = has_cp ? cpCoords[ptIdx] : nowhere;
      distances[ptIdx] =
        has_cp ? sqrt(squared_distance(qPts[ptIdx], cp)) : nodist;
      directions[ptIdx] =
        PointType(has_cp ? (cp - qPts[ptIdx]).array() : nowhere.array());
    }
  }

  // queryMeshNode.print();
  queryMeshNode.reset();

  queryMeshWrapper.saveMesh(params.distanceFile);

  if(errCount)
  {
    SLIC_INFO(axom::fmt::format(" Error exit: {} errors found.", errCount));
  }
  else
  {
    SLIC_INFO("Normal exit.");
  }

  finalizeLogger();
  MPI_Finalize();

  return errCount != 0;
}

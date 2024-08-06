// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/DiscreteShape.hpp"
#include "axom/quest/Discretize.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/klee/GeometryOperators.hpp"
#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/quest/interface/internal/QuestHelpers.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace axom
{
namespace quest
{
namespace internal
{
/*!
 * \brief Implementation of a GeometryOperatorVisitor for processing klee shape operators
 *
 * This class extracts the matrix form of supported operators and marks the operator as unvalid otherwise
 * To use, check the \a isValid() function after visiting and then call the \a getMatrix() function.
 */
class AffineMatrixVisitor : public klee::GeometryOperatorVisitor
{
public:
  AffineMatrixVisitor() : m_matrix(4, 4) { }

  void visit(const klee::Translation& translation) override
  {
    m_matrix = translation.toMatrix();
    m_isValid = true;
  }
  void visit(const klee::Rotation& rotation) override
  {
    m_matrix = rotation.toMatrix();
    m_isValid = true;
  }
  void visit(const klee::Scale& scale) override
  {
    m_matrix = scale.toMatrix();
    m_isValid = true;
  }
  void visit(const klee::UnitConverter& converter) override
  {
    m_matrix = converter.toMatrix();
    m_isValid = true;
  }

  void visit(const klee::CompositeOperator&) override
  {
    SLIC_WARNING("CompositeOperator not supported for Shaper query");
    m_isValid = false;
  }
  void visit(const klee::SliceOperator&) override
  {
    SLIC_WARNING("SliceOperator not yet supported for Shaper query");
    m_isValid = false;
  }

  const numerics::Matrix<double>& getMatrix() const { return m_matrix; }
  bool isValid() const { return m_isValid; }

private:
  bool m_isValid {false};
  numerics::Matrix<double> m_matrix;
};

}  // end namespace internal

// These were needed for linking - but why? They are constexpr.
constexpr int DiscreteShape::DEFAULT_SAMPLES_PER_KNOT_SPAN;
constexpr double DiscreteShape::MINIMUM_PERCENT_ERROR;
constexpr double DiscreteShape::MAXIMUM_PERCENT_ERROR;
constexpr double DiscreteShape::DEFAULT_VERTEX_WELD_THRESHOLD;

DiscreteShape::DiscreteShape(const axom::klee::Shape &shape,
                             const std::string& prefixPath,
                             axom::sidre::Group* parentGroup)
  : m_shape(shape)
  , m_sidreGroup(nullptr)
{
  setPrefixPath(prefixPath);
  setParentGroup(parentGroup);
}

void DiscreteShape::setParentGroup(axom::sidre::Group* parentGroup)
{
  if (parentGroup)
  {
    std::ostringstream os;
    os << this;
    std::string myGroupName = os.str();
    m_sidreGroup = parentGroup->createGroup(myGroupName);
  }
}

void DiscreteShape::clearInternalData()
{
  m_meshRep.reset();
  if (m_sidreGroup)
  {
    m_sidreGroup->getParent()->destroyGroupAndData(m_sidreGroup->getIndex());
    m_sidreGroup = nullptr;
  }
}

std::shared_ptr<mint::Mesh> DiscreteShape::createMeshRepresentation()
{
  if (m_meshRep) { return m_meshRep; }

  const axom::klee::Geometry& geometry = m_shape.getGeometry();
  const auto& geometryFormat = geometry.getFormat();

  if (geometryFormat == "memory-blueprint")
  {
    // TODO: For readability, move most of this into a function.
    // Put the in-memory geometry in m_meshRep.
    axom::sidre::Group* inputGroup = geometry.getBlueprintMesh();
    axom::sidre::Group* rootGroup = inputGroup->getDataStore()->getRoot();

    std::string modName = inputGroup->getName() + "_modified";
    while (rootGroup->hasGroup(modName)) { modName = modName + "-"; }

    axom::sidre::Group* modGroup = rootGroup->createGroup(modName);
    int allocID = inputGroup->getDefaultAllocatorID();
    modGroup->deepCopyGroup(inputGroup, allocID);

    m_meshRep.reset(axom::mint::getMesh(modGroup->getGroup(inputGroup->getName()),
                                        m_shape.getGeometry().getBlueprintTopology()));

    // Transform the coordinates of the linearized mesh.
    applyTransforms();
  }
  else if ( geometryFormat == "sphere3D" )
  {
    /*
      Discretize the sphere into m_meshRep and apply transforms:
      1. Discretize the sphere into a set of Octahedra.
      2. Split each Octahedron into 8 tets.
      3. Insert tets into tet mesh.
      4. Transform the tet mesh.
      We can reduce memory use by eliminating repeated points if it
      is a problem.
    */
    const auto& sphere = geometry.getSphere();
    using TetType = axom::primal::Tetrahedron<double, 3>;
    axom::Array<OctType> octs;
    int octCount = 0;
    axom::quest::discretize(sphere, geometry.getGenerationCount(), octs, octCount);

    constexpr int TETS_PER_OCT = 8;
    constexpr int NODES_PER_TET = 4;
    const axom::IndexType tetCount = octCount * TETS_PER_OCT;
    const axom::IndexType nodeCount = tetCount * NODES_PER_TET;

    axom::Array<Point3D> nodeCoords(nodeCount, nodeCount);
    axom::Array<axom::IndexType, 2> connectivity(
      axom::StackArray<axom::IndexType, 2>{tetCount, NODES_PER_TET});

    auto nodeCoordsView = nodeCoords.view();
    auto connectivityView = connectivity.view();
    axom::for_all<axom::SEQ_EXEC>(
      octCount,
      AXOM_LAMBDA(axom::IndexType octIdx) {
        TetType tetsInOct[TETS_PER_OCT];
        axom::primal::split(octs[octIdx], tetsInOct);
        for ( int iTet = 0; iTet < TETS_PER_OCT; ++iTet )
        {
          axom::IndexType tetIdx = octIdx*TETS_PER_OCT + iTet;
          for ( int iNode = 0; iNode < NODES_PER_TET; ++iNode )
          {
            axom::IndexType nodeIdx = tetIdx*NODES_PER_TET + iNode;
            nodeCoordsView[nodeIdx] = tetsInOct[iTet][iNode];
            connectivityView[tetIdx][iNode] = nodeIdx;
          }
        }
      });

    // TODO: Set this up with sidre::Group to have data on device.
    axom::mint::UnstructuredMesh<axom::mint::Topology::SINGLE_SHAPE>* tetMesh = nullptr;
    if (m_sidreGroup != nullptr)
    {
      tetMesh = new axom::mint::UnstructuredMesh<axom::mint::Topology::SINGLE_SHAPE>(
        3,
        axom::mint::CellType::TET,
        m_sidreGroup,
        nodeCount,
        tetCount);
    }
    else
    {
      tetMesh = new axom::mint::UnstructuredMesh<axom::mint::Topology::SINGLE_SHAPE>(
        3,
        axom::mint::CellType::TET,
        nodeCount,
        tetCount);
    }
    tetMesh->appendNodes((double*)nodeCoords.data(), nodeCount);
    tetMesh->appendCells(connectivity.data(), tetCount);
    m_meshRep.reset(tetMesh);

    applyTransforms();
  }

  if (m_meshRep)
  {
    return m_meshRep;
  }

  if(!m_shape.getGeometry().hasGeometry())
  {
    // If shape has no geometry, there's nothing to discretize.
    SLIC_WARNING(
      axom::fmt::format("Current shape '{}' of material '{}' has no geometry",
                        m_shape.getName(),
                        m_shape.getMaterial()));
    return m_meshRep;
  }

  // We handled all the non-file formats.  The rest are file formats.
  const std::string& file_format = geometryFormat;

  std::string shapePath = resolvePath();
  SLIC_INFO("Reading file: " << shapePath << "...");

  // Initialize revolved volume.
  m_revolvedVolume = 0.;

  if(utilities::string::endsWith(shapePath, ".stl"))
  {
    SLIC_ASSERT_MSG(
      file_format == "stl",
      axom::fmt::format(" '{}' format requires .stl file type", file_format));

    axom::mint::Mesh* meshRep = nullptr;
    quest::internal::read_stl_mesh(shapePath, meshRep, m_comm);
    m_meshRep.reset(meshRep);
    // Transform the coordinates of the linearized mesh.
    applyTransforms();
  }
  else if(utilities::string::endsWith(shapePath, ".proe"))
  {
    SLIC_ASSERT_MSG(
      file_format == "proe",
      axom::fmt::format(" '{}' format requires .proe file type", file_format));

    axom::mint::Mesh* meshRep = nullptr;
    quest::internal::read_pro_e_mesh(shapePath, meshRep, m_comm);
    m_meshRep.reset(meshRep);
  }
#ifdef AXOM_USE_C2C
  else if(utilities::string::endsWith(shapePath, ".contour"))
  {
    SLIC_ASSERT_MSG(
      file_format == "c2c",
      axom::fmt::format(" '{}' format requires .contour file type", file_format));

    // Get the transforms that are being applied to the mesh. Get them
    // as a single concatenated matrix.
    auto transform = getTransforms();

    // Pass in the transform so any transformations can figure into
    // computing the revolved volume.
    axom::mint::Mesh* meshRep = nullptr;
    if(m_refinementType == DiscreteShape::RefinementDynamic &&
       m_percentError > MINIMUM_PERCENT_ERROR)
    {
      quest::internal::read_c2c_mesh_non_uniform(shapePath,
                                                 transform,
                                                 m_percentError,
                                                 m_vertexWeldThreshold,
                                                 meshRep,
                                                 m_revolvedVolume,  // output arg
                                                 m_comm);
    }
    else
    {
      quest::internal::read_c2c_mesh_uniform(shapePath,
                                             transform,
                                             m_samplesPerKnotSpan,
                                             m_vertexWeldThreshold,
                                             meshRep,
                                             m_revolvedVolume,  // output arg
                                             m_comm);
    }
    m_meshRep.reset(meshRep);

    // Transform the coordinates of the linearized mesh.
    applyTransforms();
  }
#endif
  else
  {
    SLIC_ERROR(
      axom::fmt::format("Unsupported filetype for this Axom configuration. "
                        "Provided file was '{}', with format '{}'",
                        shapePath,
                        file_format));
  }

  return m_meshRep;
}

void DiscreteShape::applyTransforms()
{
  numerics::Matrix<double> transformation = getTransforms();

  // Apply transformation to coordinates of each vertex in mesh
  if(!transformation.isIdentity())
  {
    const int spaceDim = m_meshRep->getDimension();
    const int numSurfaceVertices = m_meshRep->getNumberOfNodes();
    double* x = m_meshRep->getCoordinateArray(mint::X_COORDINATE);
    double* y = m_meshRep->getCoordinateArray(mint::Y_COORDINATE);
    double* z = spaceDim > 2
      ? m_meshRep->getCoordinateArray(mint::Z_COORDINATE)
      : nullptr;

    double xformed[4];
    for(int i = 0; i < numSurfaceVertices; ++i)
    {
      double coords[4] = {x[i], y[i], (z == nullptr ? 0. : z[i]), 1.};
      numerics::matrix_vector_multiply(transformation, coords, xformed);
      x[i] = xformed[0];
      y[i] = xformed[1];
      if(z != nullptr)
      {
        z[i] = xformed[2];
      }
    }
  }
}

numerics::Matrix<double> DiscreteShape::getTransforms() const
{
  const auto identity4x4 = numerics::Matrix<double>::identity(4);
  numerics::Matrix<double> transformation(identity4x4);
  auto& geometryOperator = m_shape.getGeometry().getGeometryOperator();
  if (geometryOperator) {
    auto composite =
      std::dynamic_pointer_cast<const klee::CompositeOperator>(geometryOperator);
    if(composite)
    {
      // Concatenate the transformations

      // Why don't we multiply the matrices in CompositeOperator::addOperator()?
      // Why keep the matrices factored and multiply them here repeatedly?
      // Combining them would also avoid this if-else logic.  BTNG
      for(auto op : composite->getOperators())
      {
        // Use visitor pattern to extract the affine matrix from supported operators
        internal::AffineMatrixVisitor visitor;
        op->accept(visitor);
        if(!visitor.isValid())
        {
          continue;
        }
        const auto& matrix = visitor.getMatrix();
        numerics::Matrix<double> res(identity4x4);
        numerics::matrix_multiply(matrix, transformation, res);
        transformation = res;
      }
    }
    else
    {
      internal::AffineMatrixVisitor visitor;
      geometryOperator->accept(visitor);
      if (visitor.isValid())
      {
        transformation = visitor.getMatrix();
      }
    }
  }
  return transformation;
}

void DiscreteShape::setSamplesPerKnotSpan(int nSamples)
{
  using axom::utilities::clampLower;
  SLIC_WARNING_IF(
    nSamples < 1,
    axom::fmt::format(
      "Samples per knot span must be at least 1. Provided value was {}",
      nSamples));

  m_samplesPerKnotSpan = clampLower(nSamples, 1);
}

void DiscreteShape::setVertexWeldThreshold(double threshold)
{
  SLIC_WARNING_IF(
    threshold <= 0.,
    axom::fmt::format(
      "Vertex weld threshold should be positive Provided value was {}",
      threshold));

  m_vertexWeldThreshold = threshold;
}

void DiscreteShape::setPercentError(double percent)
{
  using axom::utilities::clampVal;
  SLIC_WARNING_IF(
    percent <= MINIMUM_PERCENT_ERROR,
    axom::fmt::format("Percent error must be greater than {}. Provided value "
                      "was {}. Dynamic refinement will not be used.",
                      MINIMUM_PERCENT_ERROR,
                      percent));
  SLIC_WARNING_IF(percent > MAXIMUM_PERCENT_ERROR,
                  axom::fmt::format(
                    "Percent error must be less than {}. Provided value was {}",
                    MAXIMUM_PERCENT_ERROR,
                    percent));
  if(percent <= MINIMUM_PERCENT_ERROR)
  {
    m_refinementType = DiscreteShape::RefinementUniformSegments;
  }
  m_percentError =
    clampVal(percent, MINIMUM_PERCENT_ERROR, MAXIMUM_PERCENT_ERROR);
}

std::string DiscreteShape::resolvePath() const
{
  const std::string& geomPath = m_shape.getGeometry().getPath();
  if(geomPath[0] == '/')
  {
    return geomPath;
  }
  if(m_prefixPath.empty())
  {
    throw std::logic_error("Relative geometry path requires a parent path.");
  }
  std::string dir;
  utilities::filesystem::getDirName(dir, m_prefixPath);
  return utilities::filesystem::joinPath(dir, geomPath);
}

}  // namespace klee
}  // namespace axom

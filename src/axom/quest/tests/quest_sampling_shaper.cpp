// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \brief Unit tests for quest's SamplingShaper class and associated replacement rules.
 */

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/klee.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"
#include "axom/sidre.hpp"
#include "axom/slic.hpp"
#include "axom/quest/SamplingShaper.hpp"

#ifndef AXOM_USE_MFEM
  #error "Quest's SamplingShaper tests on mfem meshes require mfem library."
#else
  #include "mfem.hpp"
#endif

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

#include <cmath>
#include <string>
#include <iostream>

namespace klee = axom::klee;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace sidre = axom::sidre;
namespace slic = axom::slic;

namespace
{
const std::string unit_circle_contour =
  "piece = circle(origin=(0cm, 0cm), radius=1cm, start=0deg, end=360deg)";

constexpr bool very_verbose_output = false;

/// RAII utility class to write a file at construction time and remove it
/// once the instance is out of scope
class ScopedTemporaryFile
{
public:
  ScopedTemporaryFile(const std::string& filename, const std::string& contents)
    : m_filename(filename)
  {
    std::ofstream ofs(m_filename.c_str(), std::ios::out);
    ofs << contents;
  }

  ~ScopedTemporaryFile()
  {
    axom::utilities::filesystem::removeFile(m_filename);
  }

  const std::string& getFileName() const { return m_filename; }

  std::string getFileContents() const
  {
    std::ifstream ifs(m_filename.c_str(), std::ios::in);
    std::stringstream buffer;
    buffer << ifs.rdbuf();

    return buffer.str();
  }

private:
  const std::string m_filename;
};
}  // namespace

/// Test fixture for SamplingShaper tests on MFEM meshes
class SamplingShaperTest : public ::testing::Test
{
public:
  using Point2D = primal::Point<double, 2>;
  using BBox2D = primal::BoundingBox<double, 2>;

public:
  SamplingShaperTest() : m_dc("test", nullptr, true) { }

  virtual ~SamplingShaperTest() { }

  virtual void SetUp()
  {
    const int dim = 2;
    const int polynomialOrder = 2;
    const double lo[] = {-2., -2.};
    const double hi[] = {2., 2.};
    const int celldims[] = {64, 64};

    // memory for mesh will be managed by data collection
    auto mesh =
      new mfem::Mesh(mfem::Mesh::MakeCartesian2D(celldims[0],
                                                 celldims[1],
                                                 mfem::Element::QUADRILATERAL,
                                                 true,
                                                 hi[0] - lo[0],
                                                 hi[1] - lo[1]));

    const int NE = mesh->GetNE();
    const int NV = mesh->GetNV();

    // Offset the mesh to lie w/in the bounding box
    for(int i = 0; i < NV; ++i)
    {
      double* v = mesh->GetVertex(i);
      for(int d = 0; d < dim; ++d)
      {
        v[d] += lo[d];
      }
    }

    // Set element attributes based on quadrant where centroid is located
    // These will be used later in some cases when setting volume fractions
    mfem::Array<int> v;
    for(int i = 0; i < NE; ++i)
    {
      mesh->GetElementVertices(i, v);
      BBox2D bbox;
      for(int j = 0; j < v.Size(); ++j)
      {
        bbox.addPoint(Point2D(mesh->GetVertex(v[j]), 2));
      }

      auto centroid = bbox.getCentroid();
      if(centroid[0] >= 0 && centroid[1] >= 0)
      {
        mesh->SetAttribute(i, 1);
      }
      else if(centroid[0] >= 0 && centroid[1] < 0)
      {
        mesh->SetAttribute(i, 2);
      }
      else if(centroid[0] < 0 && centroid[1] >= 0)
      {
        mesh->SetAttribute(i, 3);
      }
      else
      {
        mesh->SetAttribute(i, 4);
      }
    }

    mesh->SetCurvature(polynomialOrder);

    m_dc.SetOwnData(true);
    m_dc.SetMeshNodesName("positions");
    m_dc.SetMesh(mesh);

#ifdef AXOM_USE_MPI
    m_dc.SetComm(MPI_COMM_WORLD);
#endif
  }

  sidre::MFEMSidreDataCollection& getDC() { return m_dc; }
  mfem::Mesh& getMesh() { return *m_dc.GetMesh(); }

  /// parse and validate the Klee shapefile; fail the test if invalid
  void validateShapeFile(const std::string& shapefile)
  {
    axom::klee::ShapeSet shapeSet;

    try
    {
      shapeSet = axom::klee::readShapeSet(shapefile);
    }
    catch(axom::klee::KleeError& error)
    {
      std::vector<std::string> errs;
      for(auto verificationError : error.getErrors())
      {
        errs.push_back(
          axom::fmt::format(" - '{}': {}",
                            static_cast<std::string>(verificationError.path),
                            verificationError.message));
      }

      if(!errs.empty())
      {
        SLIC_WARNING(
          axom::fmt::format("Error during parsing klee input file '{}'. "
                            "Found the following errors:\n{}",
                            shapefile,
                            axom::fmt::join(errs, "\n")));
        FAIL();
      }
    }
  }

  /// Runs the shaping query over a shapefile using the supplie initial volume fractions
  void runShaping(const std::string& shapefile,
                  const std::map<std::string, mfem::GridFunction*>& init_vf_map = {})
  {
    SLIC_INFO_IF(very_verbose_output,
                 axom::fmt::format("Reading shape set from {}", shapefile));
    klee::ShapeSet shapeSet(klee::readShapeSet(shapefile));

    SLIC_INFO_IF(very_verbose_output, axom::fmt::format("Shaping materials..."));
    quest::SamplingShaper shaper(shapeSet, &m_dc);
    shaper.setVerbosity(very_verbose_output);

    if(!init_vf_map.empty())
    {
      shaper.importInitialVolumeFractions(init_vf_map);
    }

    if(very_verbose_output)
    {
      shaper.printRegisteredFieldNames("*** After importing volume fractions");
    }

    const auto shapeDim = shapeSet.getDimensions();
    for(const auto& shape : shapeSet.getShapes())
    {
      SLIC_INFO_IF(very_verbose_output,
                   axom::fmt::format("\tshape {} -> material {}",
                                     shape.getName(),
                                     shape.getMaterial()));

      shaper.loadShape(shape);
      shaper.prepareShapeQuery(shapeDim, shape);
      shaper.runShapeQuery(shape);
      shaper.applyReplacementRules(shape);
      shaper.finalizeShapeQuery();
    }

    shaper.adjustVolumeFractions();

    if(very_verbose_output)
    {
      shaper.printRegisteredFieldNames("*** After shaping volume fractions");
    }
  }

  BBox2D meshBoundingBox()
  {
    mfem::Vector bbmin, bbmax;
    getMesh().GetBoundingBox(bbmin, bbmax);

    return BBox2D(Point2D(bbmin.GetData()), Point2D(bbmax.GetData()));
  }

  // Computes the total volume of the associated volume fraction grid function
  double gridFunctionVolume(const std::string& name)
  {
    mfem::GridFunction* gf = m_dc.GetField(name);

    mfem::LinearForm vol_form(gf->FESpace());
    mfem::ConstantCoefficient one(1.0);
    vol_form.AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    vol_form.Assemble();

    return *gf * vol_form;
  }

  /// Registers and allocates a volume fraction grid function within the datastore
  mfem::GridFunction* registerVolFracGridFunction(const std::string& name,
                                                  int vfOrder = 2)
  {
    SLIC_ASSERT(!m_dc.HasField(name));

    auto& mesh = getMesh();
    const int dim = mesh.Dimension();

    // create grid function
    auto* coll =
      new mfem::L2_FECollection(vfOrder, dim, mfem::BasisType::Positive);
    auto* fes = new mfem::FiniteElementSpace(&mesh, coll);
    auto* vf = new mfem::GridFunction(fes);
    vf->MakeOwner(coll);

    // allocate grid function via sidre
    const int sz = fes->GetVSize();
    mfem::Vector v(m_dc.AllocNamedBuffer(name, sz)->getData(), sz);
    vf->MakeRef(fes, v, 0);

    // register grid function w/ data collection
    m_dc.RegisterField(name, vf);

    return vf;
  }

  /** 
   * \brief Initializes the values of the DOFs of a volume fraction grid function
   * using a provided lambda w/ parameters for the cell index, DOF position and cell attribute
   * The signature of DOFInitializer is [](int idx, Point2D& pt, int attribute) -> double
   */
  template <typename DOFInitializer>
  void initializeVolFracGridFunction(mfem::GridFunction* vf,
                                     DOFInitializer&& dof_initializer)
  {
    auto& mesh = this->getMesh();
    const int dim = mesh.Dimension();
    const int NE = mesh.GetNE();

    // Assume all elements have the same integration rule
    const auto* fes = vf->FESpace();
    auto* fe = fes->GetFE(0);
    auto& ir = fe->GetNodes();
    const int nq = ir.GetNPoints();

    // Get positions of DOFs
    mfem::DenseTensor pos_coef(dim, nq, NE);
    {
      const auto* geomFactors =
        mesh.GetGeometricFactors(ir, mfem::GeometricFactors::COORDINATES);

      // Rearrange positions
      for(int i = 0; i < NE; ++i)
      {
        for(int j = 0; j < dim; ++j)
        {
          for(int k = 0; k < nq; ++k)
          {
            pos_coef(j, k, i) = geomFactors->X((i * nq * dim) + (j * nq) + k);
          }
        }
      }
    }

    // Initialize volume fraction DOFs using passed in lambda based on cell index, DOF position and attribute
    mfem::Vector res(nq);
    mfem::Array<int> dofs;
    for(int idx = 0; idx < NE; ++idx)
    {
      const int attr = mesh.GetAttribute(idx);

      mfem::DenseMatrix& m = pos_coef(idx);
      for(int p = 0; p < nq; ++p)
      {
        const Point2D pt(m.GetColumn(p), dim);
        res(p) = dof_initializer(idx, pt, attr);
      }

      fes->GetElementDofs(idx, dofs);
      vf->SetSubVector(dofs, res);
    }
  }

  /// Helper to check integrated volume of a volume fraction grid fucntion
  void checkExpectedVolumeFractions(const std::string& material_name,
                                    double expected_volume,
                                    double EPS = 1e-2)
  {
    auto vf_name = axom::fmt::format("vol_frac_{}", material_name);

    EXPECT_TRUE(m_dc.HasField(vf_name)) << axom::fmt::format(
      "Did not have expected volume fraction '{:.4}' for material '{}'",
      material_name,
      vf_name);

    const double actual_volume = this->gridFunctionVolume(vf_name);
    SLIC_INFO(axom::fmt::format(
      "Shaped volume fraction of '{}' is {:.4}  (expected: {:.4})",
      material_name,
      actual_volume,
      expected_volume));

    EXPECT_NEAR(expected_volume, actual_volume, EPS);
  }

protected:
  sidre::MFEMSidreDataCollection m_dc;
};

//-----------------------------------------------------------------------------

TEST(ScopedTemporaryFile, basic)
{
  using axom::utilities::filesystem::pathExists;

  const std::string filename = "previously_nonexistent.file";
  const std::string contents = "file contents!";

  // File does not exist before entering scope
  EXPECT_FALSE(pathExists(filename));

  // File is created and exists within the scope
  {
    ScopedTemporaryFile test_file(filename, contents);

    EXPECT_EQ(test_file.getFileName(), filename);
    ASSERT_TRUE(pathExists(test_file.getFileName()));

    // check contents
    EXPECT_EQ(contents, test_file.getFileContents());
  }

  // File no longer exists outside the scope
  EXPECT_FALSE(pathExists(filename));
}

//-----------------------------------------------------------------------------

TEST_F(SamplingShaperTest, check_mesh)
{
  auto& mesh = this->getMesh();

  const int NE = mesh.GetNE();
  SLIC_INFO(axom::fmt::format("The mesh has {} elements", NE));
  EXPECT_GT(NE, 0);

  const int NV = mesh.GetNV();
  SLIC_INFO(axom::fmt::format("The mesh has {} vertices", NV));
  EXPECT_GT(NV, 0);

  const auto bbox = this->meshBoundingBox();
  SLIC_INFO(axom::fmt::format("The mesh bounding box is: {}", bbox));
}

//-----------------------------------------------------------------------------

TEST_F(SamplingShaperTest, basic_circle)
{
  const auto& testname =
    ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2

shapes:
- name: circle_shape
  material: {}
  geometry:
    format: c2c
    path: {}
)";

  const std::string circle_material = "circleMat";

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname),
                                   unit_circle_contour);

  ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                 axom::fmt::format(shape_template,
                                                   circle_material,
                                                   contour_file.getFileName()));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getFileName());
  this->runShaping(shape_file.getFileName());

  // check that the result has a volume fraction field associated with the circle material
  constexpr double expected_volume = M_PI;
  this->checkExpectedVolumeFractions(circle_material, expected_volume);

  // Save meshes and fields
  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest, disk_via_replacement)
{
  const auto& testname =
    ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2
units: cm

shapes:
- name: circle_outer
  material: outer
  geometry:
    format: c2c
    path: {0}
    units: cm
- name: void_inner
  material: inner
  geometry:
    format: c2c
    path: {0}
    units: cm
    operators:
      - scale: .5
)";

  const std::string outer_material = "outer";
  const std::string inner_material = "inner";

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname),
                                   unit_circle_contour);

  ScopedTemporaryFile shape_file(
    axom::fmt::format("{}.yaml", testname),
    axom::fmt::format(shape_template, contour_file.getFileName()));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  this->validateShapeFile(shape_file.getFileName());

  this->runShaping(shape_file.getFileName());

  // check that the result has a volume fraction field associated with the circle material
  constexpr double expected_inner_area = .5 * .5 * M_PI;
  constexpr double expected_outer_area = M_PI - expected_inner_area;
  this->checkExpectedVolumeFractions(outer_material, expected_outer_area);
  this->checkExpectedVolumeFractions(inner_material, expected_inner_area);

  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest, disk_via_replacement_with_background)
{
  using Point2D = typename SamplingShaperTest::Point2D;

  const auto& testname =
    ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2
units: cm

shapes:
- name: background
  material: {1}
  geometry:
    format: none
- name: circle_outer
  material: {2}
  geometry:
    format: c2c
    path: {0}
    units: cm
- name: void_inner
  material: {3}
  geometry:
    format: c2c
    path: {0}
    units: cm
    operators:
      - scale: .5
)";

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname),
                                   unit_circle_contour);

  // Set background material to 'void' (which is not present elsewhere)
  {
    ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                   axom::fmt::format(shape_template,
                                                     contour_file.getFileName(),
                                                     "void",
                                                     "disk",
                                                     "hole"));

    // Create an initial background material set to 1 everywhere
    std::map<std::string, mfem::GridFunction*> initialGridFunctions;
    {
      auto* vf = this->registerVolFracGridFunction("init_vf_bg");
      this->initializeVolFracGridFunction(
        vf,
        [](int, const Point2D&, int) -> double { return 1.; });
      initialGridFunctions["void"] = vf;
    }

    this->validateShapeFile(shape_file.getFileName());
    this->runShaping(shape_file.getFileName(), initialGridFunctions);

    // check that the result has a volume fraction field associated with the circle material
    constexpr double expected_hole_area = .5 * .5 * M_PI;
    constexpr double expected_disk_area = M_PI - expected_hole_area;
    const auto range = this->meshBoundingBox().range();
    const double expected_bg_area = range[0] * range[1] - M_PI;

    this->checkExpectedVolumeFractions("disk", expected_disk_area);
    this->checkExpectedVolumeFractions("hole", expected_hole_area);
    this->checkExpectedVolumeFractions("void", expected_bg_area);
  }

  // clean up data collection
  for(const auto& name :
      {"vol_frac_void", "vol_frac_hole", "vol_frac_disk", "init_vf_bg"})
  {
    this->getDC().DeregisterField(name);
  }

  // Set background and inner hole materials to 'void'
  {
    ScopedTemporaryFile shape_file(axom::fmt::format("{}.yaml", testname),
                                   axom::fmt::format(shape_template,
                                                     contour_file.getFileName(),
                                                     "void",
                                                     "disk",
                                                     "void"));

    // Create an initial background material set to 1 everywhere
    std::map<std::string, mfem::GridFunction*> initialGridFunctions;
    {
      auto* vf = this->registerVolFracGridFunction("init_vf_bg");
      this->initializeVolFracGridFunction(
        vf,
        [](int, const Point2D&, int) -> double { return 1.; });
      initialGridFunctions["void"] = vf;
    }

    this->validateShapeFile(shape_file.getFileName());
    this->runShaping(shape_file.getFileName(), initialGridFunctions);

    // check that the result has a volume fraction field associated with the circle material
    constexpr double expected_disk_area = M_PI - .5 * .5 * M_PI;
    const auto range = this->meshBoundingBox().range();
    const double expected_void_area = range[0] * range[1] - expected_disk_area;

    this->checkExpectedVolumeFractions("disk", expected_disk_area);
    this->checkExpectedVolumeFractions("void", expected_void_area);
  }

  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest, preshaped_materials)
{
  using Point2D = typename SamplingShaperTest::Point2D;

  const std::string& testname =
    ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2
units: cm

shapes:
- name: background
  material: {0}
  geometry:
    format: none
# Left replaces background void
- name: left_side
  material: {1}
  geometry:
    format: none
# Odd cells replace background void, but not left
- name: odd_cells
  material: {2}
  geometry:
    format: none
  does_not_replace: [{1}]
)";

  ScopedTemporaryFile shape_file(
    axom::fmt::format("{}.yaml", testname),
    axom::fmt::format(shape_template, "void", "left", "odds"));

  if(very_verbose_output)
  {
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  // Create an initial background material set to 1 everywhere
  std::map<std::string, mfem::GridFunction*> initialGridFunctions;
  {
    auto* vf = this->registerVolFracGridFunction("init_vf_bg");
    this->initializeVolFracGridFunction(
      vf,
      [](int, const Point2D&, int) -> double { return 1.; });
    initialGridFunctions["void"] = vf;

    // Note: element attributes were set earlier based on quadrant of cell's centroid (1, 2, 3 and 4)
    vf = this->registerVolFracGridFunction("init_vf_left");
    this->initializeVolFracGridFunction(
      vf,
      [](int, const Point2D&, int attr) -> double {
        return (attr == 3 || attr == 4) ? 1. : 0;
      });
    initialGridFunctions["left"] = vf;

    vf = this->registerVolFracGridFunction("init_vf_odds");
    this->initializeVolFracGridFunction(
      vf,
      [](int idx, const Point2D&, int) -> double {
        return idx % 2 == 1 ? 1. : 0;
      });
    initialGridFunctions["odds"] = vf;
  }

  this->validateShapeFile(shape_file.getFileName());
  this->runShaping(shape_file.getFileName(), initialGridFunctions);

  // check that the result has a volume fraction field associated with the circle material
  const auto range = this->meshBoundingBox().range();

  // Left covers half the mesh and is not replaced
  const double expected_left_area = range[0] * range[1] / 2.;
  // Odds should cover half of the right side of the mesh
  const double expected_odds_area = range[0] * range[1] / 4.;
  // The rest should be void
  const double expected_void_area =
    range[0] * range[1] - expected_left_area - expected_odds_area;

  this->checkExpectedVolumeFractions("left", expected_left_area);
  this->checkExpectedVolumeFractions("odds", expected_odds_area);
  this->checkExpectedVolumeFractions("void", expected_void_area);

  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

TEST_F(SamplingShaperTest, disk_with_multiple_preshaped_materials)
{
  using Point2D = typename SamplingShaperTest::Point2D;

  const std::string& testname =
    ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string shape_template = R"(
dimensions: 2
units: cm

shapes:
- name: background
  material: void
  geometry:
    format: none
- name: circle_outer
  material: disk
  geometry:
    format: c2c
    path: {0}
    units: cm
- name: circle_inner
  material: hole
  geometry:
    format: c2c
    path: {0}
    units: cm
    operators:
      - scale: .5
  replaces: [void, disk]
- name: left_side
  material: left
  geometry:
    format: none
  does_not_replace: [disk]
- name: odd_cells
  material: odds
  geometry:
    format: none
  replaces: [void, hole]
)";

  ScopedTemporaryFile contour_file(axom::fmt::format("{}.contour", testname),
                                   unit_circle_contour);

  ScopedTemporaryFile shape_file(
    axom::fmt::format("{}.yaml", testname),
    axom::fmt::format(shape_template, contour_file.getFileName()));

  if(very_verbose_output)
  {
    SLIC_INFO("Contour file: \n" << contour_file.getFileContents());
    SLIC_INFO("Shape file: \n" << shape_file.getFileContents());
  }

  // Create an initial background material set to 1 everywhere
  std::map<std::string, mfem::GridFunction*> initialGridFunctions;
  {
    // initial background void material is set everywhere
    auto* vf = this->registerVolFracGridFunction("init_vf_bg");
    this->initializeVolFracGridFunction(
      vf,
      [](int, const Point2D&, int) -> double { return 1.; });
    initialGridFunctions["void"] = vf;

    // initial left material is set based on mesh attributes
    // Note: element attributes were set earlier based on quadrant of cell's centroid (1, 2, 3 and 4)
    vf = this->registerVolFracGridFunction("init_vf_left");
    this->initializeVolFracGridFunction(
      vf,
      [](int, const Point2D&, int attr) -> double {
        return (attr == 3 || attr == 4) ? 1. : 0;
      });
    initialGridFunctions["left"] = vf;

    // initial "odds" material is based on the parity of the element indices
    vf = this->registerVolFracGridFunction("init_vf_odds");
    this->initializeVolFracGridFunction(
      vf,
      [](int idx, const Point2D&, int) -> double {
        return idx % 2 == 1 ? 1. : 0;
      });
    initialGridFunctions["odds"] = vf;
  }

  // For this example, we keep the full disk between radii 1 and .5
  // The 'left' material is set for all cells to the left of the y-axis
  // but does not replace the disk material
  // The interior hole is within radius .5, but is replaced by left and odds
  // The odds material is all cells w/ odd index, but not covering 'left' or 'disk'
  // The end result for the void background is everything that's left
  this->validateShapeFile(shape_file.getFileName());
  this->runShaping(shape_file.getFileName(), initialGridFunctions);

  // check that the result has the correct volume fractions
  const auto range = this->meshBoundingBox().range();
  const auto total_area = range[0] * range[1];
  const auto left_orig = total_area / 2;
  constexpr auto hole_orig = .5 * .5 * M_PI;

  constexpr double expected_disk_area = M_PI - hole_orig;
  const double expected_left_area = left_orig - expected_disk_area / 2;
  const double expected_hole_area = 0.19683;  // from program output
  const double expected_odds_area = 3.41180;  // from program output
  const double expected_void_area = 3.21497;  // from program output

  this->checkExpectedVolumeFractions("left", expected_left_area);
  this->checkExpectedVolumeFractions("disk", expected_disk_area);
  this->checkExpectedVolumeFractions("odds", expected_odds_area);
  this->checkExpectedVolumeFractions("hole", expected_hole_area);
  this->checkExpectedVolumeFractions("void", expected_void_area);

  if(very_verbose_output)
  {
    this->getDC().Save(testname, axom::sidre::Group::getDefaultIOProtocol());
  }
}

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;
#ifdef AXOM_USE_MPI
  // This is needed because of Axom's c2c reader.
  MPI_Init(&argc, &argv);
  // See if this aborts right away.
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
#endif

  ::testing::InitGoogleTest(&argc, argv);
  slic::SimpleLogger logger(slic::message::Info);

  result = RUN_ALL_TESTS();
#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif
  return result;
}

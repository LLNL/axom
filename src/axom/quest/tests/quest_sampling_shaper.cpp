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

    // Offset the mesh to lie w/in the bounding box
    for(int i = 0; i < mesh->GetNV(); ++i)
    {
      double* v = mesh->GetVertex(i);
      for(int d = 0; d < dim; ++d)
      {
        v[d] += lo[d];
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

  // parse and validate the Klee shapefile; fail the test if invalid
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

  void runShaping(const std::string& shapefile,
                  const std::map<std::string, mfem::GridFunction*>& init_vf_map = {})
  {
    SLIC_INFO(axom::fmt::format("Reading shape set from {}", shapefile));
    klee::ShapeSet shapeSet(klee::readShapeSet(shapefile));

    SLIC_INFO(axom::fmt::format("Shaping materials..."));
    quest::SamplingShaper shaper(shapeSet, &m_dc);
    shaper.setVerbosity(true);

    if(!init_vf_map.empty())
    {
      shaper.importInitialVolumeFractions(init_vf_map);
    }
    shaper.printRegisteredFieldNames("*** After importing volume fractions");

    const auto shapeDim = shapeSet.getDimensions();
    for(const auto& shape : shapeSet.getShapes())
    {
      SLIC_INFO(axom::fmt::format("\tshape {} -> material {}",
                                  shape.getName(),
                                  shape.getMaterial()));

      shaper.loadShape(shape);
      shaper.prepareShapeQuery(shapeDim, shape);
      shaper.runShapeQuery(shape);
      shaper.applyReplacementRules(shape);
      shaper.finalizeShapeQuery();
    }

    shaper.adjustVolumeFractions();

    shaper.printRegisteredFieldNames("*** After shaping volume fractions");
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
  this->getDC().Save();
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

  // Save meshes and fields
  this->getDC().Save();
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
  for(const auto& name : {"void", "hole", "disk", "init_vf_bg"})
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
    ;
    const auto range = this->meshBoundingBox().range();
    const double expected_void_area = range[0] * range[1] - expected_disk_area;

    this->checkExpectedVolumeFractions("disk", expected_disk_area);
    this->checkExpectedVolumeFractions("void", expected_void_area);
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

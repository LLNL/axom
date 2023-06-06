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
namespace slic = axom::slic;
namespace sidre = axom::sidre;
namespace quest = axom::quest;

namespace
{
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

  void runShaping(const std::string& shapefile)
  {
    SLIC_INFO(axom::fmt::format("Reading shape set from {}", shapefile));
    klee::ShapeSet shapeSet(klee::readShapeSet(shapefile));

    SLIC_INFO(axom::fmt::format("Shaping materials..."));
    quest::SamplingShaper shaper(shapeSet, &m_dc);
    shaper.setVerbosity(true);

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

TEST_F(SamplingShaperTest, check_mesh)
{
  auto& mesh = this->getMesh();

  const int NE = mesh.GetNE();
  SLIC_INFO(axom::fmt::format("The mesh has {} elements", NE));
  EXPECT_GT(NE, 0);

  const int NV = mesh.GetNV();
  SLIC_INFO(axom::fmt::format("The mesh has {} vertices", NV));
  EXPECT_GT(NV, 0);

  mfem::Vector mn, mx;
  mesh.GetBoundingBox(mn, mx);
  SLIC_INFO(
    axom::fmt::format("The mesh bounding box is:  {{min: {}, {}, max: {},{}}}",
                      mn[0],
                      mn[1],
                      mx[0],
                      mx[1]));
}

//-----------------------------------------------------------------------------

TEST_F(SamplingShaperTest, basic_circle)
{
  const auto& testname =
    ::testing::UnitTest::GetInstance()->current_test_info()->name();

  const std::string contour_contents =
    "piece = circle(origin=(0cm, 0cm), radius=1cm, start=0deg, end=360deg)";

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
                                   contour_contents);

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
  auto& dc = this->getDC();
  const std::string vf_name = axom::fmt::format("vol_frac_{}", circle_material);
  EXPECT_TRUE(dc.HasField(vf_name))
    << "Did not have material '" << vf_name << "'";

  // check shaped-in volume fraction of unit circle -- should be around PI
  SLIC_INFO("Shaped volume fraction of unit circle is "
            << this->gridFunctionVolume(vf_name));

  EXPECT_NEAR(M_PI, this->gridFunctionVolume(vf_name), 1e-2);

  //---------------------------------------------------------------------------
  // Save meshes and fields
  //---------------------------------------------------------------------------
  dc.Save();
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

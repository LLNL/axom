// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom utils
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/slic.hpp"

// Quest includes
#include "axom/quest/interface/signed_distance.hpp"
#include "quest_test_utilities.hpp"

// gtest
#include "gtest/gtest.h"

// Note: Disable death tests when MPI is enabled since they're not being properly handled
#ifdef AXOM_USE_MPI
  #define _DEATH_TESTS_LOCALLY_DISABLED 1
#else
  #define _DEATH_TESTS_LOCALLY_DISABLED 0
const char IGNORE_OUTPUT[] = ".*";
#endif

// C/C++ includes
#include <fstream>
#include <sstream>

// Aliases
namespace quest = axom::quest;
namespace mint = axom::mint;
namespace primal = axom::primal;
namespace utilities = axom::utilities;

using UnstructuredMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;

//#define WRITE_VTK_OUTPUT 1
#define REMOVE_FILES 1

//------------------------------------------------------------------------------
//  HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
#if defined(AXOM_USE_MPI) && defined(AXOM_USE_MPI3)
constexpr bool USE_MPI3_SHARED_MEMORY = true;
#endif

/*!
 * \brief Generate a mesh of 4 triangles along the XY plane.
 *
 * \param [in] file the file to write to.
 *
 * \pre file.empty() == false.
 * \note Used primarily for debugging.
 */
void generate_planar_mesh_stl_file(const std::string& file)
{
  EXPECT_FALSE(file.empty());

  std::ofstream ofs(file.c_str());
  EXPECT_TRUE(ofs.is_open());

  ofs << "solid plane" << std::endl;

  // Triangle T1
  ofs << "\t facet normal 0.0 0.0 1.0" << std::endl;
  ofs << "\t\t outer loop" << std::endl;
  ofs << "\t\t\t vertex -5.0 -5.0 0.0" << std::endl;
  ofs << "\t\t\t vertex 5.0 -5.0 0.0" << std::endl;
  ofs << "\t\t\t vertex 0.0  0.0 0.0" << std::endl;
  ofs << "\t\t endloop" << std::endl;
  ofs << "\t endfacet" << std::endl;

  // Triangle T2
  ofs << "\t facet normal 0.0 0.0 1.0" << std::endl;
  ofs << "\t\t outer loop" << std::endl;
  ofs << "\t\t\t vertex 5.0 -5.0 0.0" << std::endl;
  ofs << "\t\t\t vertex 5.0  5.0 0.0" << std::endl;
  ofs << "\t\t\t vertex 0.0  0.0 0.0" << std::endl;
  ofs << "\t\t endloop" << std::endl;
  ofs << "\t endfacet" << std::endl;

  // Triangle T3
  ofs << "\t facet normal 0.0 0.0 1.0" << std::endl;
  ofs << "\t\t outer loop" << std::endl;
  ofs << "\t\t\t vertex  5.0 5.0 0.0" << std::endl;
  ofs << "\t\t\t vertex -5.0 5.0 0.0" << std::endl;
  ofs << "\t\t\t vertex  0.0 0.0 0.0" << std::endl;
  ofs << "\t\t endloop" << std::endl;
  ofs << "\t endfacet" << std::endl;

  // Triangle T4
  ofs << "\t facet normal 0.0 0.0 1.0" << std::endl;
  ofs << "\t\t outer loop" << std::endl;
  ofs << "\t\t\t vertex -5.0  5.0 0.0" << std::endl;
  ofs << "\t\t\t vertex -5.0 -5.0 0.0" << std::endl;
  ofs << "\t\t\t vertex  0.0  0.0 0.0" << std::endl;
  ofs << "\t\t endloop" << std::endl;
  ofs << "\t endfacet" << std::endl;

  ofs << "endsolid" << std::endl;
  ofs.close();
}

/*!
 * \brief Generates a simple ASCII STL file consisting of a single triangle.
 *
 * \param [in] file the file to write to.
 * \pre file.empty() == false.
 *
 * \note Used primarily for debugging.
 */
void generate_stl_file(const std::string& file)
{
  EXPECT_FALSE(file.empty());

  std::ofstream ofs(file.c_str());
  EXPECT_TRUE(ofs.is_open());

  ofs << "solid triangle" << std::endl;
  ofs << "\t facet normal 0.0 0.0 1.0" << std::endl;
  ofs << "\t\t outer loop" << std::endl;
  ofs << "\t\t\t vertex 0.0 0.0 0.0" << std::endl;
  ofs << "\t\t\t vertex 1.0 0.0 0.0" << std::endl;
  ofs << "\t\t\t vertex 0.0 1.0 0.0" << std::endl;
  ofs << "\t\t endloop" << std::endl;
  ofs << "\t endfacet" << std::endl;
  ofs << "endsolid triangle" << std::endl;

  ofs.close();
}

/*!
 * \brief Returns the bounding box of the mesh.
 * \param [in] mesh pointer to the mesh instance.
 * \return bb bounding box of the mesh
 */
primal::BoundingBox<double, 3> getBounds(const axom::mint::Mesh* mesh)
{
  SLIC_ASSERT(mesh != nullptr);

  primal::BoundingBox<double, 3> bb;
  primal::Point<double, 3> pt;

  const int nnodes = mesh->getNumberOfNodes();
  for(int inode = 0; inode < nnodes; ++inode)
  {
    mesh->getNode(inode, pt.data());
    bb.addPoint(pt);
  }

  return (bb);
}

/*!
 * \brief Generates a uniform mesh surrounding the given triangle mesh.
 * \param [in] mesh pointer to the input mesh.
 * \param [in] umesh pointer to the uniform mesh;
 */
void getUniformMesh(const UnstructuredMesh* mesh, mint::UniformMesh*& umesh)
{
  SLIC_ASSERT(mesh != nullptr);
  SLIC_ASSERT(umesh == nullptr);

  constexpr int N = 16;  // number of points along each dimension

  primal::BoundingBox<double, 3> bb = getBounds(mesh);
  bb.expand(2.0);

  const double* lo = bb.getMin().data();
  const double* hi = bb.getMax().data();

  // construct an N x N x N grid
  umesh = new mint::UniformMesh(lo, hi, N, N, N);
}

/*!
 * \brief Tests the signed against an anlytic surface mesh input.
 * \param [in] use_shared indicates whether to use shared memory or not.
 */
void check_analytic_plane(bool use_shared = false)
{
  // STEP 0: construct uniform box mesh
  constexpr int NDIMS = 3;
  constexpr int N = 16;
  double lo[3] = {-4.0, -4.0, -4.0};
  double hi[3] = {4.0, 4.0, 4.0};
  mint::UniformMesh mesh(lo, hi, N, N, N);
  double* phi = mesh.createField<double>("phi", mint::NODE_CENTERED);
  double* err = mesh.createField<double>("err", mint::NODE_CENTERED);

  // STEP 1: generate planar STL mesh file
  const std::string file = "plane.stl";
  generate_planar_mesh_stl_file(file);

  // STEP 2: define analytic plane corresponding to the planar mesh;
  const primal::Point<double, 3> origin {0.0, 0.0, 0.0};
  const primal::Vector<double, 3> normal {0.0, 0.0, 1.0};
  primal::Plane<double, 3> analytic_plane(normal, origin);

  // STEP 2: initialize the signed distance
  quest::signed_distance_use_shared_memory(use_shared);
  quest::signed_distance_set_closed_surface(false);
  quest::signed_distance_init(file);
  EXPECT_TRUE(quest::signed_distance_initialized());

  axom::IndexType nnodes = mesh.getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < nnodes; ++inode)
  {
    double pt[NDIMS];
    mesh.getNode(inode, pt);
    primal::Point<double, NDIMS> thepoint(pt);
    phi[inode] = quest::signed_distance_evaluate(pt[0], pt[1], pt[2]);

    const double phi_expected = analytic_plane.signedDistance(thepoint);
    EXPECT_DOUBLE_EQ(phi[inode], phi_expected);
    err[inode] = utilities::abs(phi[inode] - phi_expected);
  }

#ifdef WRITE_VTK_OUTPUT
  std::ostringstream oss;
  oss << "analytic_plane_mesh_";
  oss << ((use_shared) ? "shared_test.vtk" : "test.vtk");
  mint::write_vtk(&mesh, oss.str());
#endif

  quest::signed_distance_finalize();
  EXPECT_FALSE(quest::signed_distance_initialized());

#ifdef REMOVE_FILES
  axom::utilities::filesystem::removeFile(file);
#endif
}

}  // namespace

//------------------------------------------------------------------------------
//  UNIT TESTS
//------------------------------------------------------------------------------
TEST(quest_signed_distance_interface_DeathTest, get_mesh_bounds_invalid_calls)
{
  EXPECT_FALSE(quest::signed_distance_initialized());

#if !_DEATH_TESTS_LOCALLY_DISABLED
  double lo[3];
  double hi[3];
  EXPECT_DEATH_IF_SUPPORTED(quest::signed_distance_get_mesh_bounds(lo, hi),
                            IGNORE_OUTPUT);
#endif

  constexpr int NDIMS = 3;
  constexpr double SPHERE_RADIUS = 0.5;
  constexpr int SPHERE_THETA_RES = 25;
  constexpr int SPHERE_PHI_RES = 25;
  const double SPHERE_CENTER[3] = {0.0, 0.0, 0.0};

  UnstructuredMesh* surface_mesh = new UnstructuredMesh(NDIMS, mint::TRIANGLE);
  quest::utilities::getSphereSurfaceMesh(surface_mesh,
                                         SPHERE_CENTER,
                                         SPHERE_RADIUS,
                                         SPHERE_THETA_RES,
                                         SPHERE_PHI_RES);

  quest::signed_distance_init(surface_mesh);

#if !_DEATH_TESTS_LOCALLY_DISABLED
  EXPECT_DEATH_IF_SUPPORTED(quest::signed_distance_get_mesh_bounds(nullptr, hi),
                            IGNORE_OUTPUT);

  EXPECT_DEATH_IF_SUPPORTED(quest::signed_distance_get_mesh_bounds(lo, nullptr),
                            IGNORE_OUTPUT);
#endif

  quest::signed_distance_finalize();
  delete surface_mesh;
  surface_mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_signed_distance_interface_DeathTest, call_evaluate_before_init)
{
  EXPECT_FALSE(quest::signed_distance_initialized());

#if !_DEATH_TESTS_LOCALLY_DISABLED
  EXPECT_DEATH_IF_SUPPORTED(quest::signed_distance_evaluate(0.0, 0.0),
                            IGNORE_OUTPUT);
  double x[2];
  double y[2];
  double z[2];
  double phi[2];

  EXPECT_DEATH_IF_SUPPORTED(quest::signed_distance_evaluate(x, y, z, 2, phi),
                            IGNORE_OUTPUT);
#endif
}

//------------------------------------------------------------------------------
TEST(quest_signed_distance_interface_DeathTest, set_params_after_init)
{
  EXPECT_FALSE(quest::signed_distance_initialized());

  constexpr int NDIMS = 3;
  constexpr double SPHERE_RADIUS = 0.5;
  constexpr int SPHERE_THETA_RES = 25;
  constexpr int SPHERE_PHI_RES = 25;
  const double SPHERE_CENTER[3] = {0.0, 0.0, 0.0};

  // STEP 0: initialiaze quest
  UnstructuredMesh* surface_mesh = new UnstructuredMesh(NDIMS, mint::TRIANGLE);
  quest::utilities::getSphereSurfaceMesh(surface_mesh,
                                         SPHERE_CENTER,
                                         SPHERE_RADIUS,
                                         SPHERE_THETA_RES,
                                         SPHERE_PHI_RES);

  quest::signed_distance_init(surface_mesh);

  // STEP 1: setting parameters after init() should fail
#if !_DEATH_TESTS_LOCALLY_DISABLED

  EXPECT_DEATH_IF_SUPPORTED(quest::signed_distance_set_dimension(3),
                            IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(quest::signed_distance_set_closed_surface(true),
                            IGNORE_OUTPUT);

  EXPECT_DEATH_IF_SUPPORTED(quest::signed_distance_set_verbose(true),
                            IGNORE_OUTPUT);

  EXPECT_DEATH_IF_SUPPORTED(
    quest::signed_distance_set_execution_space(quest::SignedDistExec::CPU),
    IGNORE_OUTPUT);
#endif

  // STEP 2: finalize
  quest::signed_distance_finalize();
  delete surface_mesh;
  surface_mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_signed_distance_interface, initialize)
{
  int rc = -1;

  // initializing with an empty file should return a non-zero value
  rc = quest::signed_distance_init("");
  EXPECT_TRUE(rc != 0);
  EXPECT_FALSE(quest::signed_distance_initialized());

  // initializing with a file that does not exist, should fail.
  rc = quest::signed_distance_init("badfile.stl");
  EXPECT_TRUE(rc != 0);
  EXPECT_FALSE(quest::signed_distance_initialized());

  // generate a temp STL file
  const std::string fileName = "test_triangle.stl";
  generate_stl_file(fileName);

  // test valid initialization
  quest::signed_distance_set_dimension(3);
  rc = quest::signed_distance_init(fileName);
  EXPECT_EQ(rc, 0);
  EXPECT_TRUE(quest::signed_distance_initialized());

  // finalize
  quest::signed_distance_finalize();
  EXPECT_FALSE(quest::signed_distance_initialized());

  // remove temp STL file
#ifdef REMOVE_FILES
  axom::utilities::filesystem::removeFile(fileName);
#endif
}

//------------------------------------------------------------------------------
TEST(quest_signed_distance_interface, get_mesh_bounds)
{
  constexpr int NDIMS = 3;
  constexpr double SPHERE_RADIUS = 0.5;
  constexpr int SPHERE_THETA_RES = 25;
  constexpr int SPHERE_PHI_RES = 25;
  const double SPHERE_CENTER[3] = {0.0, 0.0, 0.0};
  constexpr double EXPECTED_LO = -0.5;
  constexpr double EXPECTED_HI = 0.5;

  double lo[NDIMS];
  double hi[NDIMS];

  UnstructuredMesh* surface_mesh = new UnstructuredMesh(NDIMS, mint::TRIANGLE);
  quest::utilities::getSphereSurfaceMesh(surface_mesh,
                                         SPHERE_CENTER,
                                         SPHERE_RADIUS,
                                         SPHERE_THETA_RES,
                                         SPHERE_PHI_RES);

  quest::signed_distance_init(surface_mesh);

  quest::signed_distance_get_mesh_bounds(lo, hi);

  for(int i = 0; i < NDIMS; ++i)
  {
    EXPECT_DOUBLE_EQ(lo[i], EXPECTED_LO);
    EXPECT_DOUBLE_EQ(hi[i], EXPECTED_HI);
  }

  quest::signed_distance_finalize();

  delete surface_mesh;
  surface_mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_signed_distance_interface, analytic_plane)
{
  check_analytic_plane();

#if defined(AXOM_USE_MPI) && defined(AXOM_USE_MPI3)
  check_analytic_plane(USE_MPI3_SHARED_MEMORY);
#endif
}

//------------------------------------------------------------------------------
#if defined(AXOM_USE_MPI) && defined(AXOM_USE_MPI3)
TEST(quest_signed_distance_interface, call_twice_using_shared_memory)
{
  check_analytic_plane(USE_MPI3_SHARED_MEMORY);
  check_analytic_plane(USE_MPI3_SHARED_MEMORY);
}
#endif

//------------------------------------------------------------------------------
TEST(quest_signed_distance_interface, analytic_sphere)
{
  // STEP 0: constants
  constexpr int NDIMS = 3;
  constexpr double l1norm_expected = 6.7051997372579715;
  constexpr double l2norm_expected = 2.5894400431865519;
  constexpr double linf_expected = 0.00532092;
  constexpr double TOL = 1.e-3;

  constexpr double SPHERE_RADIUS = 0.5;
  constexpr int SPHERE_THETA_RES = 25;
  constexpr int SPHERE_PHI_RES = 25;
  const double SPHERE_CENTER[3] = {0.0, 0.0, 0.0};

  // STEP 1: create analytic sphere object to compare results with
  primal::Sphere<double, 3> analytic_sphere(SPHERE_RADIUS);

  // STEP 2: generate sphere mesh to pass to the signed distance query
  UnstructuredMesh* surface_mesh = new UnstructuredMesh(NDIMS, mint::TRIANGLE);
  quest::utilities::getSphereSurfaceMesh(surface_mesh,
                                         SPHERE_CENTER,
                                         SPHERE_RADIUS,
                                         SPHERE_THETA_RES,
                                         SPHERE_PHI_RES);

  // STEP 1: generate the uniform mesh where the signed distance will be
  //         computed
  mint::UniformMesh* umesh = nullptr;
  getUniformMesh(surface_mesh, umesh);

  // STEP 2: create node-centered fields on the uniform mesh to store
  //         the signed distance etc.
  double* phi_computed =
    umesh->createField<double>("phi_computed", mint::NODE_CENTERED);
  double* phi_expected =
    umesh->createField<double>("phi_expected", mint::NODE_CENTERED);
  double* phi_diff = umesh->createField<double>("phi_diff", mint::NODE_CENTERED);
  double* phi_err = umesh->createField<double>("phi_err", mint::NODE_CENTERED);

  // STEP 3: initialize the signed distance query
  quest::signed_distance_set_closed_surface(true);
  quest::signed_distance_init(surface_mesh);
  EXPECT_TRUE(quest::signed_distance_initialized());

  // STEP 4: Compute signed distance
  double l1norm = 0.0;
  double l2norm = 0.0;
  double linf = std::numeric_limits<double>::min();
  axom::IndexType nnodes = umesh->getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < nnodes; ++inode)
  {
    primal::Point<double, NDIMS> pt;
    umesh->getNode(inode, pt.data());

    phi_computed[inode] = quest::signed_distance_evaluate(pt[0], pt[1], pt[2]);
    phi_expected[inode] = analytic_sphere.computeSignedDistance(pt);
    EXPECT_NEAR(phi_computed[inode], phi_expected[inode], 1.e-2);

    // compute error
    phi_diff[inode] = phi_computed[inode] - phi_expected[inode];
    phi_err[inode] = utilities::abs(phi_diff[inode]);

    // update norms
    l1norm += phi_err[inode];
    l2norm += phi_diff[inode];
    linf = (phi_err[inode] > linf) ? phi_err[inode] : linf;

  }  // END for all nodes

  l2norm = std::sqrt(l2norm);

  EXPECT_NEAR(l1norm_expected, l1norm, TOL);
  EXPECT_NEAR(l2norm_expected, l2norm, TOL);
  EXPECT_NEAR(linf_expected, linf, TOL);

  // STEP 5: finalize the signed distance query
  quest::signed_distance_finalize();
  EXPECT_FALSE(quest::signed_distance_initialized());

#ifdef WRITE_VTK_OUTPUT
  mint::write_vtk(umesh, "analytic_sphere_test_mesh.vtk");
  mint::write_vtk(surface_mesh, "input_sphere_mesh.vtk");
#endif

  // STEP 6: delete mesh objects
  delete umesh;
  delete surface_mesh;

  umesh = nullptr;
  surface_mesh = nullptr;
}

//------------------------------------------------------------------------------
TEST(quest_signed_distance_interface, analytic_sphere_with_closest_pt_and_normal)
{
  // STEP 0: constants
  constexpr int NDIMS = 3;
  constexpr double l1norm_expected = 6.7051997372579715;
  constexpr double l2norm_expected = 2.5894400431865519;
  constexpr double linf_expected = 0.00532092;
  constexpr double TOL = 1.e-3;

  constexpr double SPHERE_RADIUS = 0.5;
  constexpr int SPHERE_THETA_RES = 25;
  constexpr int SPHERE_PHI_RES = 25;
  const double SPHERE_CENTER[3] = {0.0, 0.0, 0.0};

  using PointType = primal::Point<double, 3>;
  using VectorType = primal::Vector<double, 3>;
  using SphereType = primal::Sphere<double, 3>;

  // STEP 1: create analytic sphere object to compare results with
  SphereType analytic_sphere(SPHERE_RADIUS);

  // STEP 2: generate sphere mesh to pass to the signed distance query
  UnstructuredMesh* surface_mesh = new UnstructuredMesh(NDIMS, mint::TRIANGLE);
  quest::utilities::getSphereSurfaceMesh(surface_mesh,
                                         SPHERE_CENTER,
                                         SPHERE_RADIUS,
                                         SPHERE_THETA_RES,
                                         SPHERE_PHI_RES);

  // STEP 1: generate the uniform mesh where the signed distance will be
  //         computed
  mint::UniformMesh* umesh = nullptr;
  getUniformMesh(surface_mesh, umesh);

  // STEP 2: create node-centered fields on the uniform mesh to store
  //         the signed distance etc.
  double* phi_computed =
    umesh->createField<double>("phi_computed", mint::NODE_CENTERED);
  double* phi_expected =
    umesh->createField<double>("phi_expected", mint::NODE_CENTERED);
  double* phi_diff = umesh->createField<double>("phi_diff", mint::NODE_CENTERED);
  double* phi_err = umesh->createField<double>("phi_err", mint::NODE_CENTERED);

  // STEP 3: initialize the signed distance query
  quest::signed_distance_set_closed_surface(true);
  quest::signed_distance_init(surface_mesh);
  EXPECT_TRUE(quest::signed_distance_initialized());

  // STEP 4: Compute signed distance
  double l1norm = 0.0;
  double l2norm = 0.0;
  double linf = std::numeric_limits<double>::min();
  axom::IndexType nnodes = umesh->getNumberOfNodes();
  for(axom::IndexType inode = 0; inode < nnodes; ++inode)
  {
    PointType pt;
    umesh->getNode(inode, pt.data());

    PointType closest_point;
    VectorType normal;

    phi_computed[inode] = quest::signed_distance_evaluate(pt[0],
                                                          pt[1],
                                                          pt[2],
                                                          closest_point[0],
                                                          closest_point[1],
                                                          closest_point[2],
                                                          normal[0],
                                                          normal[1],
                                                          normal[2]);

    phi_expected[inode] = analytic_sphere.computeSignedDistance(pt);
    EXPECT_NEAR(phi_computed[inode], phi_expected[inode], 1.e-2);

    // Check that the computed closest point on the discretized sphere is close
    // to where it would be on an analytic sphere
    const auto cp_expected = primal::closest_point(pt, analytic_sphere);
    EXPECT_NEAR(0., primal::squared_distance(cp_expected, closest_point), 1e-2);

    // Check that the computed (pseudo)-normal on the discretized sphere is close
    // to where it would be on an analytic sphere, i.e. the dot product
    // between the two should be close to 1
    const auto normal_expected = VectorType(cp_expected).unitVector();
    EXPECT_NEAR(1., normal.dot(normal_expected), 1e-2);

    // Check that the distance (squared) is the same as the distance
    // between the query point and the returned closest point
    EXPECT_NEAR(phi_computed[inode] * phi_computed[inode],
                primal::squared_distance(pt, closest_point),
                1e-8);

    // compute error
    phi_diff[inode] = phi_computed[inode] - phi_expected[inode];
    phi_err[inode] = utilities::abs(phi_diff[inode]);

    // update norms
    l1norm += phi_err[inode];
    l2norm += phi_diff[inode];
    linf = (phi_err[inode] > linf) ? phi_err[inode] : linf;

  }  // END for all nodes

  l2norm = std::sqrt(l2norm);

  EXPECT_NEAR(l1norm_expected, l1norm, TOL);
  EXPECT_NEAR(l2norm_expected, l2norm, TOL);
  EXPECT_NEAR(linf_expected, linf, TOL);

  // STEP 5: finalize the signed distance query
  quest::signed_distance_finalize();
  EXPECT_FALSE(quest::signed_distance_initialized());

#ifdef WRITE_VTK_OUTPUT
  mint::write_vtk(umesh, "analytic_sphere_test_mesh.vtk");
  mint::write_vtk(surface_mesh, "input_sphere_mesh.vtk");
#endif

  // STEP 6: delete mesh objects
  delete umesh;
  delete surface_mesh;

  umesh = nullptr;
  surface_mesh = nullptr;
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);
#endif

  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return result;
}

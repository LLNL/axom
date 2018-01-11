

// Axom includes
#include "mint/UnstructuredMesh.hpp"
#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"
#include "quest/quest.hpp"

#include "fmt/format.h"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

// Google test include
#include "gtest/gtest.h"

typedef axom::mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;

/**
 * \brief Generates a mint Unstructured mesh representing a tetrahedron.
 * \note Allocates a UnstructuredMesh instance, which must be deleted by the user
 */
TriangleMesh* createTriangleMesh()
{
  typedef axom::mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;

  TriangleMesh* surface_mesh = new TriangleMesh(3);
  surface_mesh->insertNode( -0.000003, -0.000003, 19.999999);
  surface_mesh->insertNode(-18.213671,  4.880339, -6.666668);
  surface_mesh->insertNode(  4.880339,-18.213671, -6.666668);
  surface_mesh->insertNode( 13.333334, 13.333334, -6.666663);
  int cell[3];
  cell[0] = 0;    cell[1] = 1;    cell[2] = 2;
  surface_mesh->insertCell(cell, MINT_TRIANGLE, 3);
  cell[0] = 0;    cell[1] = 3;    cell[2] = 1;
  surface_mesh->insertCell(cell, MINT_TRIANGLE, 3);
  cell[0] = 0;    cell[1] = 2;    cell[2] = 3;
  surface_mesh->insertCell(cell, MINT_TRIANGLE, 3);
  cell[0] = 1;    cell[1] = 3;    cell[2] = 2;
  surface_mesh->insertCell(cell, MINT_TRIANGLE, 3);

  return surface_mesh;
}

TEST( quest_interface, pointer_initialize )
{
  const int IGNORE = -1;

  SLIC_INFO(fmt::format("Initializing InOutOctree over triangle mesh ..."));

#ifdef AXOM_USE_MPI
  axom::quest::initialize(MPI_COMM_WORLD, createTriangleMesh(), false, 3, IGNORE, IGNORE);
#else
  axom::quest::initialize(createTriangleMesh(), false, 3, IGNORE, IGNORE);
#endif

  EXPECT_TRUE(axom::quest::inside(3, 2, 0));
  EXPECT_TRUE(axom::quest::inside(-1, 2, -1));
  EXPECT_FALSE(axom::quest::inside(4, 4, -7));

  axom::quest::finalize();

}

int main( int argc, char** argv )
{
#ifdef AXOM_USE_MPI
  // Initialize MPI
  MPI_Init( &argc, &argv );
#endif
  ::testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();

#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return result;
}




// Axom includes
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/core/UnitTestLogger.hpp"
#include "quest_test_utilities.hpp"
#include "axom/quest/interface/quest.hpp"

#include "fmt/fmt.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

// Google test include
#include "gtest/gtest.h"

typedef axom::mint::UnstructuredMesh< axom::mint::SINGLE_SHAPE > UMesh;

TEST( quest_interface, pointer_initialize )
{
  const int IGNORE_PARAM = -1;

  SLIC_INFO(fmt::format("Initializing InOutOctree over triangle mesh ..."));

  axom::mint::Mesh* input_mesh =
    axom::quest::utilities::make_tetrahedron_mesh();

#ifdef AXOM_USE_MPI
  axom::quest::initialize(MPI_COMM_WORLD, input_mesh, false, 3, IGNORE_PARAM,
                          IGNORE_PARAM);
#else
  axom::quest::initialize(input_mesh, false, 3, IGNORE_PARAM, IGNORE_PARAM);
#endif

  EXPECT_TRUE(axom::quest::inside(3, 2, 0));
  EXPECT_TRUE(axom::quest::inside(-1, 2, -1));
  EXPECT_FALSE(axom::quest::inside(4, 4, -7));

  axom::quest::finalize();
  delete input_mesh;

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

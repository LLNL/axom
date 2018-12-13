

// Axom includes
#include "axom/mint.hpp"
#include "quest_test_utilities.hpp"

#include "axom/quest/interface/inout_query.hpp"
#include "axom/quest/interface/signed_distance.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

// Google test include
#include "gtest/gtest.h"

namespace mint = axom::mint;
namespace quest = axom::quest;

// Test initializing quest inout from a preloaded mesh
TEST( quest_initialize, inout_pointer_initialize )
{
  int rc = quest::QUEST_INOUT_SUCCESS;

  mint::Mesh* input_mesh =
    axom::quest::utilities::make_tetrahedron_mesh();

  // Note: the following call updates the input_mesh pointer
#ifdef AXOM_USE_MPI
  rc = quest::inout_init(input_mesh, MPI_COMM_WORLD);
#else
  rc = quest::inout_init(input_mesh);
#endif
  EXPECT_EQ(quest::QUEST_INOUT_SUCCESS, rc);

  EXPECT_TRUE(quest::inout_initialized());

  EXPECT_TRUE(quest::inout_inside(3, 2, 0));
  EXPECT_TRUE(quest::inout_inside(-1, 2, -1));
  EXPECT_FALSE(quest::inout_inside(4, 4, -7));

  rc = quest::inout_finalize();
  EXPECT_EQ(quest::QUEST_INOUT_SUCCESS, rc);

  delete input_mesh;
}


// Test initializing quest signed_distance from a preloaded mesh
TEST( quest_initialize, signed_distance_pointer_initialize )
{
  int rc = 0;

  mint::Mesh* input_mesh =
    axom::quest::utilities::make_tetrahedron_mesh();

#ifdef AXOM_USE_MPI
  rc = quest::signed_distance_init(input_mesh, MPI_COMM_WORLD);
#else
  rc = quest::signed_distance_init(input_mesh);
#endif
  EXPECT_EQ(0, rc);

  EXPECT_TRUE(quest::signed_distance_initialized());

  EXPECT_GT(0., quest::signed_distance_evaluate(3, 2, 0));
  EXPECT_GT(0., quest::signed_distance_evaluate(-1, 2, -1));
  EXPECT_LT(0., quest::signed_distance_evaluate(4, 4, -7));

  quest::signed_distance_finalize();

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

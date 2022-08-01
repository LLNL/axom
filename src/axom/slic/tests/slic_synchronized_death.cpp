// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core/utilities/Utilities.hpp"

#include "axom/slic/interface/slic.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/slic/streams/SynchronizedStream.hpp"

#include "gtest/gtest.h"

#include <string>
#include <sstream>

#include "mpi.h"

// namespace alias
namespace slic = axom::slic;

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
void slic_error() { SLIC_ERROR("SLIC_ERROR test!"); }
}  // namespace

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
TEST(slic_synchronized_death, test_slic_error_death)
{
  MPI_Init(NULL, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // initialize slic
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Debug);

  std::string msgfmt = "[<LEVEL>]:;;<MESSAGE>;;\n@@<FILE>\n@@<LINE>";

  slic::addStreamToAllMsgLevels(
    new slic::SynchronizedStream(&std::cout, MPI_COMM_WORLD, msgfmt));

  EXPECT_DEATH(slic_error(), "")
    << "Rank " << rank << " did not abort on SLIC_ERROR";

  // Finalize slic
  slic::finalize();

  MPI_Finalize();
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  // finalized when exiting main scope
  result = RUN_ALL_TESTS();

  return result;
}

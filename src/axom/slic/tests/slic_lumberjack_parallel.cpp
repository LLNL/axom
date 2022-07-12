// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core/utilities/Utilities.hpp"

#include "axom/slic/interface/slic.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/slic/streams/LumberjackStream.hpp"

#include "gtest/gtest.h"

#include <string>
#include <sstream>

#include "mpi.h"

// namespace alias
namespace slic = axom::slic;

namespace axom
{
namespace slic
{
namespace internal
{
std::ostringstream test_stream;

bool is_stream_empty() { return test_stream.str().empty(); }

void clear()
{
  test_stream.clear();
  test_stream.str("");
  EXPECT_TRUE(is_stream_empty());
}

}  // namespace internal
}  // namespace slic
}  // namespace axom

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
bool has_aborted = false;
void customAbortFunction() { has_aborted = true; }

}  // end anonymous namespace

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
TEST(slic_lumberjack_parallel, test_abort_error_macros)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  slic::enableAbortOnError(); /* enable abort for testing purposes */
  slic::setAbortFunction(customAbortFunction);

  EXPECT_TRUE(slic::internal::is_stream_empty());

  int val = (rank % 2) == 0 ? 42 : -42;

#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)

  SLIC_ERROR("SLIC_ERROR message is logged!");
  EXPECT_TRUE(has_aborted) << "Rank " << rank << " did not abort for SLIC_ERROR";
  has_aborted = false;
  slic::internal::clear();

  SLIC_ERROR_IF(val == 42, "SLIC_ERROR_IF message is logged!");
  EXPECT_TRUE(has_aborted) << "Rank " << rank
                           << " did not abort for SLIC_ERROR_IF";
  has_aborted = false;
  slic::internal::clear();

  axom::slic::setIsRoot(rank == 0);
  SLIC_ERROR_ROOT("SLIC_ERROR_ROOT message is logged!");
  EXPECT_TRUE(has_aborted) << "Rank " << rank
                           << " did not abort for SLIC_ERROR_ROOT";
  has_aborted = false;
  slic::internal::clear();

  SLIC_ERROR_ROOT_IF(val == 42, "SLIC_ERROR_ROOT_IF message is logged!");
  EXPECT_TRUE(has_aborted) << "Rank " << rank
                           << " did not abort for SLIC_ERROR_ROOT_IF";
  has_aborted = false;
  slic::internal::clear();

  SLIC_ASSERT(val < 0);
  EXPECT_TRUE(has_aborted) << "Rank " << rank << " did not abort for SLIC_ASSERT";
  has_aborted = false;
  slic::internal::clear();

  SLIC_ASSERT_MSG(val < 0, "val should be negative!");
  EXPECT_TRUE(has_aborted) << "Rank " << rank
                           << " did not abort for SLIC_ASSERT_MSG";
  has_aborted = false;
  slic::internal::clear();

  axom::slic::setIsRoot(true);
#else
  // SLIC_ASSERT macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  slic::disableAbortOnError(); /* disable abort for testing purposes */
  slic::setAbortFunction(axom::utilities::processAbort);
}

//------------------------------------------------------------------------------
TEST(slic_lumberjack_parallel, test_abort_warning_macros)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  slic::enableAbortOnWarning(); /* enable abort for testing purposes */
  slic::setAbortFunction(customAbortFunction);

  EXPECT_TRUE(slic::internal::is_stream_empty());

  int val = (rank % 2) == 0 ? 42 : -42;

#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)

  SLIC_WARNING("SLIC_WARNING message is logged!");
  EXPECT_TRUE(has_aborted) << "Rank " << rank
                           << " did not abort for SLIC_WARNING";
  has_aborted = false;
  slic::internal::clear();

  SLIC_WARNING_IF(val == 42, "SLIC_WARNING_IF message is logged!");
  EXPECT_TRUE(has_aborted) << "Rank " << rank
                           << " did not abort for SLIC_WARNING_IF";
  has_aborted = false;
  slic::internal::clear();

  axom::slic::setIsRoot(rank == 0);
  SLIC_WARNING_ROOT("SLIC_WARNING_ROOT message is logged!");
  EXPECT_TRUE(has_aborted) << "Rank " << rank
                           << " did not abort for SLIC_WARNING_ROOT";
  has_aborted = false;
  slic::internal::clear();

  SLIC_WARNING_ROOT_IF(val == 42, "SLIC_WARNING_ROOT_IF message is logged!");
  EXPECT_TRUE(has_aborted) << "Rank " << rank
                           << " did not abort for SLIC_WARNING_ROOT_IF";
  has_aborted = false;
  slic::internal::clear();

  SLIC_CHECK(val < 0);
  EXPECT_TRUE(has_aborted) << "Rank " << rank << " did not abort for SLIC_CHECK";
  has_aborted = false;
  slic::internal::clear();

  SLIC_CHECK_MSG(val < 0, "val should be negative!");
  EXPECT_TRUE(has_aborted) << "Rank " << rank
                           << " did not abort for SLIC_CHECK_MSG";
  has_aborted = false;
  slic::internal::clear();

  axom::slic::setIsRoot(true);

#else
  // SLIC_CHECK macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  slic::disableAbortOnWarning(); /* disable abort for testing purposes */
  slic::setAbortFunction(axom::utilities::processAbort);
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  const int RLIMIT = 8;

  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  MPI_Init(&argc, &argv);

  // initialize slic
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Debug);
  slic::disableAbortOnError(); /* disable abort for testing purposes */

  std::string msgfmt = "[<LEVEL>]:;;<MESSAGE>;;\n@@<FILE>\n@@<LINE>";

  slic::addStreamToAllMsgLevels(
    new slic::LumberjackStream(&slic::internal::test_stream,
                               MPI_COMM_WORLD,
                               RLIMIT,
                               msgfmt));

  // finalized when exiting main scope
  result = RUN_ALL_TESTS();

  // Finalize slic
  slic::finalize();

  MPI_Finalize();

  return result;
}

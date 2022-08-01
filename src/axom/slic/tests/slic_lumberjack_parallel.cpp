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

void clear_stream()
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
void check_level(const std::string& msg, const std::string& expected_level)
{
  EXPECT_FALSE(msg.empty());

  // extract level
  size_t start = msg.find("[") + 1;
  size_t end = expected_level.length();
  std::string level = msg.substr(start, end);

  // check with expected
  EXPECT_FALSE(level.empty());
  EXPECT_EQ(level, expected_level);
}

//------------------------------------------------------------------------------
void check_msg(const std::string& msg, const std::string& expected_message)
{
  EXPECT_FALSE(msg.empty());

  // extract message
  size_t start = msg.find(";;") + 2;
  size_t end = expected_message.length();
  std::string mymsg = msg.substr(start, end);

  EXPECT_FALSE(mymsg.empty());
  EXPECT_EQ(mymsg, expected_message);
}

//------------------------------------------------------------------------------
void check_file(const std::string& msg)
{
  EXPECT_FALSE(msg.empty());

  const std::string EXPECTED_FILE = std::string(__FILE__);

  // extract filename
  size_t start = msg.find("@@") + 2;
  size_t end = EXPECTED_FILE.length();
  std::string myfile = msg.substr(start, end);

  EXPECT_FALSE(myfile.empty());
  EXPECT_EQ(myfile, std::string(__FILE__));
}

//------------------------------------------------------------------------------
void check_line(const std::string& msg, int expected_line)
{
  EXPECT_FALSE(msg.empty());

  // extract line
  size_t start = msg.rfind("@@") + 2;
  std::string l = msg.substr(start);
  const int line = std::stoi(l);
  EXPECT_EQ(line, expected_line);
}

//------------------------------------------------------------------------------
bool has_aborted = false;

void custom_abort_function() { has_aborted = true; }
void reset_state()
{
  axom::slic::internal::clear_stream();
  has_aborted = false;
}

}  // end anonymous namespace

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
TEST(slic_lumberjack_parallel, test_abort_error_macros)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nranks;
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  slic::enableAbortOnError(); /* enable abort for testing purposes */
  slic::setAbortFunction(custom_abort_function);

  EXPECT_TRUE(slic::internal::is_stream_empty());

  for(int i = 0; i < nranks; i++)
  {
    int val = rank == i ? 42 : -42;

#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)

    axom::slic::setIsRoot(rank == i);

    if(rank == i)
    {
      SLIC_ERROR("SLIC_ERROR message is logged!");
      EXPECT_TRUE(has_aborted)
        << "Rank " << rank << " did not abort for SLIC_ERROR";
      check_level(slic::internal::test_stream.str(), "ERROR");
      check_msg(slic::internal::test_stream.str(),
                "SLIC_ERROR message is logged!");
      check_file(slic::internal::test_stream.str());
      check_line(slic::internal::test_stream.str(), (__LINE__ - 7));
      reset_state();

      SLIC_ERROR_IF(val == 42, "SLIC_ERROR_IF message is logged!");
      EXPECT_TRUE(has_aborted)
        << "Rank " << rank << " did not abort for SLIC_ERROR_IF";
      EXPECT_FALSE(slic::internal::is_stream_empty());
      check_level(slic::internal::test_stream.str(), "ERROR");
      check_msg(slic::internal::test_stream.str(),
                "SLIC_ERROR_IF message is logged!");
      check_file(slic::internal::test_stream.str());
      check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
      reset_state();

      SLIC_ERROR_ROOT("SLIC_ERROR_ROOT message is logged!");
      EXPECT_TRUE(has_aborted)
        << "Rank " << rank << " did not abort for SLIC_ERROR_ROOT";
      EXPECT_FALSE(slic::internal::is_stream_empty());
      check_level(slic::internal::test_stream.str(), "ERROR");
      check_msg(slic::internal::test_stream.str(),
                "SLIC_ERROR_ROOT message is logged!");
      check_file(slic::internal::test_stream.str());
      check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
      reset_state();

      SLIC_ERROR_ROOT_IF(val == 42, "SLIC_ERROR_ROOT_IF message is logged!");
      EXPECT_TRUE(has_aborted)
        << "Rank " << rank << " did not abort for SLIC_ERROR_ROOT_IF";
      EXPECT_FALSE(slic::internal::is_stream_empty());
      check_level(slic::internal::test_stream.str(), "ERROR");
      check_msg(slic::internal::test_stream.str(),
                "SLIC_ERROR_ROOT_IF message is logged!");
      check_file(slic::internal::test_stream.str());
      check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
      reset_state();

      SLIC_ASSERT(val < 0);
      EXPECT_TRUE(has_aborted)
        << "Rank " << rank << " did not abort for SLIC_ASSERT";
      EXPECT_FALSE(slic::internal::is_stream_empty());
      check_level(slic::internal::test_stream.str(), "ERROR");
      check_msg(slic::internal::test_stream.str(), "Failed Assert: val < 0");
      check_file(slic::internal::test_stream.str());
      check_line(slic::internal::test_stream.str(), (__LINE__ - 7));
      reset_state();

      SLIC_ASSERT_MSG(val < 0, "val should be negative!");
      EXPECT_TRUE(has_aborted)
        << "Rank " << rank << " did not abort for SLIC_ASSERT_MSG";
      EXPECT_FALSE(slic::internal::is_stream_empty());
      check_level(slic::internal::test_stream.str(), "ERROR");
      check_msg(slic::internal::test_stream.str(),
                "Failed Assert: val < 0\nval should be negative!");
      check_file(slic::internal::test_stream.str());
      check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
      reset_state();
    }
    axom::slic::setIsRoot(true);
#else
    // SLIC_ASSERT macros only log messages when AXOM_DEBUG is defined
    AXOM_UNUSED_VAR(val);
    EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  }  // end of for loop

  slic::disableAbortOnError(); /* disable abort for testing purposes */
  slic::setAbortFunction(axom::utilities::processAbort);
}

//------------------------------------------------------------------------------
TEST(slic_lumberjack_parallel, test_abort_warning_macros)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nranks;
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  slic::enableAbortOnWarning(); /* enable abort for testing purposes */
  slic::setAbortFunction(custom_abort_function);

  EXPECT_TRUE(slic::internal::is_stream_empty());

  for(int i = 0; i < nranks; i++)
  {
    int val = rank == i ? 42 : -42;

#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)

    axom::slic::setIsRoot(rank == i);

    if(rank == i)
    {
      SLIC_WARNING("SLIC_WARNING message is logged!");
      EXPECT_TRUE(has_aborted)
        << "Rank " << rank << " did not abort for SLIC_WARNING";
      EXPECT_FALSE(slic::internal::is_stream_empty());
      check_level(slic::internal::test_stream.str(), "WARNING");
      check_msg(slic::internal::test_stream.str(),
                "SLIC_WARNING message is logged!");
      check_file(slic::internal::test_stream.str());
      check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
      reset_state();

      SLIC_WARNING_IF(val == 42, "SLIC_WARNING_IF message is logged!");
      EXPECT_TRUE(has_aborted)
        << "Rank " << rank << " did not abort for SLIC_WARNING_IF";
      EXPECT_FALSE(slic::internal::is_stream_empty());
      check_level(slic::internal::test_stream.str(), "WARNING");
      check_msg(slic::internal::test_stream.str(),
                "SLIC_WARNING_IF message is logged!");
      check_file(slic::internal::test_stream.str());
      check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
      reset_state();

      SLIC_WARNING_ROOT("SLIC_WARNING_ROOT message is logged!");
      EXPECT_TRUE(has_aborted)
        << "Rank " << rank << " did not abort for SLIC_WARNING_ROOT";
      EXPECT_FALSE(slic::internal::is_stream_empty());
      check_level(slic::internal::test_stream.str(), "WARNING");
      check_msg(slic::internal::test_stream.str(),
                "SLIC_WARNING_ROOT message is logged!");
      check_file(slic::internal::test_stream.str());
      check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
      reset_state();

      SLIC_WARNING_ROOT_IF(val == 42, "SLIC_WARNING_ROOT_IF message is logged!");
      EXPECT_TRUE(has_aborted)
        << "Rank " << rank << " did not abort for SLIC_WARNING_ROOT_IF";
      EXPECT_FALSE(slic::internal::is_stream_empty());
      check_level(slic::internal::test_stream.str(), "WARNING");
      check_msg(slic::internal::test_stream.str(),
                "SLIC_WARNING_ROOT_IF message is logged!");
      check_file(slic::internal::test_stream.str());
      check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
      reset_state();

      SLIC_CHECK(val < 0);
      EXPECT_TRUE(has_aborted)
        << "Rank " << rank << " did not abort for SLIC_CHECK";
      EXPECT_FALSE(slic::internal::is_stream_empty());
      check_level(slic::internal::test_stream.str(), "WARNING");
      check_msg(slic::internal::test_stream.str(), "Failed Check: val < 0");
      check_file(slic::internal::test_stream.str());
      check_line(slic::internal::test_stream.str(), (__LINE__ - 7));
      reset_state();

      SLIC_CHECK_MSG(val < 0, "val should be negative!");
      EXPECT_TRUE(has_aborted)
        << "Rank " << rank << " did not abort for SLIC_CHECK_MSG";
      EXPECT_FALSE(slic::internal::is_stream_empty());
      check_level(slic::internal::test_stream.str(), "WARNING");
      check_msg(slic::internal::test_stream.str(),
                "Failed Check: val < 0\nval should be negative!");
      check_file(slic::internal::test_stream.str());
      check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
      reset_state();
    }

    axom::slic::setIsRoot(true);

#else
    // SLIC_CHECK macros only log messages when AXOM_DEBUG is defined
    AXOM_UNUSED_VAR(val);
    EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  }  // end for loop

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

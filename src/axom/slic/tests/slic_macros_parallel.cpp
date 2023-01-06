// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/core/utilities/Utilities.hpp"

#include "axom/slic/interface/slic.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/slic/streams/SynchronizedStream.hpp"
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

void clear_streams()
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
  axom::slic::internal::clear_streams();
  has_aborted = false;
}
}  // end anonymous namespace

//------------------------------------------------------------------------------
// TEST FIXTURE CLASS
//------------------------------------------------------------------------------
class SlicMacrosParallel : public ::testing::TestWithParam<std::string>
{
public:
  void SetUp() override
  {
    stream_type = GetParam();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    const int RLIMIT = 8;

    // initialize slic
    slic::initialize();
    slic::setLoggingMsgLevel(slic::message::Debug);
    slic::disableAbortOnError(); /* disable abort for testing purposes */

    std::string msgfmt = "[<LEVEL>]:;;<MESSAGE>;;\n@@<FILE>\n@@<LINE>";

    if(stream_type == "Lumberjack")
    {
      slic::addStreamToAllMsgLevels(
        new slic::LumberjackStream(&slic::internal::test_stream,
                                   MPI_COMM_WORLD,
                                   RLIMIT,
                                   msgfmt));
    }

    if(stream_type == "Synchronized")
    {
      slic::addStreamToAllMsgLevels(
        new slic::SynchronizedStream(&slic::internal::test_stream,
                                     MPI_COMM_WORLD,
                                     msgfmt));
    }
  }

  void TearDown() override { slic::finalize(); }

  std::string stream_type;

  int rank;
  int nranks;
};

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_error_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());
  SLIC_ERROR("test error message");
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "ERROR");
    check_msg(slic::internal::test_stream.str(), "test error message");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
  }
  slic::internal::clear_streams();

  SLIC_ERROR_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_ERROR_IF(true, "this message is logged!");
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "ERROR");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
  }
  slic::internal::clear_streams();

  // Check selective filtering based on root == false
  axom::slic::setIsRoot(false);
  SLIC_ERROR_ROOT_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Check selective filter based on root == true
  axom::slic::setIsRoot(true);
  SLIC_ERROR_ROOT_IF(true, "this message is logged!");
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "ERROR");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
  }
  slic::internal::clear_streams();

  // is root, but conditional is false -> no message
  axom::slic::setIsRoot(true);
  SLIC_ERROR_ROOT_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // is not root, and conditional is true -> no message
  axom::slic::setIsRoot(false);
  SLIC_ERROR_ROOT_IF(true, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Check for one rank being root
  axom::slic::setIsRoot(rank == 0);
  SLIC_ERROR_ROOT_IF(true, "this message is logged!");
  slic::flushStreams();
  if(rank == 0)
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "ERROR");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
  }
  else
  {
    EXPECT_TRUE(slic::internal::is_stream_empty());
  }
  slic::internal::clear_streams();

  // Check for more than one rank being root for SynchronizedStream
  axom::slic::setIsRoot((rank % 2) == 0);
  SLIC_ERROR_ROOT_IF(true, "this message is logged!");
  slic::flushStreams();
  if(((rank % 2) == 0 && GetParam() == "Synchronized") ||
     (rank == 0 && GetParam() == "Lumberjack"))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "ERROR");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 9));
  }
  else
  {
    EXPECT_TRUE(slic::internal::is_stream_empty());
  }
  slic::internal::clear_streams();
}

//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_warning_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());
  SLIC_WARNING("test warning message");
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "WARNING");
    check_msg(slic::internal::test_stream.str(), "test warning message");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
  }
  slic::internal::clear_streams();

  SLIC_WARNING_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_WARNING_IF(true, "this message is logged!");
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "WARNING");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
  }
  slic::internal::clear_streams();

  // Check selective filtering based on root == false
  axom::slic::setIsRoot(false);
  SLIC_WARNING_ROOT_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Check selective filter based on root == true
  axom::slic::setIsRoot(true);
  SLIC_WARNING_ROOT_IF(true, "this message is logged!");
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "WARNING");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
  }
  slic::internal::clear_streams();

  // is root, but conditional is false -> no message
  axom::slic::setIsRoot(true);
  SLIC_WARNING_ROOT_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // is not root, and conditional is true -> no message
  axom::slic::setIsRoot(false);
  SLIC_WARNING_ROOT_IF(true, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Check for one rank being root
  axom::slic::setIsRoot(rank == 0);
  SLIC_WARNING_ROOT_IF(true, "this message is logged!");
  slic::flushStreams();
  if(rank == 0)
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "WARNING");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
  }
  else
  {
    EXPECT_TRUE(slic::internal::is_stream_empty());
  }
  slic::internal::clear_streams();

  // Check for more than one rank being root for SynchronizedStream
  axom::slic::setIsRoot((rank % 2) == 0);
  SLIC_WARNING_ROOT_IF(true, "this message is logged!");
  slic::flushStreams();
  if(((rank % 2) == 0 && GetParam() == "Synchronized") ||
     (rank == 0 && GetParam() == "Lumberjack"))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "WARNING");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 9));
  }
  else
  {
    EXPECT_TRUE(slic::internal::is_stream_empty());
  }
  slic::internal::clear_streams();
}

//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_info_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());
  SLIC_INFO("test info message");
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "INFO");
    check_msg(slic::internal::test_stream.str(), "test info message");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), __LINE__ - 8);
  }
  slic::internal::clear_streams();

  SLIC_INFO_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_INFO_IF(true, "this message is logged!");
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "INFO");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
  }
  slic::internal::clear_streams();

  // Check selective filtering based on root == false
  axom::slic::setIsRoot(false);
  SLIC_INFO_ROOT_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Check selective filter based on root == true
  axom::slic::setIsRoot(true);
  SLIC_INFO_ROOT_IF(true, "this message is logged!");
  slic::flushStreams();
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "INFO");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
  }
  slic::internal::clear_streams();

  // is root, but conditional is false -> no message
  axom::slic::setIsRoot(true);
  SLIC_INFO_ROOT_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // is not root, and conditional is true -> no message
  axom::slic::setIsRoot(false);
  SLIC_INFO_ROOT_IF(true, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Check for one rank being root
  axom::slic::setIsRoot(rank == 0);
  SLIC_INFO_ROOT_IF(true, "this message is logged!");
  slic::flushStreams();
  if(rank == 0)
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "INFO");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
  }
  else
  {
    EXPECT_TRUE(slic::internal::is_stream_empty());
  }
  slic::internal::clear_streams();

  // Check for more than one rank being root
  axom::slic::setIsRoot((rank % 2) == 0);
  SLIC_INFO_ROOT_IF(true, "this message is logged!");
  slic::flushStreams();
  if(((rank % 2) == 0 && GetParam() == "Synchronized") ||
     (rank == 0 && GetParam() == "Lumberjack"))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "INFO");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 9));
  }
  else
  {
    EXPECT_TRUE(slic::internal::is_stream_empty());
  }
  slic::internal::clear_streams();
}

//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_debug_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());
  SLIC_DEBUG("test debug message");
  slic::flushStreams();
#ifdef AXOM_DEBUG
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "DEBUG");
    check_msg(slic::internal::test_stream.str(), "test debug message");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 9));
  }
  slic::internal::clear_streams();
#else
  // SLIC_DEBUG macros only log messages when AXOM_DEBUG is defined
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  SLIC_DEBUG_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_DEBUG_IF(true, "this message is logged!");
  slic::flushStreams();
#ifdef AXOM_DEBUG
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "DEBUG");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 9));
  }
  slic::internal::clear_streams();
#else
  // SLIC_DEBUG macros only log messages when AXOM_DEBUG is defined
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  // Check selective filtering based on root == false
  axom::slic::setIsRoot(false);
  SLIC_DEBUG_ROOT_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Check selective filter based on root == true
  axom::slic::setIsRoot(true);
  SLIC_DEBUG_ROOT_IF(true, "this message is logged!");
  slic::flushStreams();
#ifdef AXOM_DEBUG
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "DEBUG");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 9));
  }
  slic::internal::clear_streams();
#else
  // SLIC_DEBUG macros only log messages when AXOM_DEBUG is defined
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  // is root, but conditional is false -> no message
  axom::slic::setIsRoot(true);
  SLIC_DEBUG_ROOT_IF(false, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // is not root, and conditional is true -> no message
  axom::slic::setIsRoot(false);
  SLIC_DEBUG_ROOT_IF(true, "this message should not be logged!");
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Check for one rank being root
  axom::slic::setIsRoot(rank == 0);
  SLIC_DEBUG_ROOT_IF(true, "this message is logged!");
  slic::flushStreams();
#ifdef AXOM_DEBUG
  if(rank == 0)
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "DEBUG");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 9));
  }
  else
  {
    EXPECT_TRUE(slic::internal::is_stream_empty());
  }
  slic::internal::clear_streams();
#else
  // SLIC_DEBUG macros only log messages when AXOM_DEBUG is defined
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  // Check for more than one rank being root
  axom::slic::setIsRoot((rank % 2) == 0);
  SLIC_DEBUG_ROOT_IF(true, "this message is logged!");
  slic::flushStreams();
#ifdef AXOM_DEBUG
  if(((rank % 2) == 0 && GetParam() == "Synchronized") ||
     (rank == 0 && GetParam() == "Lumberjack"))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "DEBUG");
    check_msg(slic::internal::test_stream.str(), "this message is logged!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 10));
  }
  else
  {
    EXPECT_TRUE(slic::internal::is_stream_empty());
  }
  slic::internal::clear_streams();
#else
  // SLIC_DEBUG macros only log messages when AXOM_DEBUG is defined
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif
}

TEST_P(SlicMacrosParallel, test_abort_error_macros)
{
  const int NUM_ABORT_STATES = 2;

  slic::enableAbortOnError(); /* enable abort for testing purposes */
  slic::setAbortFunction(custom_abort_function);

  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Test for each rank
  for(int i = 0; i < nranks; i++)
  {
    // Test abort enabled/disabled
    for(int j = 0; j < NUM_ABORT_STATES; j++)
    {
      if(j == 0)
      {
        slic::enableAbortOnError();
      }
      else
      {
        slic::disableAbortOnError();
      }

      bool abort_enabled = slic::isAbortOnErrorsEnabled();

#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)

      axom::slic::setIsRoot(rank == i);

      if(rank == i)
      {
        SLIC_ERROR("SLIC_ERROR message is logged!");
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        check_level(slic::internal::test_stream.str(), "ERROR");
        check_msg(slic::internal::test_stream.str(),
                  "SLIC_ERROR message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), (__LINE__ - 7));
        reset_state();

        int val = rank == i ? 42 : -42;
        SLIC_ERROR_IF(val == 42, "SLIC_ERROR_IF message is logged!");
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::is_stream_empty());
        check_level(slic::internal::test_stream.str(), "ERROR");
        check_msg(slic::internal::test_stream.str(),
                  "SLIC_ERROR_IF message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
        reset_state();

        SLIC_ERROR_ROOT("SLIC_ERROR_ROOT message is logged!");
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::is_stream_empty());
        check_level(slic::internal::test_stream.str(), "ERROR");
        check_msg(slic::internal::test_stream.str(),
                  "SLIC_ERROR_ROOT message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
        reset_state();

        SLIC_ERROR_ROOT_IF(val == 42, "SLIC_ERROR_ROOT_IF message is logged!");
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::is_stream_empty());
        check_level(slic::internal::test_stream.str(), "ERROR");
        check_msg(slic::internal::test_stream.str(),
                  "SLIC_ERROR_ROOT_IF message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
        reset_state();

        SLIC_ASSERT(val < 0);
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::is_stream_empty());
        check_level(slic::internal::test_stream.str(), "ERROR");
        check_msg(slic::internal::test_stream.str(), "Failed Assert: val < 0");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), (__LINE__ - 7));
        reset_state();

        SLIC_ASSERT_MSG(val < 0, "val should be negative!");
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
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

      // Quiet warning about has_aborted and reset_state never being referenced
      EXPECT_EQ(has_aborted, abort_enabled);
      reset_state();

      EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

    }  //end NUM_ABORT_STATES loop
  }    // end nranks loop

  slic::disableAbortOnError(); /* disable abort for testing purposes */
  slic::setAbortFunction(axom::utilities::processAbort);
}

//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_abort_warning_macros)
{
  const int NUM_ABORT_STATES = 2;

  slic::enableAbortOnWarning(); /* enable abort for testing purposes */
  slic::setAbortFunction(custom_abort_function);

  EXPECT_TRUE(slic::internal::is_stream_empty());

  // Test for each rank
  for(int i = 0; i < nranks; i++)
  {
    // Test abort enabled/disabled
    for(int j = 0; j < NUM_ABORT_STATES; j++)
    {
      if(j == 0)
      {
        slic::enableAbortOnWarning();
      }
      else
      {
        slic::disableAbortOnWarning();
      }

      bool abort_enabled = slic::isAbortOnWarningsEnabled();

#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)

      axom::slic::setIsRoot(rank == i);

      if(rank == i)
      {
        SLIC_WARNING("SLIC_WARNING message is logged!");
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::is_stream_empty());
        check_level(slic::internal::test_stream.str(), "WARNING");
        check_msg(slic::internal::test_stream.str(),
                  "SLIC_WARNING message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
        reset_state();

        int val = rank == i ? 42 : -42;
        SLIC_WARNING_IF(val == 42, "SLIC_WARNING_IF message is logged!");
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::is_stream_empty());
        check_level(slic::internal::test_stream.str(), "WARNING");
        check_msg(slic::internal::test_stream.str(),
                  "SLIC_WARNING_IF message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
        reset_state();

        SLIC_WARNING_ROOT("SLIC_WARNING_ROOT message is logged!");
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::is_stream_empty());
        check_level(slic::internal::test_stream.str(), "WARNING");
        check_msg(slic::internal::test_stream.str(),
                  "SLIC_WARNING_ROOT message is logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
        reset_state();

        SLIC_WARNING_ROOT_IF(val == 42, "SLIC_WARNING_ROOT_IF msg logged!");
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::is_stream_empty());
        check_level(slic::internal::test_stream.str(), "WARNING");
        check_msg(slic::internal::test_stream.str(),
                  "SLIC_WARNING_ROOT_IF msg logged!");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), (__LINE__ - 8));
        reset_state();

        SLIC_CHECK(val < 0);
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
        EXPECT_FALSE(slic::internal::is_stream_empty());
        check_level(slic::internal::test_stream.str(), "WARNING");
        check_msg(slic::internal::test_stream.str(), "Failed Check: val < 0");
        check_file(slic::internal::test_stream.str());
        check_line(slic::internal::test_stream.str(), (__LINE__ - 7));
        reset_state();

        SLIC_CHECK_MSG(val < 0, "val should be negative!");
        slic::outputLocalMessages();
        EXPECT_EQ(has_aborted, abort_enabled);
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
      EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

    }  //end NUM_ABORT_STATES loop
  }    // end nranks loop

  slic::disableAbortOnWarning(); /* disable abort for testing purposes */
  slic::setAbortFunction(axom::utilities::processAbort);
}

//------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_assert_macros)
{
  slic::internal::clear_streams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  constexpr int val = 42;
  SLIC_ASSERT(val < 0);
  slic::flushStreams();
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "ERROR");
    check_msg(slic::internal::test_stream.str(), "Failed Assert: val < 0");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 9));
  }
  slic::internal::clear_streams();
#else
  // SLIC_ASSERT macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  SLIC_ASSERT(val > 0);
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_ASSERT_MSG(val < 0, "val should be negative!");
  slic::flushStreams();
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "ERROR");
    check_msg(slic::internal::test_stream.str(),
              "Failed Assert: val < 0\nval should be negative!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 10));
  }
  slic::internal::clear_streams();
#else
  // SLIC_ASSERT macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif
}

// ------------------------------------------------------------------------------
TEST_P(SlicMacrosParallel, test_check_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());

  constexpr int val = 42;
  SLIC_CHECK(val < 0);
  slic::flushStreams();
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "WARNING");
    check_msg(slic::internal::test_stream.str(), "Failed Check: val < 0");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 9));
  }
  slic::internal::clear_streams();
#else
  // SLIC_CHECK macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  SLIC_CHECK(val > 0);
  slic::flushStreams();
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_CHECK_MSG(val < 0, "val should be negative!");
  slic::flushStreams();
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  if(GetParam() == "Synchronized" || (GetParam() == "Lumberjack" && rank == 0))
  {
    EXPECT_FALSE(slic::internal::is_stream_empty());
    check_level(slic::internal::test_stream.str(), "WARNING");
    check_msg(slic::internal::test_stream.str(),
              "Failed Check: val < 0\nval should be negative!");
    check_file(slic::internal::test_stream.str());
    check_line(slic::internal::test_stream.str(), (__LINE__ - 10));
  }
  slic::internal::clear_streams();
#else
  // SLIC_CHECK macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif
}

//------------------------------------------------------------------------------
const std::string parallel_streams[] = {"Synchronized", "Lumberjack"};
INSTANTIATE_TEST_SUITE_P(core_memory_management,
                         SlicMacrosParallel,
                         ::testing::ValuesIn(parallel_streams));

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  MPI_Init(&argc, &argv);

  // finalized when exiting main scope
  result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}

// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// axom includes
#include "axom/config.hpp"

// slic includes
#include "axom/slic/interface/slic.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"

// gtest includes
#include "gtest/gtest.h"  // for gtest macros

// C/C++ includes
#include <string>   // for C++ string
#include <sstream>  // for std::ostringstream

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

}  // end anonymous namespace

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(slic_macros, test_error_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());
  SLIC_ERROR("test error message");
  EXPECT_FALSE(slic::internal::is_stream_empty());
  check_level(slic::internal::test_stream.str(), "ERROR");
  check_msg(slic::internal::test_stream.str(), "test error message");
  check_file(slic::internal::test_stream.str());
  check_line(slic::internal::test_stream.str(), (__LINE__ - 5));
  slic::internal::clear();

  SLIC_ERROR_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_ERROR_IF(true, "this message is logged!");
  EXPECT_FALSE(slic::internal::is_stream_empty());
  check_level(slic::internal::test_stream.str(), "ERROR");
  check_msg(slic::internal::test_stream.str(), "this message is logged!");
  check_file(slic::internal::test_stream.str());
  check_line(slic::internal::test_stream.str(), (__LINE__ - 5));
  slic::internal::clear();
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_warning_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());
  SLIC_WARNING("test warning message");
  EXPECT_FALSE(slic::internal::is_stream_empty());
  check_level(slic::internal::test_stream.str(), "WARNING");
  check_msg(slic::internal::test_stream.str(), "test warning message");
  check_file(slic::internal::test_stream.str());
  check_line(slic::internal::test_stream.str(), (__LINE__ - 5));
  slic::internal::clear();

  SLIC_WARNING_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_WARNING_IF(true, "this message is logged!");
  EXPECT_FALSE(slic::internal::is_stream_empty());
  check_level(slic::internal::test_stream.str(), "WARNING");
  check_msg(slic::internal::test_stream.str(), "this message is logged!");
  check_file(slic::internal::test_stream.str());
  check_line(slic::internal::test_stream.str(), (__LINE__ - 5));
  slic::internal::clear();
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_info_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());
  SLIC_INFO("test info message");
  EXPECT_FALSE(slic::internal::is_stream_empty());
  check_level(slic::internal::test_stream.str(), "INFO");
  check_msg(slic::internal::test_stream.str(), "test info message");
  check_file(slic::internal::test_stream.str());
  check_line(slic::internal::test_stream.str(), __LINE__ - 5);
  slic::internal::clear();

  SLIC_INFO_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_INFO_IF(true, "this message is logged!");
  EXPECT_FALSE(slic::internal::is_stream_empty());
  check_level(slic::internal::test_stream.str(), "INFO");
  check_msg(slic::internal::test_stream.str(), "this message is logged!");
  check_file(slic::internal::test_stream.str());
  check_line(slic::internal::test_stream.str(), (__LINE__ - 5));
  slic::internal::clear();
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_debug_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());
  SLIC_DEBUG("test debug message");
#ifdef AXOM_DEBUG
  EXPECT_FALSE(slic::internal::is_stream_empty());
  check_level(slic::internal::test_stream.str(), "DEBUG");
  check_msg(slic::internal::test_stream.str(), "test debug message");
  check_file(slic::internal::test_stream.str());
  check_line(slic::internal::test_stream.str(), (__LINE__ - 6));
  slic::internal::clear();
#else
  // SLIC_DEBUG macros only log messages when AXOM_DEBUG is defined
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  SLIC_DEBUG_IF(false, "this message should not be logged!");
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_DEBUG_IF(true, "this message is logged!");
#ifdef AXOM_DEBUG
  EXPECT_FALSE(slic::internal::is_stream_empty());
  check_level(slic::internal::test_stream.str(), "DEBUG");
  check_msg(slic::internal::test_stream.str(), "this message is logged!");
  check_file(slic::internal::test_stream.str());
  check_line(slic::internal::test_stream.str(), (__LINE__ - 6));
  slic::internal::clear();
#else
  // SLIC_DEBUG macros only log messages when AXOM_DEBUG is defined
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_assert_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());

  constexpr int val = 42;
  SLIC_ASSERT(val < 0);
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  EXPECT_FALSE(slic::internal::is_stream_empty());
  check_level(slic::internal::test_stream.str(), "ERROR");
  check_msg(slic::internal::test_stream.str(), "Failed Assert: val < 0");
  check_file(slic::internal::test_stream.str());
  check_line(slic::internal::test_stream.str(), (__LINE__ - 6));
  slic::internal::clear();
#else
  // SLIC_ASSERT macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  SLIC_ASSERT(val > 0);
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_ASSERT_MSG(val < 0, "val should be negative!");
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  EXPECT_FALSE(slic::internal::is_stream_empty());
  check_level(slic::internal::test_stream.str(), "ERROR");
  check_msg(slic::internal::test_stream.str(),
            "Failed Assert: val < 0\nval should be negative!");
  check_file(slic::internal::test_stream.str());
  check_line(slic::internal::test_stream.str(), (__LINE__ - 7));
  slic::internal::clear();
#else
  // SLIC_ASSERT macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif
}

//------------------------------------------------------------------------------
TEST(slic_macros, test_check_macros)
{
  EXPECT_TRUE(slic::internal::is_stream_empty());

  constexpr int val = 42;
  SLIC_CHECK(val < 0);
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  EXPECT_FALSE(slic::internal::is_stream_empty());
  check_level(slic::internal::test_stream.str(), "WARNING");
  check_msg(slic::internal::test_stream.str(), "Failed Check: val < 0");
  check_file(slic::internal::test_stream.str());
  check_line(slic::internal::test_stream.str(), (__LINE__ - 6));
  slic::internal::clear();
#else
  // SLIC_CHECK macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif

  SLIC_CHECK(val > 0);
  EXPECT_TRUE(slic::internal::is_stream_empty());

  SLIC_CHECK_MSG(val < 0, "val should be negative!");
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)
  EXPECT_FALSE(slic::internal::is_stream_empty());
  check_level(slic::internal::test_stream.str(), "WARNING");
  check_msg(slic::internal::test_stream.str(),
            "Failed Check: val < 0\nval should be negative!");
  check_file(slic::internal::test_stream.str());
  check_line(slic::internal::test_stream.str(), (__LINE__ - 7));
  slic::internal::clear();
#else
  // SLIC_CHECK macros only log messages when AXOM_DEBUG is defined
  AXOM_UNUSED_VAR(val);
  EXPECT_TRUE(slic::internal::is_stream_empty());
#endif
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  // initialize slic
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Debug);
  slic::disableAbortOnError(); /* disable abort for testing purposes */

  std::string msgfmt = "[<LEVEL>]:;;<MESSAGE>;;\n@@<FILE>\n@@<LINE>";

  slic::addStreamToAllMsgLevels(
    new slic::GenericOutputStream(&slic::internal::test_stream, msgfmt));

  // finalized when exiting main scope
  result = RUN_ALL_TESTS();

  // Finalize slic
  slic::finalize();

  return result;
}

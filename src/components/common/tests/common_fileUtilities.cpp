/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#include "gtest/gtest.h"
#include <fstream>

#include "common/FileUtilities.hpp"

TEST(gtest_common_fileUtilities,getCWD_smoke)
{
  // This test just checks that we can call the getCWD function
  // It does not in any way confirm the results

  std::cout<<"Checking that we can call getCWD()" << std::endl;

  std::string cwd = asctoolkit::utilities::filesystem::getCWD();

  std::cout <<" CWD is: " << cwd << std::endl;

  EXPECT_TRUE(true);
}



TEST(gtest_common_fileUtilities,common_fileUtil_joinPath)
{
  using namespace asctoolkit::utilities::filesystem;

  std::string fdir = "abc";
  std::string fdirWithSlash = "abc/";
  std::string fname = "def";

  std::string fullfile = "abc/def";

  std::cout<< "Testing joinPath file utility" << std::endl;

  EXPECT_EQ( fullfile, joinPath( fdir,fname) );

  EXPECT_EQ( fullfile, joinPath( fdirWithSlash,fname) );


  std::string fnameWithSubdir = "def/ghi";
  EXPECT_EQ( "abc/def/ghi", joinPath( fdir,fnameWithSubdir) );

}


TEST(gtest_common_fileUtilities,common_fileUtil_pathExists)
{
  using namespace asctoolkit::utilities::filesystem;

  std::cout<<"Testing pathExists on file that we know is present (the cwd)."<< std::endl;
  const std::string missingFile = "m_i_s_s_i_n_g__f_i_l_e";

  std::string cwd = asctoolkit::utilities::filesystem::getCWD();
  EXPECT_TRUE( pathExists(cwd) );
  EXPECT_FALSE( pathExists( joinPath(cwd,missingFile) ) );
}


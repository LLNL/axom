/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "gtest/gtest.h"
#include <fstream>

#include "axom/core/utilities/FileUtilities.hpp"

TEST(core_fileUtilities,getCWD_smoke)
{
  // This test just checks that we can call the getCWD function
  // It does not in any way confirm the results

  std::cout<<"Checking that we can call getCWD()" << std::endl;

  std::string cwd = axom::utilities::filesystem::getCWD();

  std::cout <<" CWD is: " << cwd << std::endl;

  EXPECT_TRUE(true);
}



TEST(core_fileUtilities,joinPath)
{
  using namespace axom::utilities::filesystem;

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


TEST(core_fileUtilities,pathExists)
{
  using namespace axom::utilities::filesystem;

  std::cout<<"Testing pathExists on file "
           << "that we know is present (the cwd)."<< std::endl;

  const std::string missingFile = "m_i_s_s_i_n_g__f_i_l_e";

  std::string cwd = axom::utilities::filesystem::getCWD();
  EXPECT_TRUE( pathExists(cwd) );
  EXPECT_FALSE( pathExists( joinPath(cwd,missingFile) ) );
}

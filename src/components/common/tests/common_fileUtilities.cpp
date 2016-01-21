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

namespace {

    // Use a file that we know is available in the git repo (as of Jan 2016).
    // It is a data file from slam's Unstructured mesh example
    // TODO: This path should be modified once
    //       we move the data files into their own repository

    const std::string presentFile = "ball_1.vtk";
    const std::string presentDir = "src/components/slam/data";
    const std::string presentDirWithSlash = presentDir + "/";

    const std::string missingFile = "m_i_s_s_i_n_g__f_i_l_e";
}


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

    std::cout<< "Testing pathExists file utility on file that we know is present." << std::endl;

    // not sure which directory code will be called from
    // check src and two parents
    std::string fullFile = presentDirWithSlash + presentFile;

    bool valid = pathExists( fullFile );
    for(int i=0; i < 3; ++i)
    {
        fullFile = "../" + fullFile;
        valid = valid || pathExists( fullFile );
    }

    EXPECT_TRUE( valid);

    EXPECT_FALSE( pathExists(missingFile) );
}


//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}

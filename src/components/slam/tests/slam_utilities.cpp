/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/**
 * \file testStaticVariableRelation.cxx
 *
 *  Created on: Apr 29, 2015
 *      Author: weiss27
 */


#include <iostream>
#include <iterator>
#include <fstream>

#include "gtest/gtest.h"


#include "common/FileUtilities.hpp"
#include "slam/Utilities.hpp"


namespace {

    // Use a file that we know is available in the git repo (as of Jan 2016).
    // It is a data file from slam's Unstructured mesh example
    // TODO: This path should be modified once
    //       we move the data files into their own repository
    // See also: common_fileUtilities.cpp test in common directory

    const std::string presentFile = "ball_1.vtk";
    const std::string presentDir = "src/components/slam/data";
    const std::string presentDirWithSlash = presentDir + "/";

    const std::string missingFile = "m_i_s_s_i_n_g__f_i_l_e";
}

TEST(gtest_slam_utilities,findingAncestorPaths)
{
  std::cout << "\n****** Testing funtion that recursively finds a valid path in cwd ancestors." << std::endl;

    using namespace asctoolkit::slam::util;
    using namespace asctoolkit::utilities::filesystem;

    std::string path1 = joinPath(presentDir, presentFile);
    std::string ps1 = findFileInAncestorDirs(path1);

    std::ifstream fileStreamP1( ps1.c_str() );
    EXPECT_TRUE( fileStreamP1);
    fileStreamP1.close();
    std::cout <<"For file: " << path1 << " valid path was: " << ps1 << std::endl;


    /// Second test -- check if code works when 'slash' is present in path
    std::string path2 = joinPath(presentDirWithSlash, presentFile);
    std::string ps2 = findFileInAncestorDirs(path2);

    std::ifstream fileStreamP2( ps2.c_str() );
    EXPECT_TRUE( fileStreamP2);
    fileStreamP2.close();
    std::cout <<"For file: " << path2 << " valid path was: " << ps2 << std::endl;

    /// Second test -- check if code works when 'slash' is present in path
    std::string path3 = joinPath("../" + presentDir, presentFile);
    std::string ps3 = findFileInAncestorDirs(path3);

    std::ifstream fileStreamP3( ps3.c_str() );
    EXPECT_TRUE( fileStreamP3);
    fileStreamP3.close();
    std::cout <<"For file: " << path3 << " valid path was: " << ps3 << std::endl;



    /// Third test -- check for file that should not be there.

    std::string path4 = joinPath(presentDir, missingFile);

    std::cout<<"\n** We are about to look for a nonexistent file."
             <<"Please ignore the log warning. ***"
             << std::endl;

    std::string ps4 = findFileInAncestorDirs(path4);
    EXPECT_EQ( path4, ps4);

    std::ifstream fileStreamP4( ps4.c_str() );
    EXPECT_FALSE( fileStreamP4);
    fileStreamP4.close();
    std::cout <<"There was no valid path for file: " << path4 << ".\n"
              << "Function returned: " << ps4 << std::endl;

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

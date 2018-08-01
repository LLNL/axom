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


/**
 * \file slam_utilities
 *
 * \brief Unit tests for utility functions in Slam
 */


#include <iostream>
#include <iterator>
#include <fstream>

#include "gtest/gtest.h"


#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/slam/Utilities.hpp"
#include "axom/slic/interface/slic.hpp"

#ifdef WIN32
    #include <direct.h>
    #define ChangeCurrentDir _chdir
#else
    #include <unistd.h>
    #define ChangeCurrentDir  chdir
#endif


namespace
{

// Use a file that we know is available in the git repo (as of Jan 2016).
// It is a data file from slam's Unstructured mesh example
// TODO: This path should be modified once
//       we move the data files into their own repository

const std::string presentFile = "ball_1.vtk";
const std::string presentDir = "src/axom/slam/data";
const std::string presentDirWithSlash = presentDir + "/";

const std::string missingFile = "m_i_s_s_i_n_g__f_i_l_e";
}

TEST(slam_utilities,findingAncestorPaths)
{
  SLIC_INFO("Testing function that recursively finds"
            << " a valid path in cwd ancestors.");

  using namespace axom::slam::util;
  using namespace axom::utilities::filesystem;

  std::string path1 = joinPath(presentDir, presentFile);
  std::string ps1 = findFileInAncestorDirs(path1);

  std::ifstream fileStreamP1( ps1.c_str() );
  //
  EXPECT_TRUE( fileStreamP1.is_open());
  fileStreamP1.close();
  SLIC_INFO("Valid path for file: " << path1 << " was: " << ps1);


  /// Second test -- check if code works when 'slash' is present in path
  std::string path2 = joinPath(presentDirWithSlash, presentFile);
  std::string ps2 = findFileInAncestorDirs(path2);

  std::ifstream fileStreamP2( ps2.c_str() );
  EXPECT_TRUE( fileStreamP2.is_open());
  fileStreamP2.close();
  SLIC_INFO("Valid path for file: " << path2 << " was: " << ps2 );

  /// Second test -- check if code works when 'slash' is present in path
  std::string path3 = joinPath("../" + presentDir, presentFile);
  std::string ps3 = findFileInAncestorDirs(path3);

  std::ifstream fileStreamP3( ps3.c_str() );
  EXPECT_TRUE( fileStreamP3.is_open());
  fileStreamP3.close();
  SLIC_INFO("Valid path for file: " << path3 << " was: " << ps3 );



  /// Third test -- check for file that should not be there.

  std::string path4 = joinPath(presentDir, missingFile);

  SLIC_INFO("** About to look for a nonexistent file. "
            << " Ignore the log warning. ***");
  std::string ps4 = findFileInAncestorDirs(path4);
  EXPECT_EQ( path4, ps4);

  std::ifstream fileStreamP4( ps4.c_str() );
  EXPECT_FALSE( fileStreamP4.is_open());
  fileStreamP4.close();
  SLIC_INFO("No valid path for file '"
            << path4 << "'. Function returned " << ps4 );

}

//----------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // Change the directory to one that we know will contain the desired file
  SLIC_ERROR_IF(argc !=2,
                "slam_utilities requires a parameter for the working directory");
  int err = ChangeCurrentDir( argv[1] );
  SLIC_ERROR_IF( err != 0, "chdir failed");

  result = RUN_ALL_TESTS();

  return result;
}

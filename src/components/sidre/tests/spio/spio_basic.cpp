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

#include "axom/config.hpp"   // for AXOM_USE_HDF5
#include "sidre/sidre.hpp"
#include "sidre/IOManager.hpp"

using axom::sidre::Group;
using axom::sidre::DataStore;
using axom::sidre::IOManager;

using axom::sidre::detail::sidre_int64;


//------------------------------------------------------------------------------

// Check that the mapping between sidre and relay protocols is accurate
TEST(spio_basic, root_name)
{
  // Set up map from sidre protocols to expected relay protocols
  std::map<std::string, std::string> protocolMap;
  protocolMap["sidre_hdf5"]         = "hdf5";
  protocolMap["conduit_hdf5"]       = "hdf5";
  protocolMap["sidre_json"]         = "json";
  protocolMap["conduit_bin"]        = "json";
  protocolMap["json"]               = "json";
  protocolMap["sidre_conduit_json"] = "conduit_json";

  typedef std::map<std::string, std::string>::const_iterator MapIt;
  for(MapIt it = protocolMap.begin() ; it != protocolMap.end() ; ++it)
  {
    const std::string& sidreProtocol = it->first;
    const std::string& expRelayProtocol = it->second;

    std::string relayProtocol =
      IOManager::correspondingRelayProtocol(sidreProtocol);

    EXPECT_EQ(expRelayProtocol, relayProtocol);
  }
}

//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}

// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/core.hpp"
#include "axom/curves.hpp"
#include "axom/slic.hpp"

// namespace aliases
namespace slic = axom::slic;

int main ( int AXOM_NOT_USED(argc), char** AXOM_NOT_USED(argv) )
{
  // STEP 0: initialize the slic logging environment
  slic::initialize();
  slic::setLoggingMsgLevel( slic::message::Debug );
  slic::addStreamToAllMsgLevels(
      new slic::GenericOutputStream( &std::cout, "[<LEVEL>]: <MESSAGE>\n" ) );

  SLIC_INFO( axom::curves::dummy() );

  slic::finalize();
  return 0;
}




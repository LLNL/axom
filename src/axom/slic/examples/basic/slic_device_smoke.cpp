// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/CLI11.hpp"

#include <vector>
#include <map>

// RAJA
#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"
#endif

using namespace axom;

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);

  //----------------------------------------------------------------------------
  // STEP 0: Initialize logger
  //----------------------------------------------------------------------------
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Debug);
  slic::disableAbortOnError();

  //----------------------------------------------------------------------------
  // STEP 1: Create log streams
  //----------------------------------------------------------------------------

  std::string msgfmt = "[<LEVEL>]:;;<MESSAGE>;;\n@@<FILE>\n@@<LINE>";
  slic::addStreamToAllMsgLevels(new slic::GenericOutputStream(&std::cout, msgfmt));

  //----------------------------------------------------------------------------
  // STEP 2: Call SLIC macro on device
  //----------------------------------------------------------------------------

  constexpr int BLOCK_SIZE = 256;

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(AXOM_USE_CUDA)
  axom::for_all<axom::CUDA_EXEC<BLOCK_SIZE>>(
    1,
    AXOM_LAMBDA(axom::IndexType i) { printf("Test print on CUDA device"); });
  #endif

  #if defined(AXOM_USE_HIP)
  axom::for_all<axom::HIP_EXEC<BLOCK_SIZE>>(
    1,
    AXOM_LAMBDA(axom::IndexType i) { printf("Test print on HIP device"); });
  #endif
#endif

  // Finalize slic
  slic::finalize();

  return 0;
}
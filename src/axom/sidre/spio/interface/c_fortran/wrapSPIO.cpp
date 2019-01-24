// wrapSPIO.cpp
// This is generated code, do not edit
//
// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
//
// Produced at the Lawrence Livermore National Laboratory
//
// LLNL-CODE-741217
//
// All rights reserved.
//
// This file is part of Axom.
//
// For details about use and distribution, please read axom/LICENSE.
//
#include <stdlib.h>
#include "axom/sidre/spio/IOManager.hpp"
#include "typesSPIO.h"

// splicer begin CXX_definitions
// splicer end CXX_definitions

extern "C" {

// splicer begin C_definitions
// splicer end C_definitions

// Release C++ allocated memory.
void SPIO_SHROUD_memory_destructor(SPI_SHROUD_capsule_data* cap)
{
  void* ptr = cap->addr;
  switch (cap->idtor)
  {
  case 0:     // --none--
  {
    // Nothing to delete
    break;
  }
  case 1:     // axom::sidre::IOManager
  {
    axom::sidre::IOManager* cxx_ptr =
      reinterpret_cast<axom::sidre::IOManager*>(ptr);
    delete cxx_ptr;
    break;
  }
  default:
  {
    // Unexpected case in destructor
    break;
  }
  }
  cap->addr = NULL;
  cap->idtor = 0;    // avoid deleting again
}

}  // extern "C"

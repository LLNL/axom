/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/**
 *  \file SidreWrapperHelpers.hpp
 *
 *  \brief File used to contain helper functions for Fortran/C API wrappers.
 *         User code should not include this file.
 *
 */

// Standard C++ headers
#include <string>

// Other toolkit project headers
#include "slic/slic.hpp"

// SiDRe project headers
#include "SidreWrapperHelpers.hpp"


namespace asctoolkit
{
namespace sidre
{

static int global_int;
static void *global_void;

/*!
 * \brief Return DataView for a Fortran allocatable.
 *
 * The Fortran allocatable array is the buffer for the DataView.
 */
void *register_allocatable(DataGroup *group,
			   const std::string &name,
			   void *array, int atk_type, int rank)
{
  SLIC_ASSERT( !name.empty() );
  SLIC_ASSERT_MSG( group->hasView(name) == false, "name == " << name );

#if 0
  if ( name.empty() || group->hasView(name) ) 
  {
    return ATK_NULLPTR;
  }
  else 
  {
    DataBuffer * buff = group->getDataStore()->createBuffer();
    DataView * const view = new DataView( name, this, buff);
    buff->attachView(view);
    return attachView(view);
  }
#endif
    global_void = array;
    global_int = atk_type;
    global_int = rank;
    return NULL;
}

} /* end namespace sidre */
} /* end namespace asctoolkit */


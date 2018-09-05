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
#ifndef MINT_MESH_HELPERS_HPP_
#define MINT_MESH_HELPERS_HPP_


#include "axom/core/Macros.hpp" // for AXOM_NOT_USED
#include "axom/core/Types.hpp"  // for nullptr

namespace axom
{
namespace mint
{
namespace internal
{

inline int dim( const double* AXOM_NOT_USED(x),
                const double* y,
                const double* z )
{
  return ( ( z != nullptr ) ? 3 : ( (y != nullptr ) ? 2 : 1 ) );
}


} /* namespace internal */
} /* namespace mint */
} /* namespace axom */

#endif /* MINT_MESH_HELPERS_HPP_ */

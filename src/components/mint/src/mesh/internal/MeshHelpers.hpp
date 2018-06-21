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


#include "axom/Macros.hpp" // for AXOM_NOT_USED
#include "axom/Types.hpp"  // for AXOM_NULLPTR

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
  return ( ( z != AXOM_NULLPTR ) ? 3 : ( (y != AXOM_NULLPTR ) ? 2 : 1 ) );
}

//------------------------------------------------------------------------------
inline int dim ( const IndexType& AXOM_NOT_USED( Ni ),
                 const IndexType& Nj,
                 const IndexType& Nk  )
{
  return ( (Nk >= 1) ? 3 : ( (Nj >= 1) ? 2 : 1 ) );
}

} /* namespace internal */
} /* namespace mint */
} /* namespace axom */

#endif /* MINT_MESH_HELPERS_HPP_ */

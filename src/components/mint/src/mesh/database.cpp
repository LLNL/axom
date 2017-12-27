/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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
#include "mint/database.hpp"

// Axom includes
#include "axom/Types.hpp" // for AXOM_NULLPTR
#include "slic/slic.hpp"  // for SLIC macros

#ifdef MINT_USE_SIDRE
#include "sidre/sidre.hpp"
#else
namespace axom
{
namespace sidre
{
class DataStore
{ };
class View
{ };
} /* namespace sidre */
} /* namespace axom */
#endif

namespace axom
{

static sidre::DataStore* s_database = AXOM_NULLPTR;

namespace mint
{
namespace database
{

//------------------------------------------------------------------------------
void set( sidre::DataStore* db )
{
#ifdef MINT_USE_SIDRE
  SLIC_ASSERT( db != AXOM_NULLPTR );
  s_database = db;
#else
  s_database = AXOM_NULLPTR;
  SLIC_WARNING( "Mint is not compiled with Sidre!" );
#endif
}

//------------------------------------------------------------------------------
sidre::View* get_view( const std::string& path )
{
#ifdef MINT_USE_SIDRE
  SLIC_ASSERT( s_database != AXOM_NULLPTR );
  return s_database->getRoot()->getView( path );
#else
  SLIC_WARNING( "Mint is not compiled with Sidre!" );
  return AXOM_NULLPTR;
#endif
}

//------------------------------------------------------------------------------
sidre::Group* get_group( const std::string& path )
{
#ifdef MINT_USE_SIDRE
  SLIC_ASSERT( s_database != AXOM_NULLPTR );
  return s_database->getRoot()->getGroup( path );
#else
  SLIC_WARNING( "Mint is not compiled with Sidre!" );
  return AXOM_NULLPTR;
#endif
}

} /* namespace database */
} /* namespace mint */
} /* namespace axom */

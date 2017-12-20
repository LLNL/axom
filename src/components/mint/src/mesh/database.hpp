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

#ifndef MINT_DATABASE_HPP_
#define MINT_DATABASE_HPP_

#ifndef MINT_USE_SIDRE
#define MINT_USE_SIDRE
#endif

#include <string> // for C++ string

namespace axom
{

// Forward Declarations
namespace sidre
{
  class View;
  class Group;
  class DataStore;
}

namespace mint
{

namespace database
{

/*!
 * \brief Sets the datastore instance that will be used.
 * \param [in] db pointer to the datastore
 * \pre db != AXOM_NULLPTR
 */
void set( sidre::DataStore* db );

/*!
 * \brief Returns a sidre::View instance along the given path.
 * \param [in] path the path of the requested view.
 * \return pointer to the view.
 * \note If the view does not exist, AXOM_NULLPTR is returned.
 */
sidre::View* get_view( const std::string& path );

/*!
 * \brief Returns a sidre::Group instance along the given path.
 * \param [in] path the path of the requested group.
 * \return pointer to the group.
 * \note If the group does not exist, AXOM_NULLPTR is returned.
 */
sidre::Group* get_group( const std::string& path );

} /* namespace database */
} /* namespace mint */
} /* namespace axom */


#endif /* MINT_DATABASE_HPP_ */

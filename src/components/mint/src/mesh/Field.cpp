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

#include "Field.hpp"

#include "mint/FieldTypes.hpp"
#include "axom/Types.hpp"

#include <cstddef>

namespace axom
{
namespace mint
{

Field::Field() :
  m_name(""),
  m_num_tuples(0),
  m_num_components(0),
  m_type(UNDEFINED_FIELD_TYPE)

{}

//------------------------------------------------------------------------------
Field::Field( const std::string& name,
              int size,
              int num_components ) :
  m_name( name ),
  m_num_tuples( size ),
  m_num_components( num_components ),
  m_type(UNDEFINED_FIELD_TYPE)
{}

//------------------------------------------------------------------------------
Field::~Field()
{}

//------------------------------------------------------------------------------
double* Field::getDoublePtr()
{
  return AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
const double* Field::getDoublePtr() const
{
  return AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
int* Field::getIntPtr()
{
  return AXOM_NULLPTR;
}

//------------------------------------------------------------------------------
const int* Field::getIntPtr() const
{
  return AXOM_NULLPTR;
}

} /* namespace mint */
} /* namespace axom */

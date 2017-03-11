/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "Field.hpp"

#include "mint/FieldTypes.hpp"
#include "axom/CommonTypes.hpp"

#include <cstddef>

namespace axom {
namespace mint {

Field::Field():
  m_name(""),
  m_num_tuples(0),
  m_num_components(0),
  m_type(UNDEFINED_FIELD_TYPE)

{}

//------------------------------------------------------------------------------
Field::Field( const std::string& name,
              int size,
              int num_components ):
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
int* Field::getIntPtr()
{
  return AXOM_NULLPTR;
}

} /* namespace mint */
} /* namespace axom */

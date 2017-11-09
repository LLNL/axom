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

#ifndef FIELDTYPES_HPP_
#define FIELDTYPES_HPP_

#include <string>

namespace axom
{
namespace mint
{

enum FieldTypes
{
  UNDEFINED_FIELD_TYPE =-1,

  DOUBLE_FIELD_TYPE,
  INTEGER_FIELD_TYPE,

  NUMBER_OF_FIELD_TYPES
};

template < typename FieldType >
struct field_of
{
  static const int type = UNDEFINED_FIELD_TYPE;
};


template < >
struct field_of< double >
{
  static const int type = DOUBLE_FIELD_TYPE;
};


template < >
struct field_of< int >
{
  static const int type = INTEGER_FIELD_TYPE;
};


} /* namespace mint */
} /* namespace axom */

#endif /* FIELDTYPES_HPP_ */

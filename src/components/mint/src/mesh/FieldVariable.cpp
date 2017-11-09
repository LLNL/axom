#include "FieldVariable.hpp"
#include "mint/Vector.hpp"


namespace axom
{
namespace mint
{


//------------------------------------------------------------------------------
template <>
double * FieldVariable< double >::getDoublePtr()
{
  return m_data.getData();
}

//------------------------------------------------------------------------------
template <>
const double * FieldVariable< double >::getDoublePtr() const
{
  return m_data.getData();
}

//------------------------------------------------------------------------------
template <>
int * FieldVariable< int >::getIntPtr()
{
  return m_data.getData();
}

//------------------------------------------------------------------------------
template <>
const int * FieldVariable< int >::getIntPtr() const
{
  return m_data.getData();
}


}   /* end namespace mint */
}   /* end namespace axom */

//
// DataUserType.cpp
//

#include "DataUserType.hpp"

namespace DataStore
{

void DataUserType::AddMember( const std::string& name )
{
    DataMember mem(name); 

    m_members.push_back(mem);
}

}

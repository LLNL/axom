/*
 * DataUserType.hpp
 */

#ifndef DATAUSERTYPE_HPP_
#define DATAUSERTYPE_HPP_

#include <string>
#include <vector>
#include "DataMember.hpp"

namespace DataStore
{

/**
 * \class DataUserType
 *
 * \brief Model a C struct, C++ class, Fortran derived type
 */
class DataUserType
{
public:
    /// non-callable default constructor
    DataUserType() = delete;

    DataUserType( const std::string& name ) :
	m_name(name)
    { }

    /**
     * \brief default destructor
     */
    ~DataUserType() {}

    inline int GetBPI()
    {
	return m_bpi;
    }

    void AddMember( const std::string& name );

#if 0
    const DataMembervoid GetMembers() const
    {
	return m_members;
    }
#endif

private:
    std::string m_name;
    int m_bpi;            // bytes per item (sum of members + alignment if struct)
    std::vector<DataMember> m_members;
};


}

#endif



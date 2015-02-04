/*
 * DataMember.hpp
 *
 * A DataMember models a member of DataType:
 *
 *  struct {
 *     int i1;
 *     int i2[10];
 *     int *i3;   // GetShape can compute length
 */

#ifndef DATAMEMBER_HPP_
#define DATAMEMBER_HPP_

#include "DatastoreInterface.hpp"

namespace DataStore
{

class DataUserType;

/**
 * \class DataMember
 *
 * \brief Model a member of a DataType
 */
class DataMember
{
public:
    /// non-callable default constructor
    DataMember() = delete;

    DataMember( const std::string& name ) :
	m_name(name)
    { }


    /**
     * \brief default destructor
     */
    ~DataMember() {}

    void GetShape();   // Shape may be computed from other members

private:
    std::string m_name;
    int m_offset;            // Offset to start of member
    DataUserType *m_type;
    DataShape *m_dataShape;
};


}

#endif



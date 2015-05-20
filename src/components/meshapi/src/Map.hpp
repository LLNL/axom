/**
 * \file Map.h
 *
 * \brief Basic API for a map from each element of a set to some domain
 *
 */

#ifndef MESHAPI_MAP_HPP_
#define MESHAPI_MAP_HPP_

#include <vector>
#include <sstream>
#include <iostream>

#include "common/Types.hpp"
#include "common/Utilities.hpp"
#include "meshapi/OrderedSet.hpp"

namespace asctoolkit {
namespace meshapi    {

    template<typename DataType>
    class Map
    {
    public:
        typedef MeshIndexType                                          Index;
        typedef MeshSizeType                                           size_type;

        typedef asctoolkit::meshapi::OrderedSet                        SetType;

        typedef std::vector<DataType>                                  OrderedMap;

    public:
        Map(SetType const* theSet = NULL) : m_set(theSet)
        {
            if(m_set) { m_data.resize( m_set->size() ); }
        }

        Map(SetType const* theSet, DataType defaultValue) : m_set(theSet)
        {
            if(m_set) { m_data.resize( m_set->size(), defaultValue ); }
        }

        ~Map(){}

        DataType const& operator[](Index setIndex) const
        {
            verifyIndex(setIndex);
            return m_data[setIndex];
        }

        DataType & operator[](Index setIndex)
        {
            verifyIndex(setIndex);
            return m_data[setIndex];
        }


        SetType const* set() const { return m_set; }


        //* Placeholder for function that returns the (pointer to) underlying data **/
        OrderedMap      & data()        { return m_data; }
        //* Placeholder for function that returns the (const pointer to) underlying data **/
        OrderedMap const& data() const  { return m_data; }


        size_type size() const { return m_set ? m_set->size() : size_type(); }

        bool isValid(bool verboseOutput = false) const;

    private:
        inline void  verifyIndex(Index setIndex)       const { ATK_ASSERT( m_set && (setIndex < m_set->size() ) ); }

    private:
        SetType const*  m_set;
        OrderedMap         m_data;
    };




    template<typename DataType>
    bool Map<DataType>::isValid(bool verboseOutput) const
    {
        bool bValid = true;

        std::stringstream errStr;

        if(!m_set)
        {
            if(! m_data.empty() )
            {
                if(verboseOutput)
                {
                    errStr << "\n\t* the underlying set was never set, but its associated data is not empty"
                        <<" , data has size " << m_data.size();
                }

                bValid = false;
            }
        }
        else
        {
            if( m_data.size() != m_set->size())
            {
                if(verboseOutput)
                {
                    errStr << "\n\t* the underlying set and its associated mapped data have different sizes"
                        <<" , underlying set has size " << m_set->size()
                        <<" , data has size " << m_data.size();
                    ;
                }

                bValid = false;
            }
        }


        if(verboseOutput)
        {
            std::stringstream sstr;

            sstr<<"\n*** Detailed results of isValid on the map.\n";
            if(bValid)
            {
                sstr<<"Map was valid."<< std::endl;
            }
            else
            {
                sstr<<"Map was NOT valid.\n"
                         << sstr.str()
                         << std::endl;
            }

            if(!m_set)
            {
                sstr<<"\n** map is empty.";
            }
            else
            {
                sstr<< "\n** underlying set has size " << m_set->size() <<": ";

                sstr<< "\n** Mapped data:";
                for(Index idx = 0; idx < this->size(); ++idx)
                {
                    sstr<<"\n\telt[" << idx << "]:\t" << (*this)[idx];
                }
            }
            std::cout<< sstr.str() << std::endl;
        }

        return bValid;
    }




} // end namespace meshapi
} // end namespace asctoolkit



#endif // MESHAPI_MAP_HPP_

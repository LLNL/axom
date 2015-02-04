/*
 * DataFunction.hpp
 */

#ifndef DATAFUNCTION_HPP_
#define DATAFUNCTION_HPP_

namespace DataStore
{

/**
 * \class DataFunction
 *
 * \brief Abstract class to provide interface for functions.
 */
class DataFunction
{
public:
    /// non-callable default constructor
    //    DataFunction() = delete;

    /**
     * \brief default destructor
     */
    //    virtual ~DataFunction();  // XXX undefined reference to `vtable for DataStore::DataFunction'
    virtual ~DataFunction() {};


    virtual void AddArgument();
    virtual void Call() = 0;
};


}

#endif



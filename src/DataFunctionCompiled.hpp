/*
 * DataFunctionCompiled.hpp
 */

#ifndef DATAFUNCTIONCOMPILED_HPP_
#define DATAFUNCTIONCOMPILED_HPP_

#include "DataFunction.hpp"

namespace DataStore
{

/**
 * \class DataFunctionCompiled
 *
 * \brief Concrete class to provide interface for compiled functions.
 */
class DataFunctionCompiled : DataFunction
{
public:
    DataFunctionCompiled( void *(*fcn)()) :
	m_fcn(fcn)
    {}
    //    ~DataFunctionCompiled() {};


    void AddArgument() override;
    void Call() override;

private:
    void *(*m_fcn)(void);
};


}

#endif



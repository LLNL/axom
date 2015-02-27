/*
 * DataFunctionPython.hpp
 */

#ifndef DATAFUNCTIONPYTHON_HPP_
#define DATAFUNCTIONPYTHON_HPP_

#include "DataFunction.hpp"
#include "Python.h"

namespace DataStore
{

/**
 * \class DataFunctionPython
 *
 * \brief Concrete class to provide interface for python functions.
 */
class DataFunctionPython : DataFunction
{
public:
    DataFunctionPython( PyObject *fcn ) :
	m_fcn(fcn)
    {
	//int PyCallable_Check(PyObject *o)
	Py_INCREF(fcn);
    }

    ~DataFunctionPython()
    {
	Py_DECREF(m_fcn);
    };

    void AddArgument() override;
    void Call() override;

private:
    PyObject *m_fcn;
};


}

#endif



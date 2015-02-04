/*
 * DataFunctionPython.hpp
 */

#include "DataFunctionPython.hpp"

namespace DataStore
{

#if 0
void DataFunctionPython::AddArgument()
{
    return;
}
#endif

void DataFunctionPython::Call()
{
    PyObject *args = PyTuple_New(0);
    PyObject *kw = NULL;
    PyObject* rv = PyObject_Call(m_fcn, args, kw);
    Py_DECREF(args);
    return;
}

}

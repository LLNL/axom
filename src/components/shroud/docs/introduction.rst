Introduction
============

Shroud is a tool for exposing a C++ API to Fortran.
It does this by using C as the lingua franca.
C++ can communicate with C via a common heritage and the ``extern "C"`` keyword.
Fortran uses the interoperability with C features introduced in Fortran 2003
including the ``iso_c_binding`` module and the ``bind`` and ``value`` keywords.
A C API for C++ API is produced as a by product of the Fortran wrapping.

Goals
-----

  * Simplify the creating of wrapper for Fortran.
  * Preserves the object-oriented style of C++ classes.
  * Create a Fortran idiomatic API from the C++ API.


Using a C++ API to create an object and call a method::

    Instance * inst = new Instance;
    inst->method(1);

In Fortran this becomes::

    type(instance) inst
    inst = instance_new()
    call inst%method(1)

.. note :: The ability to generate C++ wrappers for Fortran is not supported.

Issues
------

There is a long history of ad-hoc solutions to provide C and Fortran interoperability.
Any solution must address several problems:

  * Name mangling of externals.  This includes namespaces and operator overloading in C++.
  * Call-by-reference vs call-by-value differences
  * Length of string arguments.
  * Blank filled vs null terminated strings.

The 2003 Fortran standard added several features for interoperability with C:

  * iso_c_binding - intrinsic module which defines fortran kinds for matching with C's types.
  * ``BIND`` keyword to control name mangling of externals.
  * ``VALUE`` attribute to allow pass-by-value.

In addition, Fortran 2003 provides some object oriented programming facilities:

   * Type extension
   * Procedure Polymorphism with Type-Bound Procedures
   * Enumerations compatible with C

A Fortran pointer is similar to to a C++ class.  It not only has the address of 
the memory but also contains meta-data such as the type, kind and shape of the array.

   * GNU Fortran http://gcc.gnu.org/wiki/ArrayDescriptorUpdate
   * Intel 15.0 https://software.intel.com/en-us/node/525356

Fortran provides a **pointer** and **allocatable** attributes which are not
directly supported by C.  Each vendor has their own pointer struct.
The Chasm library from LANL can be used to write portable access to Fortran pointers.
http://chasm-interop.sourceforge.net
Eventually this will be supported in Fortran via the Further Interoperability of Fortran and C - Technical Specification ISO/IEC TS 29113:2012
http://www.iso.org/iso/iso_catalogue/catalogue_tc/catalogue_detail.htm?csnumber=45136


Requirements
------------

Fortran wrappers are generated as free-form source and require a Fortran 2003 compiler.

Shroud is written in Python and requires the modules:
  * PyYAML https://pypi.python.org/pypi/PyYAML/3.11
  * Parsley https://pypi.python.org/pypi/Parsley

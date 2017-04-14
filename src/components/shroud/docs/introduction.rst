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

A Fortran pointer is similar to a C++ instance in that it not only has the address of 
the memory but also contains meta-data such as the type, kind and shape of the array.  Some vendors document the struct used to store the metadata for an array.

   * GNU Fortran http://gcc.gnu.org/wiki/ArrayDescriptorUpdate
   * Intel 17 https://software.intel.com/en-us/node/678452

..   * Intel 15.0 https://software.intel.com/en-us/node/525356

Fortran provides a **pointer** and **allocatable** attributes which are not
directly supported by C.  Each vendor has their own pointer struct.
Eventually this will be supported in Fortran via the Further Interoperability of Fortran and C -
`Technical Specification ISO/IEC TS 29113:2012 <http://www.iso.org/iso/iso_catalogue/catalogue_tc/catalogue_detail.htm?csnumber=45136>`_

Overview
--------

Input is read from a YAML file which describes the types, functions,
and classes to wrap.  This file must be created by the user.  Shroud
does not parse C++ code to extract the API. That was considered a
large task and not needed for the size of the API of the library that
inspired Shroud's development. In addition, there is a lot of semantic
information which must be provided by the user that may be difficult
to infer from the source alone.  However, the task of created the
input file is simplified since the C++ declaration can be
cut-and-pasted into the YAML file.

In some sense, Shroud can be thought of as a fancy macro processor.
It takes the function declarations from the YAML file and break them
down into a series of contexts (library, class, function, argument) and defines
a dictionary of format macros of the form key=value.  There are then a
series of macro templates which are expanded to create the wrapper
functions. Some name templates can be specified as options.  But the
overall structure of the generated code is defined by the classes and
functions in the YAML file as well as the requirements of C++ and
Fortran syntax.

Requirements
------------

Fortran wrappers are generated as free-form source and require a Fortran 2003 compiler.

Shroud is written in Python and has been tested with version 2.7 and 3.6.
It requires the modules:

  * PyYAML https://pypi.python.org/pypi/PyYAML/3.11
  * Parsley https://pypi.python.org/pypi/Parsley


XKCD
----

`XKCD <https://xkcd.com/1319>`_

.. image:: automation.png

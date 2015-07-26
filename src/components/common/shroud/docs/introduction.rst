Introduction
============

Shroud is a tool for exposing a C++ API to Fortran.
It does this by using C as the lingua franca.
C++ can communicate with C via a common heritage and the **extern "C"** keyword.
Fortran uses the interoperability with C features introduced in Fortran 2003
including the **iso_c_binding** module and the **bind** keyword.
A C API for C++ API is produced as a by produced of the Fortran wrapping.

The Fortran wrappers preserves the object-oriented style using the
object-oriented features of Fortran 2003.

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

  * Name mangling of externals
  * Call-by-reference vs call-by-value differences
  * Length of string arguments.
  * Blank filled vs null terminated strings.

The 2003 Fortran standard added several features for interoperability with C:

  * iso_c_binding - intrinsic module which defines fortran kinds for matching with C's types.
  * BIND keyword to control name mangling of externals.
  * VALUE attribute to allow pass-by-value.

In addition, Fortran 2003 provides some object oriented programming facilities:

   * Type extension
   * Procedure Polymorphism with Type-Bound Procedures
   * Enumerations compatible with C

A Fortran pointer is similar to to a C++ class.  It not only has the address of 
the memory but also contains meta-data such as the type, kind and shape of the array.

Fortran provides a **pointer** and **allocatable** attributes which are not
directly supported by C.  Each vendor has their own pointer struct.
The Chasm library from LANL can be used to write portable access to Fortran pointers.


Requirements
------------

Shroud is written in Python an requires the modules:
  * PyYAML https://pypi.python.org/pypi/PyYAML/3.11
  * Parsley https://pypi.python.org/pypi/Parsley

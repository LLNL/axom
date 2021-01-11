.. _inlet_advanced_types_label:

Advanced Types
==============

In addition to Inlet's primitive types (bool, int, double, string), user-defined types
and functions can also be defined as part of an input file.

Adding User-Defined Types to a Schema
-------------------------------------

To add a single (i.e., not array) user-defined type to the input file, use the ``addTable``
function of the Inlet or Table classes to add a Table (collection of Fields and sub-Tables)
that will represent the fields of the struct.

Consider a simple struct that contains only primitive types, whose definition in Lua might look like:

.. code-block:: Lua

  mesh = {
    filename = "/data/star.mesh",
    serial = 1,
    parallel = 1
  }

Its Inlet schema can be defined as follows:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_simple_start
   :end-before: _inlet_userdef_simple_end
   :language: C++

This would be used by defining a table for the ``Mesh`` instance and then defining the struct
schema on that subtable, e.g.:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_simple_usage_start
   :end-before: _inlet_userdef_simple_usage_end
   :language: C++

The definition of a static ``defineSchema`` member function is not required, and is just used
for convenience.  The schema definition for a class or struct could also be implemented as a
free function for third-party types, or even in the same place as the sub-table declaration. 

Arrays of user-defined types are also supported in Inlet.  First, use the ``addGenericArray``
function to create a subtable, then define the schema on that table:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_array_usage_start
   :end-before: _inlet_userdef_array_usage_end
   :language: C++

Retrieving User-Defined Types from an Input File
------------------------------------------------

In order to have Inlet extract from the input file into a struct, it needs to know how that struct can be built.
This is accomplished by a specializing of the struct ``FromInlet<T>`` for your type ``T``.  It must define a
single member function with the signature ``T operator()(const inlet::Table&)``.  This function should return
an instance of type ``T`` with its members populated with the corresponding fields in the input file.
For the simple ``Mesh`` example whose schema is defined above, the specialization might look like:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_simple_frominlet_start
   :end-before: _inlet_userdef_simple_frominlet_end
   :language: C++

The ``Table::operator[]`` is used to extract fields from the input file and is automatically converted to
the correct member variable types when the function's return value is constructed.  This conversion does
not happen automatically for user-defined types.  If a ``Mesh`` object as defined above is located at the
path "mesh" within the input file, it can be retrieved as follows:

.. code-block:: C++

  Mesh mesh = inlet["mesh"].get<Mesh>();

As with the schema definition, the ``FromInlet`` specialization for a user-defined type will work for both
single instances of the type and arrays of the type.

Arrays in Inlet are implemented as ``std::unordered_map<int, T>``, so retrieving an array of user-defined
types is as follows, in this case for an array of the ``BoundaryCondition`` struct:

.. code-block:: C++

  auto mesh = inlet["bcs"].get<std::unordered_map<int, BoundaryCondition>>();

Adding Functions to a Schema
----------------------------

For input file types that support functions, e.g., Lua, functions can also be read from the input file
into a ``std::function``, the wrapper for callables provided by the C++ standard library.  This is accomplished
by calling ``addFunction`` on an Inlet or Table object.

Consider the following Lua function that accepts a three-dimensional vector (split into components x, y, and z)
and returns a double:

.. code-block:: Lua

  coef = function (x, y, z)
    return x + (y * 0.5) + (z * 0.25)
  end

The schema for this function would be defined as follows:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_func_coef_start
   :end-before: _inlet_userdef_func_coef_end
   :language: C++

Note that a single type tag is passed for the return type, while a vector of tags is passed
for the argument types.  Currently a maximum of two arguments are supported, with possible argument
types ``Double`` or ``Vec3D``.  These correspond to the C++ types ``double`` and
``axom::primal::Vector3D``, respectively.

.. note::  The function retrieval implementation for Lua will automatically expand vector arguments into three
  arguments and contract vector returns from three scalars.  That is, a function whose Inlet schema contains a ``Vec3D``
  argument should accept three scalar arguments in its place, and a function whose Inlet schema contains a ``Vec3D``
  return value should return three scalar arguments in its place.

Retrieving Functions from an Input File
---------------------------------------

To retrieve a function, both the implicit conversion and ``get<T>`` syntax is supported.  For example,
a function can be retrieved as follows:

.. literalinclude:: ../../examples/mfem_coefficient.cpp
   :start-after: _inlet_mfem_coef_simple_retrieve_start
   :end-before: _inlet_mfem_coef_simple_retrieve_end
   :language: C++

It can also be assigned directly to a ``std::function`` without the need to use ``get<T>``:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_bc_struct_start
   :end-before: _inlet_userdef_bc_struct_end
   :language: C++

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_bc_struct_retrieve_start
   :end-before: _inlet_userdef_bc_struct_retrieve_end
   :language: C++

Additionally, if a function does not need to be stored, the overhead of a copy can be eliminated
by calling it directly:

.. code-block:: C++

  double result = inlet["coef"].call<double>(axom::inlet::FunctionType::Vec3D{3, 5, 7});

.. note::  Using ``call<ReturnType>(ArgType1, ArgType2, ...)`` requires both that the return type
  be explicitly specified and that argument types be passed with the exact type as used in the 
  signature defined as part of the schema.  This is because the arguments do not participate in
  overload resolution.

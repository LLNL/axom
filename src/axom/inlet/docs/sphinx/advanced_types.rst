.. _inlet_advanced_types_label:

##############
Advanced Types
##############

In addition to Inlet's primitive types (bool, int, double, string), user-defined types
and functions can also be defined as part of an input file.

*****************
Individal Structs
*****************

Defining And Storing
--------------------

To add a single (i.e., not array) user-defined type to the input file, use the ``addStruct``
function of the Inlet or Table classes to add a Table (collection of Fields and sub-Tables)
that will represent the fields of the struct.

Consider a simple Lua table that contains only primitive types, whose definition might look like:

.. code-block:: Lua

  car = {
      make = "BestCompany",
      seats = 2,
      horsepower = 200,
      color = "blue"
  }

or in YAML, something like:

.. code-block:: YAML
  
  car:
    make: BestCompany
    seats: 2
    horsepower: 200
    color: blue

Its Inlet schema can be defined as follows:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_simple_start
   :end-before: _inlet_userdef_simple_end
   :language: C++

This would be used by creating an ``inlet::Table`` for the ``Car`` instance and then defining the struct
schema on that subtable, e.g.:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_simple_usage_start
   :end-before: _inlet_userdef_simple_usage_end
   :language: C++

.. note:: 
  The definition of a static ``defineSchema`` member function is not required, and is just used
  for convenience.  The schema definition for a class or struct could also be implemented as a
  free function for third-party types, or even in the same place as the sub-table declaration.

Accessing
---------

In order to convert from Inlet's internal representation to a custom C++ ``struct``, you must provide 
deserialization logic.  This is accomplished by a specialization of the ``FromInlet<T>`` functor for your
type ``T``, which implements a member function with the signature ``T operator()(const inlet::Table&)``.
This function should return an instance of type ``T`` with its members populated with the corresponding
fields in the input file. For the simple ``Car`` example whose schema is defined above,the specialization
might look like:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_simple_frominlet_start
   :end-before: _inlet_userdef_simple_frominlet_end
   :language: C++

In the above snippet, ``Table::operator[]`` is used to extract data from Inlet's internal representation
which is automatically converted to the correct member variable types when the function's return value is
constructed.  This conversion does not happen automatically for user-defined types.
If a ``Car`` object as defined above is located at the path "car" within the input file, it can be retrieved as follows:

.. code-block:: C++

  Car car = inlet["car"].get<Car>();

**********************************
Arrays and Dictionaries of Structs
**********************************

Arrays of user-defined types are also supported in Inlet.  

Defining And Storing
--------------------

Consider a collection of cars, described in Lua as:

.. code-block:: Lua

  fleet = {
    {
      make = "Globex Corp",
      seats = 3,
      horsepower = 155,
      color = "green"
    },
    {
      make = "Initech",
      seats = 4,
      horsepower = 370
      -- uses default value for color
    },
    {
      make = "Cyberdyne",
      seats = 1,
      horsepower = 101,
      color = "silver"
    }
  }

or in YAML, as:

.. code-block:: YAML

  fleet:
    - make: Globex Corp
      seats: 3
      horsepower: 155
      color: green
    - make: Initech
      seats: 4
      horsepower: 370
      # uses default value for color
    - make: Cyberdyne
      seats: 1
      horsepower: 101
      color: silver

First, use the ``addStructArray`` function to create a subtable, then define the schema on that table:

.. literalinclude:: ../../examples/user_defined_type.cpp
   :start-after: _inlet_userdef_collection_usage_start
   :end-before: _inlet_userdef_collection_usage_end
   :language: C++

.. note::
  The schema definition logic for a struct is identical between individual instances of structs
  and arrays of structs.  The distinction is made by ``Table`` on which the struct schema is
  defined - specifically, whether it is obtained via ``addStruct`` or ``addStructArray``.

Associative arrays are also supported, using string keys or a mixture of string and integer keys.
The ``addStructDictionary`` function can be used analogously to the ``addStructArray`` function
for these associative arrays.

.. note::
  Although many of Inlet's input file languages do not distinguish between a "dictionary" type
  and a "record" type, Inlet treats them differently for type safety reasons:

  *Dictionaries* use arbitrary strings or integers for their keys, and their values (entries)
  can only be retreived as a homogenous type.  In other words, dictionaries must map to
  ``std::unordered_map<Key, Value>`` for fixed key and value types.

  *Structs* contain a fixed set of named fields, but these fields can be of any type.
  As the name suggests, these map to ``structs`` in C++.

Accessing
---------

As with the schema definition, the ``FromInlet`` specialization for a user-defined type will work for both
single instances of the type and arrays of the type.

To retrieve an array of structs as a contiguous array of user-defined type, use ``std::vector``:

.. code-block:: C++

  auto fleet = base["fleet"].get<std::vector<Car>>();

Some input file languages support non-contiguous array indexing, so you can also retrieve arrays 
as ``std::unordered_map<int, T>``:

.. code-block:: C++

  auto fleet = inlet["fleet"].get<std::unordered_map<int, Car>>();

.. note::
  If a non-contiguous array is retrieved as a (contiguous) ``std::vector``, the elements will be
  ordered by increasing index.

String-keyed dictionaries are implemented as ``std::unordered_map<std::string, T>`` and can be retrieved
in the same way as the array above.  For dictionaries with a mix of string and integer keys, the
``inlet::VariantKey`` type can be used, namely, by retrieving a ``std::unordered_map<inlet::VariantKey, T>``.

*********
Functions
*********

For input file types that support functions, e.g., Lua, functions can also be read from the input file
into a ``std::function``, the wrapper for callables provided by the C++ standard library.

Defining And Storing
--------------------

This is accomplished by calling ``addFunction`` on an Inlet or Table object.

Consider the following Lua function that accepts a vector in **R**\ :sup:`2` or **R**\ :sup:`3` and returns a double:

.. code-block:: Lua

  coef = function (v)
    if v.dim == 2 then
      return v.x + (v.y * 0.5)
    else
      return v.x + (v.y * 0.5) + (v.z * 0.25)
    end
  end

The schema for this function would be defined as follows:

.. literalinclude:: ../../examples/mfem_coefficient.cpp
   :start-after: _inlet_mfem_func_coef_start
   :end-before: _inlet_mfem_func_coef_end
   :language: C++

The return type and argument types are described with the ``inlet::FunctionTag`` enumeration, which has the following members:

  * ``Double`` - corresponds to a C++ ``double``
  * ``String`` - corresponds to a C++ ``std::string``
  * ``Vector`` - corresponds to a C++ ``inlet::InletVector``
  * ``Void`` - corresponds to C++ ``void``, should only be used for functions that don't return a value

Note that a single type tag is passed for the return type, while a vector of tags is passed
for the argument types.  Currently a maximum of two arguments are supported. 
To declare a function with no arguments, simply leave the list of argument types empty.

.. note::  The ``InletVector`` type (and its Lua representation) are statically-sized vectors with
  a maximum dimension of three.  That is, they can also be used to represent two-dimensional vectors.

In Lua, the following operations on the ``Vector`` type are supported (for ``Vector`` s ``u``, ``v``, and ``w``):

1. Construction of a 3D vector: ``u = Vector.new(1, 2, 3)``
#. Construction of a 2D vector: ``u = Vector.new(1, 2)``
#. Construction of an empty vector (default dimension is 3): ``u = Vector.new()``
#. Vector addition and subtraction: ``w = u + v``, ``w = u - v``
#. Vector negation: ``v = -u``
#. Scalar multiplication: ``v = u * 0.5``, ``v = 0.5 * u``
#. Indexing (1-indexed for consistency with Lua): ``d = u[1]``, ``u[1] = 0.5``
#. L2 norm and its square: ``d = u:norm()``, ``d = u:squared_norm()``
#. Normalization: ``v = u:unitVector()``
#. Dot and cross products: ``d = u:dot(v)``, ``w = u:cross(v)``
#. Dimension retrieval: ``d = u.dim``
#. Component retrieval: ``d = u.x``, ``d = u.y``, ``d = u.z``

Accessing
---------

To retrieve a function, both the implicit conversion and ``get<T>`` syntax is supported.  For example,
a function can be retrieved as follows:

.. literalinclude:: ../../examples/mfem_coefficient.cpp
   :start-after: _inlet_mfem_coef_simple_retrieve_start
   :end-before: _inlet_mfem_coef_simple_retrieve_end
   :language: C++

It can also be assigned directly to a ``std::function`` without the need to use ``get<T>``:

.. code-block:: C++
  std::function<double(FunctionType::Vector)> coef = inlet["coef"];

Additionally, if a function does not need to be stored, the overhead of a copy can be eliminated
by calling it directly:

.. code-block:: C++

  double result = inlet["coef"].call<double>(axom::inlet::FunctionType::Vector{3, 5, 7});

.. note::  Using ``call<ReturnType>(ArgType1, ArgType2, ...)`` requires both that the return type
  be explicitly specified and that argument types be passed with the exact type as used in the 
  signature defined as part of the schema.  This is because the arguments do not participate in
  overload resolution.

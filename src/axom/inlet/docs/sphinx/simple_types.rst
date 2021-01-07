############
Simple Types
############

Inlet supports four simple types in the following groups: scalars, tables, arrays,
and dictionaries. Inlet refers to individual items as Fields. Tables can contain multiple
Fields, a single array, and a single dictionar.

.. note::  There is an inherent global table that holds all global level Fields.  This can be
  accessed via your `Inlet` class instance.


*******
Scalars
*******

There are four supported scalar types: `bool`, `int`, `double`, and `string`.

In this example we will be using the following part of an input file:

.. literalinclude:: ../../examples/simple_types.lua
   :start-after: _inlet_simple_types_scalar_input_start
   :end-before: _inlet_simple_types_scalar_input_end
   :language: lua


Defining And Storing
--------------------

This example shows how to add the four simple scalar types with descriptions to the
input file schema and add their values, if present in the input file, to the Sidre
DataStore to be accessed later.

.. literalinclude:: ../../examples/simple_types.cpp
   :start-after: _inlet_simple_types_scalar_add_start
   :end-before: _inlet_simple_types_scalar_add_end
   :language: C++

You can also add default values to Fields to fall back to if they are not defined
in your input file. The last added Field was intentionally not present in the input file.  Not all
fields need to be present, unless they are marked required, like `a_simple_double`.

Accessing
---------

Accessing field values stored in Inlet can be accessed via their name with the `[]` operator.

.. literalinclude:: ../../examples/simple_types.cpp
   :start-after: _inlet_simple_types_scalar_access_start
   :end-before: _inlet_simple_types_scalar_access_end
   :language: C++

Here is an example of how to check if the input file contained an optional field:

.. literalinclude:: ../../examples/simple_types.cpp
   :start-after: _inlet_simple_types_scalar_contains_start
   :end-before: _inlet_simple_types_scalar_contains_end
   :language: C++


******
Tables
******

Coming soon!

Defining And Storing
--------------------

Coming soon!

Accessing
---------

Coming soon!

******
Arrays
******

Coming soon!

Defining And Storing
--------------------

Coming soon!

Accessing
---------

Coming soon!

************
Dictionaries
************

Coming soon!

Defining And Storing
--------------------

Coming soon!

Accessing
---------

Coming soon!

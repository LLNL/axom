############
Simple Types
############

To help structure your input file, Inlet categories the information into
two types: Fields and Tables.

Fields refer to the individual scalar values that are either at the global level or
that is contained inside of a Table.

Tables can contain multiple Fields, as well as, a single array, and a single dictionary.

.. note::  There is an inherent global Table that holds all global level Fields.  This can be
  accessed via your `Inlet` class instance.


*******
Fields
*******

There are four supported field types: `bool`, `int`, `double`, and `string`.

In this example we will be using the following part of an input file:

.. literalinclude:: ../../examples/fields.cpp
   :start-after: _inlet_simple_types_fields_input_start
   :end-before: _inlet_simple_types_fields_input_end
   :language: lua


Defining And Storing
--------------------

This example shows how to add the four simple field types with descriptions to the
input file schema and add their values, if present in the input file, to the Sidre
DataStore to be accessed later.

.. literalinclude:: ../../examples/fields.cpp
   :start-after: _inlet_simple_types_fields_add_start
   :end-before: _inlet_simple_types_fields_add_end
   :language: C++

You can also add default values to Fields to fall back to if they are not defined
in your input file. The last added Field was intentionally not present in the input file.  Not all
fields need to be present, unless they are marked required, like `a_simple_double`.

Accessing
---------

Accessing field values stored in Inlet can be accessed via their name with the `[]` operator.

.. literalinclude:: ../../examples/fields.cpp
   :start-after: _inlet_simple_types_fields_access_start
   :end-before: _inlet_simple_types_fields_access_end
   :language: C++

Here is an example of how to check if the input file contained an optional field:

.. literalinclude:: ../../examples/fields.cpp
   :start-after: _inlet_simple_types_fields_contains_start
   :end-before: _inlet_simple_types_fields_contains_end
   :language: C++


******
Tables
******

Tables help with grouping associated data together into a single container. Tables can
contain multiple individually named Fields, as well as, a single array, and a single
dictionary.    

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

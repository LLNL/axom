############
Simple Types
############

To help structure your input file, Inlet categorizes the information into
two types: Fields and Containers.

Fields refer to the individual scalar values that are either at the global level or
that are contained inside of a Container.

Containers can contain multiple Fields, other sub-Containers, as well as a single array or a single dictionary.

.. note::  There is a global Container that holds all top-level Fields.  This can be
  accessed via your `Inlet` class instance.


*******
Fields
*******

In Inlet, Fields represent an individual scalar value of primitive type.
There are four supported field types: ``bool``, ``int``, ``double``, and ``string``.

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
fields need to be present, unless they are marked required, like ``a_simple_double``.

Accessing
---------

Accessing field values stored in Inlet can be accessed via their name with the ``[]`` operator
or through the templated ``get<T>`` function.  The ``[]`` operator is more streamlined but
can lead to compile time ambiquity depending on how it is used.  The example below shows
an example of this.

Prior to accessing optional fields, you should verify they were provided by the user via
the ``contains`` function.

The ``contains`` function returns ``true`` if the field was *either*
provided by the user or via a default.  To check if the field was provided by the user (and not
via a default), you can use the ``isUserProvided`` method, which returns ``true`` if the value
was provided by the user in the input file.

Accessing a value that was not provided by the user, or 
a default value, will result in a runtime error.

.. literalinclude:: ../../examples/fields.cpp
   :start-after: _inlet_simple_types_fields_access_start
   :end-before: _inlet_simple_types_fields_access_end
   :language: C++

.. note:: The field ``does_not_exist`` was purposefully left this out of the
   user-provided input file to show no warnings/errors are thrown during runtime
   for defining optional fields in the schema.


**********
Containers
**********

Containers help with grouping associated data together into a single collection. Containers can
contain multiple individually named Fields, multiple sub-Containers, as well as a single
array or a single dictionary.

In this example, we will be using the following part of an input file:

.. literalinclude:: ../../examples/containers.cpp
   :start-after: _inlet_simple_types_containers_input_start
   :end-before: _inlet_simple_types_containers_input_end
   :language: lua


Defining And Storing
--------------------

This example shows how to add a Container with a nested Container to the
input file schema and add the underlying field values to the Sidre
DataStore to be accessed later.

.. literalinclude:: ../../examples/containers.cpp
   :start-after: _inlet_simple_types_containers_add_start
   :end-before: _inlet_simple_types_containers_add_end
   :language: C++

This example also shows that the ``color`` Field that was not given in the
input file but used the default value that was specified in the schema.

.. note:: Inlet also has an ``addStruct`` member for defining more complex types,
   such as nested structures. See :ref:`Advanced Types <inlet_advanced_types_label>`
   for more details


Accessing
---------

Field values stored inside a container can be accessed via their name with the ``[]`` operator.
They can be accessed from the Inlet class instance with their fully qualified name or you
can get the Container instance first, then access it with the relative name.

.. literalinclude:: ../../examples/containers.cpp
   :start-after: _inlet_simple_types_containers_access_start
   :end-before: _inlet_simple_types_containers_access_end
   :language: C++


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

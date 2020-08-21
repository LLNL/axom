.. _inlet_quick_start_label:

Quick Start
===========

Inlet's workflow is broken into the three following steps:  

* Defining the schema of your input deck
* Verifying that the user given input deck is valid
* Accessing the read and stored data in your program

.. _inlet_defining_schema_label:

Defining Schema
---------------

The first step in using Inlet is to defining the schema of your input deck.
Inlet defines an input deck into two basic classes: Tables and Fields. Basically
Fields are individual values and Tables hold groups of Fields and Tables.

Define the schema by using the following functions, on either the main Inlet class, for
global Tables and Fields, or on individual Table classes, for Tables and Fields under that Table:

========================= ===================
Name                      Description
========================= ===================
addTable                  Adds a Table to the input deck schema with the given name.
addBool                   Adds a boolean Field to the global or parent Table with the given name.
addDouble                 Adds a double Field to the global or parent Table with the given name.
addInt                    Adds a integer Field to the global or parent Table with the given name.
addString                 Adds a string Field to the global or parent Table with the given name.
========================= ===================

All possible Tables and Fields that are can be found in the input deck must be defined
at this step.  The value of the Field is read and stored into the Sidre datastore when you call the appropriate
add function. Use the ``required`` class member function on the Table and Field class to indicate that
they have to present in the given input deck. You can also set a default value to each field via the type-safe
``Field::defaultValue()`` member functions. Doing so will populate the corresponding Fields value
if the specific Field is not present in the input deck. The following example shows these concepts:

.. literalinclude:: ../../examples/verification.cpp
   :start-after: _inlet_workflow_defining_schema_start
   :end-before: _inlet_workflow_defining_schema_end
   :language: C++


.. _inlet_verification_label:

Verification
------------

This step helps ensure that the given input deck follows the rules expected by the code.  These
rules are not verified until you call ``Inlet::verify()``.  Doing so will return true/false and
output SLIC warnings to indicate which Field or Table violated which rule.

As shown above, both Tables and Fields can be marked as ``required``. Fields have two additional
basic rules that can be enforced with the following ``Field`` class member functions:

========================= ===================
Name                      Description
========================= ===================
validValues               Indicates the Field can only be set to a valid values.
range                     Indicates the Field can only be set to inclusively between two values.
========================= ===================

Inlet also provides functionality to write your own custom rules via callable lambda verifiers.
Fields and Tables can both register one lambda each via their ``registerVerifier()`` member functions.
The following example adds a custom verifier that simply verifies that the dimension of the simulation
match up with the dimensions of a given vector:

.. literalinclude:: ../../examples/verification.cpp
   :start-after: _inlet_workflow_verification_start
   :end-before: _inlet_workflow_verification_end
   :language: C++

.. note::  ``Inlet::getGlobalTable()->registerVerifier()`` can be used to add a verifier to apply rules
  to the Fields at the global level.


.. _inlet_accessing_data_label:

Accessing Data
--------------

After the input deck has been read and verified by the previous steps.  You can access the data by name
via ``Inlet::get()`` functions.  These functions are type-safe, fill the given variable with what is found,
and return a boolean whether the Field was present in the input deck or had a default value it could fall
back on.  Variable are named in a language agnostic way and are converted from Inlet's representation 
to the language specific version inside of the appropriate ``Reader``. For example, Inlet refers to the
Lua variable ``vector={x=3}`` as ``vector/x``.

For example, given the previous verificiation example, this access previously read values:

.. literalinclude:: ../../examples/verification.cpp
   :start-after: _inlet_workflow_accessing_data_start
   :end-before: _inlet_workflow_accessing_data_end
   :language: C++

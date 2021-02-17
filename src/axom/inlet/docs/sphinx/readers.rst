#######
Readers
#######

Inlet has built-in support for three input file languages: JSON, Lua, and YAML.
Due to language features, not all readers support all Inlet features.
Below is a container that lists supported features:


.. list-container:: Supported Language Features
   :header-rows: 1

   * - 
     - JSON
     - Lua
     - YAML
   * - Primitive Types
     - bool, double, int, string
     - bool, double, int, string
     - bool, double, int, string
   * - Dictionaries
     - X
     - X
     - X
   * - Arrays
     - X
     - X
     - X
   * - Non-contiguous Arrays
     - 
     - X
     - 
   * - Mixed-typed key Arrays
     - 
     - X
     - 
   * - Callback Functions
     - 
     - X
     - 

***********************
Extra Lua functionality
***********************

Inlet opens four Lua libraries by default: ``base``, ``math``, ``string``, ``package``.  All libraries are
documented `here <https://sol2.readthedocs.io/en/v2.20.6/api/state.html?highlight=open_libraries#enumerations>`_. 

For example, you can add the `io` library by doing this:

.. literalinclude:: ../../examples/lua_library.cpp
   :start-after: _inlet_io_library_add_start
   :end-before: _inlet_io_library_add_end
   :language: C++

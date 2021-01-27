
Inlet User Guide
================

.. note:: Inlet, and this guide, is under heavy development.

Inlet, named because it provides a place of entry, is a C++ library that
provides an easy way to read, store, and access input files for computer
simulations in a variety of input languages.


API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/inlettop.html>`_


Introduction
------------

Inlet provides an easy and extensible way to handle input files for a simulation code.
We provide readers for JSON, Lua, and YAML. Additional
languages can be supported via an implementation of Inlet's Reader interface.

All information read from the input file is stored via Inlet into a user-provided Sidre
DataStore.  This allows us to utilize the functionality of Sidre, such as simulation
restart capabilities.

Inlet is used to define the schema of the information expected in your input file.
That data is then read via a Reader class into the Sidre Datastore.  You can then verify
that the input file met your criteria and use that information later in your code.


Requirements
------------

* Sidre - Inlet stores all data from the input file in the Sidre DataStore
* (Optional) Lua - Inlet provides a Lua reader class that assists in Lua input files


.. toctree::
   :caption: Contents
   :maxdepth: 2

   quick_start
   readers
   simple_types
   advanced_types

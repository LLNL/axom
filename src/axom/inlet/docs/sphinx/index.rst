
Inlet User Documentation
=============================

.. note:: Inlet, and this documentation, is under heavy development.

Inlet, named because it provides a place of entry, is a C++ library that
provides an easy way to read, store, and access information from an input deck.


Introduction
------------

Inlet provides an easy and extensible way to handle input decks for simulation code.
We provide Lua functionality but any language can be used via an inherited Reader class.
Inlet is used to define the structure of the information expected in your input deck.
That data is then read via a Reader class into the Sidre Datastore.  You can then verify
that the input deck met your criteria and use that information later in your code.


Requirements
------------

* Sidre - Inlet stores all data from the input deck in the Sidre DataStore
* (Optional) Lua - Inlet provides a Lua reader class that assists in Lua input decks

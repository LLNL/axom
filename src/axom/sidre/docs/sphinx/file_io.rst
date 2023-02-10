.. ## Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Serial File I/O
******************************************************

Sidre provides file I/O to HDF5 files and file handles and to JSON files.
Serial I/O is accomplished using methods of the ``Group`` class; parallel I/O 
is done using the ``IOManager`` class.

HDF5 is an optional dependency for Sidre. Sidre APIs that rely on HDF5 will be
available only if Sidre was compiled with HDF5 support. The CMake variable
that controls dependency on HDF5 is ``AXOM_USE_HDF5``, defined in the
file axom/config.hpp.

.. _sidre-serial-io:

File I/O using the Group class
-------------------------------

The ``Group`` class provides methods: ``save()``, ``load()``, and
``loadExternalData()`` methods. Each method can be called with a file
name or an HDF5 handle.  The ``save()`` and ``load()`` methods allow users 
to specify a protocol, which indicates how the operation should be performed 
and what file format to use. The ``loadExternalData()`` method takes no 
protocol argument since it is only used with the ``sidre_hdf5`` protocol. 
The ``save()`` method can also take a pointer to an ``Attribute`` object. 
If the attribute pointer is null, all views are saved. If a non-null attribute 
pointer is given, only views with that attribute set will be saved.

The ``load()`` method retrieves the hierarchical data structure stored in the 
file and creates groups, views, and attributes to represent the hierarchy
with its root in the group on which the method is called. By default, the 
contents of the group are destroyed prior to reading the file contents. This 
can be suppressed by passing a Boolean value of true as the 
``preserve_contents`` argument to the ``load()`` method, resulting in the
current group subtree being merged with the file contents.

Usage of the ``save()`` and ``load()`` methods is shown in the following example.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _serial_io_save_start
   :end-before: _serial_io_save_end
   :language: C++

The ``loadExternalData()`` method is used to read "external" data from an HDF5 
file created with the ``sidre_hdf5`` protocol.  This is data referred to by 
*external* views. Such views refer to data that is not stored in Sidre buffers 
but in an external data allocation.

Overloads of the ``save()`` and ``load()`` methods that take HDF5 handles and 
the ``loadExternalData()`` method are used to implement parallel I/O through 
the Sidre ``IOManager`` class.

Please see `Sidre API Documentation <../../../../doxygen/html/sidretop.html>`_ 
for more detailed description of these ``Group`` class methods and the Sidre 
I/O protocols that can be used with them.

Protocols for file I/O
----------------------

``sidre_hdf5``

    This is the protocol that supports the most features. It preserves the
    entire Sidre hierarchy of Groups and Views, with Buffers saved separately
    from the Views that use them. Each Buffer is saved exactly once, to avoid
    duplicating data when a Buffer is used by more than one View. The schemas
    describing the datatypes in each View are saved, and all numerical data is
    saved with full binary accuracy.
     
    ``sidre_hdf5`` is the only protocal that supports the
    ``loadExternalData()`` method and is the only protocoal that allows
    data to be added to root files produced by ``IOManager`` after they are
    first created.

``sidre_conduit_json``

    This protocol saves data in the same layout as ``sidre_hdf5``, but in plain
    text in the JSON format. The datatypes in all Views are saved, but
    floating-point values saved in plain text may not load with exact binary
    equality.

``sidre_json``

    This uses the same hierarchical format as the above protocols and saves
    in the JSON format, but it does not save datatype information for
    scalar and string Views. When loading scalars or strings, the library will
    determine implicitly what type to use, which may not be the same type
    used in the original objects that were saved. For example, the value of
    a 32-bit integer will be saved in plain text with no information about
    its original bit width. When loading, the library may read that value
    into a 64-bit integer.

``conduit_hdf5``

    This saves group as a conduit node hierarchy. Datatypes are preserved, but
    no information about Buffers or external arrays is saved. Data arrays for
    each View are saved; if multiple Views use data from the same Buffer, the
    values from that Buffer will duplicated.

``conduit_json``

    This is the same as ``conduit_hdf5`` except it uses JSON text. The same
    caveats about floating-point accuracy apply here as in the other JSON
    protocols.

``conduit_bin``

    This protocol writes data in conduit's native binary format, and also
    produces an anccompanying JSON file that contains the schema describing
    how to parse the binary data. The hierarchal layout of the data is the
    same as ``conduit_hdf5`` and ``conduit_json``. The binary accuracy of
    floating-point values is preserved. 

``json``

    This protocol saves the smallest amount of information, as it Writes all
    data in a basic JSON hierarchy with no type information beyond
    what can be interpreted implicitly. When loading, the library will allocate 
    the largest supported bit width for integer or floating-point values.

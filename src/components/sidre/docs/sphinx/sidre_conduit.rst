******************************************************
Sidre Interaction with Conduit
******************************************************

.. Again, should the link point to github or to LC conduit home?

Internally, Sidre uses the in-memory data description capabilities of 
`Conduit <https://github.com/LLNL/conduit>`_.  Sidre also leverages Conduit
to facilitate data exchange, demonstrated here as applied to visualization.
The following discussion gives a basic overview of Sidre's capabilities
when combined with Conduit.  Please see the reference documentation for 
more details.

Mesh Blueprint
--------------

The `Mesh Blueprint <http://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html>`_
is a data exchange protocol supported by Conduit, consisting of a
properly-structured Datastore saved as an HDF5 or JSON file and a Conduit index
file.  The Blueprint can accomodate structured or unstructured meshes, with
node- or element-centered fields.  The following example shows how to create a
Blueprint-conforming Datastore containing two unstructured adjacent hexahedrons,
one element-centered field, and one node-centered field.

A simulation will organize its Sidre data as the code design dictates.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _tiny_create_start
   :end-before: _tiny_create_end
   :language: C++

To use the Mesh Blueprint, make a new Datastore conforming to the protocol.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_restructure_start
   :end-before: _blueprint_restructure_end
   :language: C++

Once the Datastore is properly structured, save it, then use Conduit to save the 
index file (ending with `.root`).  This toy data set is small enough that we can
choose to save it as JSON.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_save_start
   :end-before: _blueprint_save_end
   :language: C++

Data Visualization
------------------

The `VisIt <http://visitusers.org/>`_ tool can read in a Blueprint, interpret
the index file, and sensibly display the data contained in the data file.
VisIt produced the following image from the Mesh Blueprint file saved above.

.. image:: visit_conduit.png

.. Cyrus, can you fill in a paragraph here with a really simple recipe?
   The following statement is lame but true and summarizes my personal knowledge.

Configuring and using VisIt is beyond the scope of this document.

Conduit is also a foundational building block for the 
`Ascent <http://ascent.readthedocs.io/>`_ project, which provides a powerful
data analytics and visualization facility (without copying memory) to 
distributed-memory simulation codes.

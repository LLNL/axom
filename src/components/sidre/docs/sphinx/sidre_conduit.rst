.. _sidre-conduit:
******************************************************
Sidre Interaction with Conduit
******************************************************

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
Blueprint-conforming Datastore containing two unstructured adjacent hexahedrons
with one node-centered field and one element-centered field.  In the diagram,
nodes are labeled in black, the node-centered field values are in blue, and
the element-centered field values are in green.

.. image:: figs/tiny_mesh.png
   :width: 300px

A simulation organizes its Sidre data as the code design dictates.
Here is a simple example.

.. image:: figs/ds.png
   :width: 300px

Here is the code to create that Dataset ``ds``.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _tiny_create_start
   :end-before: _tiny_create_end
   :language: C++

To use the Mesh Blueprint, make a new Datastore ``cds`` conforming to the protocol
that points at the data of the original Datastore.  The structure of the
conforming Datastore is shown below (summarizing the
`Mesh Blueprint <http://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html>`_
documentation).

First build top-level groups required by the Blueprint.

.. image:: figs/cds.png
   :width: 650px

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_restructure_toplevel_start
   :end-before: _blueprint_restructure_toplevel_end
   :language: C++

Add the node coordinates.  The Groups in ``cds`` will use the Buffers within ``ds``
as external pointers that have a description: while Sidre and Conduit will use the
array type and shape information, ``cds`` will not deallocate the memory.

.. image:: figs/cdscoords.png
   :width: 650px

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_restructure_coords_start
   :end-before: _blueprint_restructure_coords_end
   :language: C++

Arrange the nodes into elements.  Each simulation has its own knowledge of
topology.  This tiny example doesn't encode topology in ``ds``,
so we must explicitly specify it.

.. image:: figs/cdstopo.png
   :width: 650px

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_restructure_topo_start
   :end-before: _blueprint_restructure_topo_end
   :language: C++

Link the fields into ``cds``.  As with the node positions, ``cds`` gets an external
pointer to the field data.

.. image:: figs/cdsfields.png
   :width: 650px

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_restructure_field_start
   :end-before: _blueprint_restructure_field_end
   :language: C++

Conduit includes a ``verify`` method to test if the structure
of the ``cds`` conforms to the Mesh Blueprint.  This is valuable for writing and 
debugging data adapters.
Once the Datastore is properly structured, save it, then use Conduit to save the 
index file (ending with `.root`).  This toy data set is small enough that we can
choose to save it as JSON.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _blueprint_restructure_save_start
   :end-before: _blueprint_restructure_save_end
   :language: C++

The code listed above produces the files `tinymesh.json` and `tinymesh.root`.
Any code that uses Mesh Blueprint can open and use this pair of files.

Data Visualization
------------------

The `VisIt <http://visitusers.org/>`_ tool can read in a Blueprint, interpret
the index file, and sensibly display the data contained in the data file.
Starting from version 2.13.1, VisIt can open a `.root` file just like any other
data file.  VisIt produced the following image from the Mesh Blueprint file saved above.

.. image:: figs/tiny_mesh_rendered.png
   :width: 600px

Conduit is also a foundational building block for the 
`Ascent <http://ascent.readthedocs.io/>`_ project, which provides a powerful
data analytics and visualization facility (without copying memory) to 
distributed-memory simulation codes.

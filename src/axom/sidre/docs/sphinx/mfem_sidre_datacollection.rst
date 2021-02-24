.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _mfem-sidre-datacollection:

******************************************************
Using Sidre with MFEM
******************************************************

The ``MFEMSidreDataCollection`` class implements `MFEM <https://mfem.org>`_'s 
``DataCollection`` interface for collecting data for a given simulation.
It combines fields and the mesh on which they are defined.  

``MFEMSidreDataCollection`` internally uses the Mesh Blueprint, a standardized
 representation, to store its data.

See the :ref:`Conduit page <sidre-conduit>` for more information on the Mesh Blueprint.

Uses
--------------

The data in an instance of the ``MFEMSidreDataCollection`` class can be saved to a file using a variety of formats.  
This provides visualization and restart capabilities, as these files can also be
loaded back into an instance of the class.

Current options for output formats (protocols) include:

   - ``sidre_hdf5``
   - ``sidre_conduit_json``
   - ``sidre_json``
   - ``conduit_hdf5``

.. Note::
   The ``AXOM_USE_HDF5`` build option must be enabled to use the HDF5-based formats.

A typical visualization example might look like the following:

.. code-block::

   mfem::Mesh mesh;
   // ...read the mesh in from a file
   MFEMSidreDataCollection dc("data_name", &mesh);

   mfem::GridFunction soln;
   // ...initialize with FiniteElementSpace, etc

   // Note: Any number of fields can be registered
   dc.RegisterField("solution", &soln);

   // Save the initial state
   dc.SetCycle(0); // Iteration counter
   dc.SetTime(0.0); // Simulation time
   // Filename and protocol, both of which are optional
   dc.Save("data_name", "sidre_hdf5");

   for (int i = 0; i < n_iter; i++)
   {
      // Calculate the next iteration of the solution field, then...
      
      dc.SetCycle(i);
      dc.SetTime(dt * i);
      dc.Save("data_name", "sidre_hdf5");
   }

See the ``sidre_mfem_datacollection_vis`` example for a more thorough example of the above functionality.

.. Note::
   The ``owns_mesh_data`` option must be set to true when constructing an instance of the class for the 
   mesh to be read back in properly when a restart occurs.

.. Warning::
   Although the ``mfem::DataCollection`` interface provides functionality for collection quadrature fields,
   this is not supported by ``MFEMSidreDataCollection``.

Experimental support for complete reconstruction of a simulation's mesh and fields is also provided by
``MFEMSidreDataCollection``.  That is, when an output file is read in using ``MFEMSidreDataCollection::Load``,
the data read in will be used to reconstruct MFEM objects than can be accessed with the ``GetField`` and
``GetMesh`` methods.  

This functionality is demonstrated in the ``sidre_mfem_datacollection_restart`` example, which is a stripped-down
example of how a simulation code might utilize the automatic reconstruction logic when loading in a datastore.

.. Note::
  The mesh/field reconstruction logic requires that the save file was created with the ``MFEMSidreDataCollection``
  class. In `Mesh Blueprint <http://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html>`_ terms, the
  following constraints are imposed on the structure of the data:

  * There must be a coordinate set named ``coords``
  * There must be a topology named ``mesh`` with corresponding attributes stored in a field named ``mesh_material_attribute``
  * There must be a topology named ``boundary`` with corresponding attributes stored in a field named ``boundary_material_attribute``

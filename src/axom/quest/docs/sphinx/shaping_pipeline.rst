.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. pipeline:

Shaping Overview
================

Shaping is the process of overlaying additional detail onto a mesh by converting
shape geometry into materials described as volume fractions within each mesh zone.
Axom;s Klee component describes the models used for shaping. Klee shapes include
a material name, a file path that contains the shape geometry, replacement rules,
and transforms that can be applied to the shape geometry. Axom's Quest component
contains shaping infrastructure that takes the shapes from Klee and generates the
volume fractions for the shapes on a target mesh.

Shaping Pipeline
=================

Quest shaping classes provide an interface that corresponds to several stages in
a shaping pipeline. At the high level, a collection of shapes is read using Klee
and these shapes are run through the shaping pipeline to produce volume fractions
on the mesh. This process is repeated for each shape until at the end of shaping,
the target mesh contains a set of volume fraction fields that encode the shaped-in
materials.

We include relevant Axom headers

.. code-block:: yaml

  #include <axom/klee.hpp>
  #include <axom/quest.hpp>
  #include <axom/sidre.hpp>
  #include <axom/slic.hpp>
  #include <memory>


Axom shaping classes operate on MFEM meshes and use grid functions to represent the
volume fraction fields that encode each shaped-in material. MFEM meshes can be loaded
or constructed in memory.

 * More information on MFEM is covered at the [MFEM Examples page](https://mfem.org/features/#extensive-examples).

The MFEM mesh also needs an associated data collection, *(shapingDC)*, to contain the
grid functions. Axom provides *MFEMSidreDataCollection*, a derived class of MFEM's
DataCollection class that can interoperate with Axom's Sidre component.

.. literalinclude:: ../../examples/shaping_driver.cpp
   :start-after: _load_mesh_start
   :end-before: _load_mesh_end
   :language: C++


We create the desired shaper _(SamplingShaper shown here)_ and set its parameters.
Note that some settings are common to both sampling and intersection shaping while
some settings are available only in derived classes.

.. code-block:: yaml

  using quest = axom::quest;
  auto shaper = std::make_unique<quest::SamplingShaper>();
  shaper.setSamplesPerKnotSpan(25);
  shaper.setVerbosity(true);


The shaping driver can optionally be initialized with volume fractions from the
calling code. This step can be skipped if the volume fractions are to be produced
solely using shapes from Klee. Axom shaping classes operate on MFEM meshes with
grid functions that encode the volume fractions. The names of the volume fraction
grid functions all have a *"vol_frac_"* prefix, followed by the name of the material.

.. literalinclude:: ../../examples/shaping_driver.cpp
   :start-after: import_volume_fractions_start
   :end-before: import_volume_fractions_end
   :language: C++


After all shaping input data have been read, the actual *shaping pipeline* can begin. This is
where each shape is processed within the shaper. The *loadShape()* method is called to
make the shape be loaded from its geometry file. The *prepareShapeQuery()* method builds
internal data structures that aid in shaping. The *runShapeQuery()* method builds the
query mesh and intersects it with the target mesh, creating the volume fractinons for
the shape. The *applyReplacementRules()* method incorporates the shape's volume fractions
into the existing volume fraction grid functions, subject to the replacement rules defined
for the shape. Finally, the finalizeShapeQuery()* method performs internal cleanup in the
shaper so it is ready to process the next shape.

.. literalinclude:: ../../examples/shaping_driver.cpp
   :start-after: _shaping_pipeline_begin
   :end-before: _shaping_pipeline_end
   :language: C++

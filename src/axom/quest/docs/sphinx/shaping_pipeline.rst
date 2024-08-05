.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. pipeline:

Shaping Overview
================

Shaping is the process of overlaying additional detail into a mesh by converting
shape geometry into materials described as volume fractions within each mesh zone.
Shaping is used when it is not feasible or practical to directly build features
into the mesh.

.. figure:: figs/shaping_overview.png
   :width: 600px

   Shaping permits details to be added into meshes.


Axom\'s Klee component describes the models used for shaping. Klee shapes include
a material name, a file path that contains the shape geometry, replacement rules,
and transforms that can be applied to the shape geometry. Axom\'s Quest component
contains shaping infrastructure that takes the shapes from Klee and generates the
volume fractions for the shapes on a target mesh.

Shaping Pipeline
=================

Shaping involves creating a target mesh and data collection, reading a shape set,
creating a shaper, and then iterating over the shapes to pass each shape into the
shaper. The shaper is responsible for determining overlap between the shape and
the target mesh and producing grid functions that contain that overlap, or volume
fraction.

First, we include relevant Axom headers:

.. code-block:: c++

  #include <axom/klee.hpp>
  #include <axom/quest.hpp>
  #include <axom/sidre.hpp>
  #include <axom/slic.hpp>

  using quest = axom::quest;
  using klee = axom::klee;
  using slic = axom::slic;
  using sidre = axom::sidre;

Axom shaping classes operate on MFEM meshes and use grid functions to represent the
volume fraction fields that encode each shaped-in material. MFEM meshes can be loaded
or constructed in memory.

 * More information on MFEM is covered at the `MFEM Examples page <https://mfem.org/features/#extensive-examples>`_.

The MFEM mesh also needs an associated data collection, *(shapingDC)*, to contain the
grid functions. Axom provides *MFEMSidreDataCollection*, a derived class of MFEM's
DataCollection class that can interoperate with Axom\'s Sidre component.

.. literalinclude:: ../../examples/shaping_driver.cpp
   :start-after: _load_mesh_start
   :end-before: _load_mesh_end
   :language: C++


We create the desired shaper *(SamplingShaper shown here)* and set its parameters.
Note that some parameters are common to each shaper type while others are specific to
sampling or intersection shaping classes.

.. code-block:: yaml

  auto shaper = new quest::SamplingShaper;
  shaper->setSamplesPerKnotSpan(25);
  shaper->setVerbosity(true);


The shaper will operate on shapes, which can be obtained by reading a Klee shape set.

.. code-block:: c++

  auto shapeSet = klee::readShapeSet("/path/to/klee/file");


The shaper can **optionally** be pre-initialized with volume fractions from the
calling code. This step can be skipped if the volume fractions are to be produced
solely using shapes from Klee. Volume fraction field grid functions have a *"vol_frac_"* prefix,
followed by the name of the material. Grid functions are registered with the
*shapingDC* data collection and the shaper's *importInitialVolumeFractions()*
method aids in importing the initial volume fractions onto the target mesh.

.. literalinclude:: ../../examples/shaping_driver.cpp
   :start-after: import_volume_fractions_start
   :end-before: import_volume_fractions_end
   :language: C++


After all shaping input data have been read, the actual *shaping pipeline* can begin. This is
where each shape is processed within the shaper and this procedure applies to both sampling
and intersection shaping. The *loadShape()* method is called to
make the shape be loaded from its geometry file. The *prepareShapeQuery()* method builds
internal data structures that aid in shaping. The *runShapeQuery()* method builds the
query mesh and intersects it with the target mesh, creating the volume fractinons for
the shape. The *applyReplacementRules()* method incorporates the shape's volume fractions
into the existing volume fraction grid functions, subject to the replacement rules defined
for the shape. Finally, the *finalizeShapeQuery()* method performs internal cleanup in the
shaper so it is ready to process the next shape.

.. literalinclude:: ../../examples/shaping_driver.cpp
   :start-after: _shaping_pipeline_begin
   :end-before: _shaping_pipeline_end
   :language: C++

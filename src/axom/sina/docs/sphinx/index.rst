.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

===================
Sina C++ User Guide
===================

The Sina C++ library can read and write JSON files in the Sina schema. It
can be used by simulation applications to summarize run data to be ingested
into a database using the Sina tool suite.

The top-level object in the Sina schema is the Document. It contains lists
of Record and Relationship objects. The example below shows the basics.
For more details, see the `Tutorial <tutorial.rst>`_.

.. literalinclude:: ../../examples/sina_basic.cpp
   :language: cpp

After running the above, the file "MySinaData.json" will contain the
following:

.. code:: json

    {
        "records": [
            {
                "application": "My Sim Code",
                "local_id": "run1",
                "type": "run",
                "user": "jdoe",
                "version": "1.2.3"
            }
        ],
        "relationships": []
    }

.. toctree::
   :caption: Contents
   :maxdepth: 2

   tutorial
   core_concepts
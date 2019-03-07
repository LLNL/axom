.. ##
.. ## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory
.. ##
.. ## LLNL-CODE-741217
.. ##
.. ## All rights reserved.
.. ##
.. ## This file is part of Axom.
.. ##
.. ## For details about use and distribution, please read axom/LICENSE.
.. ##

**********************
Implementation details
**********************

.. _policy-label:

Policy-based design
-------------------

Handling the combinatorial explosion of features; avoid paying for what we don't need

* SizePolicy, StridePolicy, OffsetPolicy (compile time vs. runtime)
* IndirectionPolicy (none, C-array, std::vector, custom, e.g. mfem::Array)
* SubsettingPolicy (none, virtual parent, concrete parent)
* OwnershipPolicy (local, sidre, other repository)


Feature diagram of OrderedSet policies (subset).

.. figure:: figs/orderedset_feature_diagram.png
   :figwidth: 100%
   :alt: Feature diagram for slam's ordered set

The figure shows how certain these policies interact with the subscript operator.


.. _setup-label:

Simplifying mesh setup
----------------------

* Builder classes
    * Chained initialization using named-parameter idiom
* Generator classes to simplify types

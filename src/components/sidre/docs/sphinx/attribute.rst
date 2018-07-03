.. ##
.. ## Copyright (c) 2017-18, Lawrence Livermore National Security, LLC.
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

.. _attribute-label:

==========
Attribute
==========

Attributes provide storage for View metadata.  When each
Attribute is constructed in the DataStore, it gets a name and a default value.
Each View inherits all of its DataStore's Attributes at their default values.
A program may explicitly set any of these Attributes for each
View.  The program may also query the value of a View's Attribute, query whether
the Attribute was explicitly set, or clear the Attribute back to its default
value.

Attribute values are available for a program to use in its own logic.  If a
program provides an Attribute pointer to :code:`Group::save()` (discussed in the next
section), only Views with that Attribute explicitly set will be saved.  Further
extensions to Sidre that use Attributes and their values are planned.

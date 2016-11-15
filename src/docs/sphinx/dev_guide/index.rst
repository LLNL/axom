.. ##
.. ## Copyright (c) 2016, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory.
.. ##
.. ## All rights reserved.
.. ##
.. ## This file cannot be distributed without permission and
.. ## further review from Lawrence Livermore National Laboratory.
.. ##

CS Toolkit Developer Guide
=========================================================

This guide describes key aspects of the software development process for
the CS Toolkit project. The guide is intended for all team members and 
contributors. It is especially helpful for familiarizing new individuals
about how the team works. The overarching goal is that all CS Toolkit 
development follows similar practices to ensure robustness, consistency,
ease of use, and comprehensive testing. Everyone who contributes to the 
CS Toolkit should be aware of and follow these guidelines. However, the
benefits of uniformity should be balanced with allowances for individual 
preferences, which may be superior to rigid conventions in certain situations.

The guidelines should not be viewed as fixed for all time. They should evolve 
with project needs and be improved when processes can be improved. Changes 
should be agreed to by team members after assessing their merits using 
their collective professional judgment. When changes are made, this guide
should be updated accordingly.

This guide is generated using Sphinx and is written in the *reStructuredText* 
markup language. The document source lives in the CS Toolkit source code 
repository. You can build it from source by using the proper make system 
target::

   $ make dev_guide_docs

If you need to edit the guide and need information about reStructuredText, 
please see the `Sphinx Guide <http://www.sphinx-doc.org/en/stable/index.html>`_. 


**Contents:**

.. toctree::
   :maxdepth: 3

   dev_model
   dev_tools
   add_new_component
   testing
   misc_tasks



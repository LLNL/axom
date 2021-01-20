.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _continuous_integration-label:

*******************************
Continuous Integration 
*******************************

The Axom project uses two CI tools,
`Azure Pipelines <https://azure.microsoft.com/en-us/services/devops/pipelines/>`_
via Github and `Bamboo <https://www.atlassian.com/software/bamboo>`_ 
on the LC Restricted Zone (RZ).

.. _azure_pipelines-label:

===============
Azure Pipelines 
===============

Every Pull Request created on Github is automatically run through a series of
CI jobs to ensure that the Axom source builds and passes our unit tests.
These configurations mimic the LC compilers and systems as closely as possible
via Docker containers that have our third-party libraries pre-built on them.


.. _bamboo-label:

==========
RZ Bamboo 
==========

We use the `Axom RZ Bamboo project <https://rzlc.llnl.gov/bamboo/browse/ASC>`_ 
primarily for testing the develop branch against the various LC System Types.
There are two types of Bamboo Plans which are automatically triggered to build
and run tests against the develop branch.  The first is triggered nightly and
on any repository change on Github that builds and tests the current Axom source
against previously built third-party libraries.  The second is triggered nightly
and builds and tests the complete third-party library build and then the Axom source
against those libraries.

This plan may be run manually at any time by selecting the plan and clicking
on 'Run plan' as described above. Each member of the team receives an email 
notification every morning about the current state of all jobs.


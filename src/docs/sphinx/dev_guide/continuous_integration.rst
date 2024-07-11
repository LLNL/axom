.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _continuous_integration-label:

*******************************
Continuous Integration 
*******************************

The Axom project uses two CI tools,
`Azure Pipelines <https://azure.microsoft.com/en-us/services/devops/pipelines/>`_
via GitHub and `GitLab CI <https://docs.gitlab.com/ee/ci/>`_ 
on the LLNL LC Collaboration Zone (CZ).

.. _azure_pipelines-label:

===============
Azure Pipelines 
===============

Every Pull Request created on GitHub is automatically run through a series of
CI jobs to ensure that the Axom source builds and passes our unit tests.
These configurations mimic the LC compilers and systems as closely as possible
via Docker containers that have our third-party libraries pre-built on them.

Axom's GitHub project is also configured to require pull requests to pass checks 
from our LC GitLab CI (as described below).


.. _gitlab-label:

==========
LC GitLab 
==========

We also maintain a mirror of the `Axom project on LLNL's LC GitLab instance <https://lc.llnl.gov/gitlab/axom/axom>`_
primarily for testing Axom pull requests against the various LC System Types and compilers.

There are two types of GitLab plans.
The first is triggered automatically by pull requests on GitHub,
while the second runs nightly and tests
Axom's ``develop`` branch against a new build of our third-party library stack.

Our GitLab CI configuration also allows manual runs. To initiate a new run, 
navigate to the `CI/CD page, <https://lc.llnl.gov/gitlab/axom/axom/-/pipelines>`_
click on the "Run pipeline" button and select the branch to test.


.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

****************************************
Axom Development Process Summary
****************************************

This section provides a high-level overview of key Axom software development
topics and includes links to more detailed discussion.


======================================================
Software Development Cycles
======================================================

The Axom team uses a sprint-based development process. We collect
and track issues (bugs, feature requests, tasks, etc.) using ``Github``
and define a set of development tasks (i.e., issues) to complete for each 
sprint. While the team meets to discuss issues and plan which ones will be 
worked in each sprint, developers of individual Axom components may plan and 
schedule work in any way that works for them as long as this is coordinated
with other team efforts. Work performed in each sprint work period is tracked 
as a single unified sprint encompassing activities for the entire project.



======================================================
Software Releases and Version Numbers
======================================================

Typically, Axom releases are done when it makes sense to make new features
or other changes available to users. A release may coincide with the completion
of a sprint cycle or it may not.

See :ref:`release-label` for a description of the Axom release process.

The Axom team follows the **semantic versioning** scheme for assigning
release numbers. Semantic versioning conveys specific meaning about 
the code and modifications from version to version by the way version
numbers are constructed.

See :ref:`semver-label` for a description of semantic versioning.


======================================================
Branch Development
======================================================

The Axom project has a ``Github`` project space and the team follows 
the **Gitflow** branching model for software development and reviews. Gitflow 
is a common workflow centered around software releases. It makes clear which 
branches correspond to which phases of development and those phases are 
represented explicitly in the structure of the source code repository. As 
in other branching models, developers develop code locally and push their 
work to a central repository.

See :ref:`gitflow-label` for a detailed description of how we use Gitflow.


======================================================
Code Reviews and Acceptance
======================================================

Before any code is merged into one of our main Gitflow branches (i.e., develop 
or main), it must be adequately tested, documented, and reviewed 
for acceptance by other team members. The review process is initiated via 
a *pull request* on the Axom Github project.

See :ref:`pullrequest-label` for a description of our review process and 
how we use pull requests.


======================================================
Testing and Code Health
======================================================

Comprehensive software testing processes and use of code health tools (e.g., 
static analysis, memory checkers, code coverage) are essential ingredients 
in the Axom development process.

See :ref:`testing-label` for a description of our software testing process,
including *continuous integration*.


======================================================
Software Development Tools
======================================================

In addition to the tools listed above, we use a variety of other tools to help
manage and automate our software development processes. The *tool philosophy*
adopted by the Axom project focuses on three central tenets:

  * Employ robust, commonly-used tools and don't re-invent something that already exists.
  * Apply tools in ways that non-experts find easy to use.
  * Strive for automation and reproducibility.

The main interaction hub for Axom developers is the **Atlassian
tool suite** on the Livermore Computing Collaboration Zone (CZ). These tools
can be accessed through the `MyLC Portal <https://lc.llnl.gov>`_.
Developer-level access to Axom project spaces in these tools requires 
membership in the LC group 'axomdev'. If you are not in this group, and need 
to be, please send an email request to 'axom-dev@llnl.gov'.

The main tools we use are listed below. Please navigate the links
provided for details about how we use them and helpful information about 
getting started with them.

* **Confluence.**  We use the `Axom Confluence space <https://lc.llnl.gov/confluence/display/ASCT>`_ for team discussion (e.g., hashing out design ideas), maintaining meeting notes, etc.

* **Github.** We use the `Axom Github project <https://github.com/LLNL/axom>`_ to manage our issues and Git repository which contains the Axom source code, build configurations, scripts, test suites, documentation, etc.

  * See :ref:`github-label` for more information about how we use Git and Github.

* **Bamboo.** We use Bamboo for continuous integration to ensure code quality on our LC systems.:  `Axom RZ Bamboo project <https://rzlc.llnl.gov/bamboo/browse/ASC>`_

  * See :ref:`bamboo-label` for more information about how we use Bamboo.

* **Azure Pipelines.** We use Azure Pipelines for continuous integration to ensure every code change passes a
  level of quality before being merged.:  `Azure Pipelines <https://azure.microsoft.com/en-us/services/devops/pipelines/>`_

  * See :ref:`azure_pipelines-label` for more information about how we use Azire Pipelines.

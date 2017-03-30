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

****************************************
Axom Development Process Summary
****************************************

This section provides a brief summary of various aspects of Axom software 
development. 

The main interaction hub for Axom developers is the Atlassian 
tool suite on the Livermore Computing Collaboration Zone (CZ). These tools 
can be accessed through the `MyLC Portal <https://lc.llnl.gov>`_. 
Developer-level access to Axom project spaces requires membership in the LC 
group 'axomdev'. If you are not in this group and need to be, please send 
an email request to 'axom-dev@llnl.gov'.


======================================================
Development and Release Cycles
======================================================

The Axom team uses an agile, sprint-based development process. 
We target a set of development tasks for each quarterly (3 month) releease 
cycles for the project as a whole. Developers of individual software 
components in Axom may plan work on a finer time-scale if they choose to do so.
However, project-wide schedule planning for a release is done every three 
months. Work performed in each release cycle is tracked as a single unified
sprint encompassing the entire project.

See :ref:`releasecycle-label` for more information about how we do release 
planning and tracking. 


======================================================
Release Versioning
======================================================

The Axom team follows the *semantic* versioning scheme for assigning
release numbers. Semantic versioning is a methodology for assigning version 
numbers to software releases in a way that conveys specific meaning about 
the code and modifications from version to version. 

See :ref:`semanticversioning-label` for details on how we apply semantic 
versioning.


======================================================
Branch Development
======================================================

For software development and reviews, the Axom team follows the 'Gitflow' 
branch development model. Gitflow is a branching model centered around 
software releases. It is a simple workflow that makes clear which branches 
correspond to which phases of development and those phases are represented 
explicitly in the structure of the source code repository. As in other 
branching models, developers develop code locally and push their work to 
a central repository.

See :ref:`gitflow-label` for a detailed description of how we use Gitflow.


======================================================
Code Reviews and Acceptance
======================================================

Before any code is merged into one of our main Gitflow branches (i.e., develop 
or master), it must be adequately tested and documented. Then it is reviewed 
for acceptance by other team members. The review process is initiated via 
a *pull request* on the Axom Bitbucket project.

See :ref:`pullrequest-label` for a description of our review process and use of
pull requests.


======================================================
Testing and Code Health
======================================================

Sound software testing processes and use of code health tools (e.g., static
analysis, memory checkers, code coverage) are essential ingredient in the
Axom development process.

See :ref:`testing-label` for a description of our software testing process,
including *continuous integration*.


======================================================
Software Development Tools
======================================================

We use a variety of tools for software development. Our tool philosophy has
three main tenets:

  * Employ robust, commonly-used tools. Don't re-invent something that already exists
  * Apply tools in ways that are easy for non-experts
  * Strive for automation and reproducibility

The main tools we use are listed below. Details about how we use
them and helpful information about getting started can be found via the 
provided links.

* We use the `Axom Confluence space <https://lc.llnl.gov/confluence/display/ASCT>`_ for team discussion, planning, maintaining meeting notes, etc.
* We use the `Axom Bitbucket project <https://lc.llnl.gov/bitbucket/projects/ATK>`_ to manage our Git repository which contains the Axom source code, build configurations, scripts, test suites, documentation, etc.
* We use the `Axom JIRA project <https://lc.llnl.gov/jira/projects/ATK>`_ for issue tracking.
* We use the `Axom Bamboo project <https://lc.llnl.gov/bamboo/browse/ASC>`_ for 
continuous integration and automated testing.


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


.. _pullrequest-label:


======================================================
Pull Requests and Code Reviews
======================================================

Before any code is merged into the develop or master branches, it
must be tested, reviewed, and accepted. Submitting a pull request on
the Axom Bitbucket project to merge a branch into develop or master 
initiates the test and review processes. All builds and tests must pass 
for a pull request to be approved. Also, it is expected that unit tests 
are constructed to exercise any new functionality that is introduced. This 
will be assessed by reviewers of each pull request. See :ref:`testing-label` 
for more information about testing.

Code changes in a pull request must be accepted by at least one member
of the Axom development team other than the originator of the pull
request. It is recommended to have several team members review pull 
requests, especially when changes affect APIs. Pull request reviewers can be 
selected on Bitbucket when the pull request is created. Changes reviewed by 
the team are accepted, rejected, or commented on for improvement; e.g., 
issues to be addressed, suggested changes, etc. Pull requests can be updated
with additional changes as needed. When a pull request is approved, it can 
be merged. If the merged branch is no longer needed for development, it 
should be deleted.

In addition to successful compilation and passing tests, changes to the 
develop and master branches should be scrutinized in other ways and using 
other tools. In particular :

* The code should compile cleanly at the highest warning level with the 
  main compilers supported by the project. All warnings **must** be 
  understood and eliminated if possible. Reducing a compiler warning 
  level to eliminate warning messages **is not** acceptable.

  Compiler warnings, while seemingly innocuous at times, often indicate
  problems that do not appear until later or until specific run-time
  conditions are encountered.

* Static analysis tools **should** be applied to the code using tools such
  as `cppcheck`, etc. to identify potential implementation issues.

* Runtime memory checking, using a  tool such as Valgrind, **should** be 
  performed to verify that there are no leaks or other memory issues. 


---------------------
Pull Request Summary
---------------------

  #. When code is ready to be considered for acceptance, create a pull request
     on the Axom Bitbucket project. Identify the appropriate reviewers 
     and add them to the pull request.

  #. Code must build successfully and all relevant tests must pass, including
     new tests required for new functionality.

  #. All issues (build failures, test failures, reviewer requests) must be 
     addressed before a pull request is accepted.

  #. Pull requests must be approved by at least one member of development 
     team other than the pull request originator.

  #. When a pull request is approved it may be merged. If the merged branch is
     no longer needed, it should be deleted. This can be done when merging
     with Bitbucket. 


---------------------------
Code Review Checklist
---------------------------

Beyond build and test correctness, we also want to ensure that code follows
common conventions before acceptance. The following list summarizes concerns 
we want to identify during pull request reviews and resolve before a pull 
request is approved for merging. Please see the Axom coding guidelines
**(insert link here)**
document for details. 

 #. A new file or directory must be located in in the proper location; e.g.,
    in the same directory with existing files supporting related functionality.
 #. File contents must be organized clearly and structure must be consistent 
    with conventions. 
 #. Namespace and other scoping conventions must be followed. 
 #. Names (files, types, methods, variables, etc.) must be clear, easily
    understood by others, and consistent with usage in other parts of the code.
    Terminology must be constrained; i.e., don't introduce a new term for 
    something that already exists and don't use the same term for different 
    concepts.
 #. Documentation must be clear and follow conventions. Minimal, but adequate, 
    documentation is preferred.
 #. Implementations must be correct, robust, portable, and understandable to
    other developers.
 #. Adequate tests (unit and performance) tests must be added for new 
    functionality.



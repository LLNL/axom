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

======================================================
The Code
======================================================

Our Git repository contains the Toolkit source code, documentation, test 
suites and all files and scripts used for configuring and building the code.
The repo lives in our 
`CZ Bitbucket project <https://lc.llnl.gov/bitbucket/projects/ATK>`_.

We use our `JIRA project space <https://lc.llnl.gov/jira/browse/ATK>`_ for 
issue tracking.


--------------------------------
Getting the Code
--------------------------------

Access to our repository and Atlassian tools requires membership in the LC group 
``toolkit``. If you're not in the group, please send email to 
'asctoolkit-dev@llnl.gov' and request to be added.

SSH keys
^^^^^^^^^

If you have not used Bitbucket before, you will need to
`create an SSH key <https://confluence.atlassian.com/bitbucketserver/creating-ssh-keys-776639788.html>`_ 
and `add the key to your Bitbucket profile <https://confluence.atlassian.com/bitbucketserver/ssh-user-keys-for-personal-use-776639793.html>`_.

Cloning the repo
^^^^^^^^^^^^^^^^^^

To clone the repo into your local working space, type the following::

  $ git clone --recursive ssh://git@cz-bitbucket.llnl.gov:7999/atk/asctoolkit.git

Important notes:

  * You don't need to remember the URL for the Toolkit repo above. It can be
    found by going to the CS Toolkit repo on our Bitbucket project and
    clicking on the 'Clone' action button that appears when you hover your
    mouse cursor over the ellipses on the top left of the web page.
  * The ``--recursive`` argument above is needed to pull in our build system,
    called *BLT*. BLT is standalone product that is used by other code projects;
    thus, it lives in its own reposiotory. It is a Git sub-module of the 
    Toolkit.


--------------------
Repository Layout
--------------------

If you need to look through the repository, this section explains how it is
organized. If you do not need this information and just want to build the
code, please continue on to the next section.

The top-level CS Toolkit directory contains three directories:

  scripts
    Scripts that we maintain to simplify development and usage tasks
  host-configs
    Detailed configuration information for platforms and 
    compilers we support (more about these later)
  src
    The bulk of the repo contents.
    Within the **src** directory, you will find the following directories:
    
    blt
      BLT build system submodule is cloned here
    cmake
      Toolkit-specific CMake customizations to BLT build system
    components
      Files for individual Toolkit components
    docs
      General Toolkit documentation files
    thirdparty
      Tests to make sure TPLs are built properly

In the **components** directory, you will find a directory for each of the
Toolkit components. Although there are dependencies among them, each is 
developed and maintained in a largely self-contained fashion. Toolkit 
component dependencies are essentially treated as library dependencies.
Each component directory contains subdirectories for: user documentation,
examples, tests, and the component header and implementation files.

 

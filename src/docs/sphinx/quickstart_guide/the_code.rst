.. ##
.. ## Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory.
.. ##
.. ## LLNL-CODE-741217
.. ##
.. ## All rights reserved.
.. ##
.. ## This file is part of Axom.
.. ##
.. ## For details about use and distribution, please read axom/LICENSE.
.. ##

======================================================
The Code
======================================================

Our Git repository contains the Axom source code, documentation, test 
suites and all files and scripts used for configuring and building the code.
The repository lives in our 
`CZ Bitbucket project <https://lc.llnl.gov/bitbucket/projects/ATK>`_.

We use our `JIRA project space <https://lc.llnl.gov/jira/browse/ATK>`_ for 
issue tracking. Please report issues, feature requests, etc. there or send 
email to the Axom development team.


--------------------------------
Getting the Code
--------------------------------

Access to our repository and Atlassian tools requires membership in the LC 
group ``axom``. If you're not in the group, please send email to 
'axom-dev@llnl.gov' and request to be added.

SSH keys
^^^^^^^^^

If you have not used Bitbucket before, you will need to
`create an SSH key <https://confluence.atlassian.com/bitbucketserver/creating-ssh-keys-776639788.html>`_ 
and `add the key to your Bitbucket profile <https://confluence.atlassian.com/bitbucketserver/ssh-user-keys-for-personal-use-776639793.html>`_.

Cloning the repo
^^^^^^^^^^^^^^^^^^

To clone the repo into your local working space, type the following::

  $ git clone --recursive ssh://git@cz-bitbucket.llnl.gov:7999/atk/axom.git

Important notes:

  * You don't need to remember the URL for the Axom repo above. It can be
    found by going to the Axom repo on our Bitbucket project and
    clicking on the 'Clone' action button that appears when you hover your
    mouse cursor over the ellipses at top of the web page on the left.
  * The ``--recursive`` argument above is needed to pull in our build system,
    called *BLT*, which is a standalone product that lives in its own repository.
    It is a Git sub-module in Axom.


--------------------
Repository Layout
--------------------

If you need to look through the repository, this section explains how it is
organized. If you do not need this information and just want to build the
code, please continue on to the next section.

The top-level Axom directory contains three directories:

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
      Axom-specific CMake customizations to BLT build system
    components
      Files for individual Axom components (see below)
    docs
      General Axom documentation files
    thirdparty
      Tests to make sure TPLs are built properly

In the **components** directory, you will find a directory for each of the
Axom components. Although there are dependencies among them, each is 
developed and maintained in a largely self-contained fashion. Axom 
component dependencies are essentially treated as library dependencies.
Each component directory contains subdirectories for the component header
and implementation files, as well as user documentation, examples and tests.

 

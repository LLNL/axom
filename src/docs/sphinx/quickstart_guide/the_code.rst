.. ##
.. ## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

The Code
========

Our Git repository contains the Axom source code, documentation, test
suites and all files and scripts used for configuring and building the code.
The repository lives in our
`CZ Bitbucket project <https://lc.llnl.gov/bitbucket/projects/ATK>`_.

We use our `JIRA project space <https://lc.llnl.gov/jira/browse/ATK>`_ for
issue tracking. Please report issues, feature requests, etc. there or send
email to the Axom development team.


Getting the Code
----------------

Access to our repository and Atlassian tools requires membership in the LC
group ``axom``. If you're not in the group, please send email to
'axom-dev@llnl.gov' and request to be added.

SSH keys
^^^^^^^^

If you have not used Bitbucket before, you will need to
`create an SSH key <https://confluence.atlassian.com/bitbucketserver/creating-ssh-keys-776639788.html>`_
and `add the key to your Bitbucket profile <https://confluence.atlassian.com/bitbucketserver/ssh-user-keys-for-personal-use-776639793.html>`_.

Cloning the repo
^^^^^^^^^^^^^^^^

To clone the repo into your local working space, type the following::

  $ git clone --recursive ssh://git@cz-bitbucket.llnl.gov:7999/atk/axom.git

Important notes:

  * You don't need to remember the URL for the Axom repo above. It can be
    found by going to the Axom repo on our Bitbucket project and
    clicking on the 'Clone' action button that appears when you hover your
    mouse cursor over the ellipses at top of the web page on the left.
  * The ``--recursive`` argument above is needed to pull in our build system,
    called *BLT*, which is a standalone product that lives in
    `its own repository <https://github.com/llnl/blt>`_.
    It is a Git sub-module in Axom.  Documentation for BLT can be
    found `here <https://llnl-blt.readthedocs.io/en/latest/>`_.
  * If you forget to pass the ``--recursive`` argument to the git clone command,
    the following commands can be typed after cloning:

    .. code:: bash

      $ cd axom
      $ git submodule init
      $ git submodule update


Repository Layout
-----------------

If you need to look through the repository, this section explains how it is
organized. If you do not need this information and just want to build the
code, please continue on to the next section.

The top-level Axom directory contains three directories:

  host-configs
    Detailed configuration information for platforms and compilers we support.

    See :ref:`hostconfig-label` for more information.
  scripts
    Scripts that we maintain to simplify development and usage tasks
  src
    The bulk of the repo contents.

    Within the **src** directory, you will find the following directories:

    axom
      Directories for individual Axom components (see below)
    cmake
      Axom's build system lives here.

      The BLT submodule is cloned into the **blt** subdirectory.
    docs
      General Axom documentation files
    examples
      Example programs that utilize Axom in their build systems
    thirdparty
      Built-in third party libraries with tests to ensure they are built properly.

In the **axom** directory, you will find a directory for each of the
Axom components. Although there are dependencies among them, each is
developed and maintained in a largely self-contained fashion. Axom
component dependencies are essentially treated as library dependencies.
Each component directory contains subdirectories for the component header
and implementation files, as well as user documentation, examples and tests.

Axom has the following built-in third party libraries:

  `fmt <http://fmtlib.net/latest/index.html>`_
      BSD-licensed string formatting library
  `sparsehash <https://github.com/sparsehash/sparsehash>`_
      BSD-licenced associative containers for C++

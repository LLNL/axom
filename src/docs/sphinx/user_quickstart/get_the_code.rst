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
Accessing the Code
======================================================

The CS Toolkit Git repository is housed on our
`CZ Bitbucket project <https://https://lc.llnl.gov/bitbucket/projects/ATK>`_
The repo contains the Toolkit source code, build configurations, scripts,
test suites, documentation, etc. Repository access requires membership in
the LC group 'toolkit. If you're not in the group, send email to
'asctoolkit-dev@llnl.gov' and request to be added.

-------------
SSH keys
-------------

If you have not used Bitbucket before, you will need to
`create an SSH key <https://confluence.atlassian.com/bitbucketserver/creating-ssh-keys-776639788.html>`_ and `add the key to your Bitbucket profile <https://confluence.atlassian.com/bitbucketserver/ssh-user-keys-for-personal-use-776639793.html>`_.

--------------------------------
Getting a local copy of the code
--------------------------------

To use the code, you must clone the repo into your working space. This is done 
by typing::

  $ git clone --recursive ssh://git@cz-bitbucket.llnl.gov:7999/atk/asctoolkit.git

Important notes:

  * You don't need to remember the URL for the Toolkit repo above. It can be
    found by going to the CS Toolkit repo on our Bitbucket project and
    clicking on the 'Clone' Action button that appears when you hover your
    mouse cursor over the ellipses on the top left of the web page.
  * The '--recursive' argument above is needed to pull the BLT build system into
    your local copy of the repo. It is a Git sub-module of the Toolkit.

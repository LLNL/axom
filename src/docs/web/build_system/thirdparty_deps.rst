Managing Third-party Libraries 
------------------------------
* We are using Spack `<https://github.com/scalability-llnl/spack>`_ to build and manage external dependences for the toolkit.
* To fully automate the process (including obtaining Spack), we are using a small python script "uberenv.py"

Uberenv
^^^^^^^

* This is a thin veneer around Spack, basically a python script.
* Clones spack from githup
* It setup up spack compilers and applies a few spack tweaks, including custom package scripts.
* Calls spack to build a meta-package that builds all of our third-party dependencies and generates a configuration file.

.. image:: ../images/Uberenv.jpg

To create the thirdparty libraries use::

    python scripts/uberenv/uberenv.py --prefix {path to dest} --spec spec  [ --mirror {path} ]

Mirrors
^^^^^^^

To avoid downloading source files for each build, a mirror of the distibution
files can be created for use by spack::

    spack mirror create -d {directory} --dependencies uberenv-asctoolkit

.. https://lc.llnl.gov/confluence/display/ASCT/Third+Party+Dependencies
   https://lc.llnl.gov/confluence/display/ASCT/Toolkit+Development+Environment

Managing Third-party Libraries 
------------------------------
* We are using Spack `<https://github.com/scalability-llnl/spack>`_ to build and manage external dependences for the toolkit.
* To fully automate the process (including obtaining Spack), we are using a small python script "uberenv.py"

Uberenv:
--------
* This is a thin veneer around Spack, basically a python script.
* Clones spack from githup
* It setup up spack compilers and applies a few spack tweaks, including custom package scripts.
* Calls spack to build a meta-package that builds all of our third-party dependencies and generates a configuration file.

.. image:: ../images/Uberenv.jpg




=============================
Current Third Party Packages:
=============================

This page contains a list of the 3rd party packages that we are using in the toolkit + Conduit. 

 

------
C/C++
------

------------------------------------
gperftools: Google Performance tools 
------------------------------------

`<https://github.com/LLNL/conduit/tree/master/src/thirdparty_builtin/gperftools-2.2.1>`_

---------------------------------
gtest: Google unit test framework 
---------------------------------

`<https://github.com/LLNL/conduit/tree/master/src/thirdparty_builtin/gtest-1.7.0>`_

----------------------------
libb64: Base64 encoding
----------------------------

`<https://github.com/LLNL/conduit/tree/master/src/thirdparty_builtin/libb64-1.2.1>`_

-------------------------------
rapidjson: JSON Parsing Library
-------------------------------

`<https://github.com/LLNL/conduit/tree/master/src/thirdparty_builtin/rapidjson>`_

----------------------------
civetweb: Embedded Webserver
----------------------------

`<https://github.com/LLNL/conduit/tree/master/src/thirdparty_builtin/civetweb>`_

(NOTE: The webserver will not be started on the LC without a separate security solution, but we wanted it compiled into conduit)

------------------------------------
google benchmark: micro benchmarking
------------------------------------

`<https://github.com/google/benchmark>`_

----------------------------
boost: C++ utilities 
----------------------------

`<http://www.boost.org>`_

`<http://downloads.sourceforge.net/project/boost/boost/1.58.0/boost_1_58_0.tar.bz2>`_

-------------------------------
google sparse hash: hashing imp
-------------------------------

`<https://code.google.com/p/sparsehash>`_

`<https://sparsehash.googlecode.com/files/sparsehash-2.0.2.tar.gz>`_

----------------------------
lcov: code coverage reports
----------------------------

`<https://ltp.sourceforget.net/coverage/lcov.php>`_

`<http://downloads.sourceforge.net/ltp/lcov-1.11.tar.gz>`_

------
ADIOS
------

`<https://www.olcf.ornl.gov/center-projects/adios/>`_

Download 1.9.0


--------------
Java Script
--------------

----------------------------
d3:  js data vis
----------------------------

`<https://github.com/LLNL/conduit/tree/master/src/libs/conduit_io/web_clients/rest_client/resources/d3>`_

---------------------------
jquery: js utils
---------------------------

`<https://github.com/LLNL/conduit/tree/master/src/libs/conduit_io/web_clients/wsock_test/resources>`_

--------------------------------------
fattable:  js streaming display tables
--------------------------------------

`<https://github.com/LLNL/conduit/tree/master/src/libs/conduit_io/web_clients/rest_client/resources/fattable>`_

--------
Fortran
--------

----------------------------
fruit: fortran unit testing
----------------------------

`<https://github.com/LLNL/conduit/tree/master/src/thirdparty_builtin/fruit-3.3.9>`_


------
Python
------
 

----------------------------
sphinx (for documentation)
----------------------------

`<http://sphinx-doc.org/>`_

`<https://pypi.python.org/packages/source/S/Sphinx/Sphinx-1.3.1.tar.gz#md5=8786a194acf9673464c5455b11fd4332>`_


----------------------------
breathe (doxygen to sphinx)
----------------------------

`<https://breathe.readthedocs.org/en/latest/>`_

`<https://pypi.python.org/packages/source/b/breathe/breathe-4.0.0.tar.gz#md5=32316d5a890a3124ea3e8a9e0b2b3b97>`_

----------------------------
parsely: lang parsing tools
----------------------------

`<https://pypi.python.org/pypi/Parsley>`_

`<https://pypi.python.org/packages/source/P/Parsley/Parsley-1.2.tar.gz>`_


----------------------------
pyyaml: yaml parsing tools
----------------------------

`<http://pyyaml.org/>`_

`<http://pyyaml.org/download/pyyaml/PyYAML-3.11.tar.gz>`_

--------------------------
cog: python code gen
--------------------------

`<http://nedbatchelder.com/code/cog>`_

`<https://pypi.python.org/packages/source/c/cogapp/cogapp-2.4.tar.gz>`_

---------------------------------
setuptools: python package helper
---------------------------------

`<https://pypi.python.org/pypi/setuptools>`_

`<https://pypi.python.org/packages/source/s/setuptools/setuptools-11.3.tar.gz>`_

------------------------------
pybindgen: python c++ bindings 
------------------------------

`<https://launchpad.net/pybindgen>`_

`<https://pypi.python.org/packages/source/P/PyBindGen/PyBindGen-0.17.0.tar.gz>`_
 
==============================
Managing Third-party Libraries 
==============================
* We are using Spack `<https://github.com/scalability-llnl/spack>`_ to build and manage external dependences for the toolkit.
* To fully automate the process (including obtaining Spack), we are using a small python script "uberenv.py"

========
Uberenv:
========
* This is a thin veneer around Spack, basically a python script.
* Clones spack from githup
* It setup up spack compilers and applies a few spack tweaks, including custom package scripts.
* Calls spack to build a meta-package that builds all of our third-party dependencies and generates a configuration file.

.. image:: ../images/Uberenv.jpg



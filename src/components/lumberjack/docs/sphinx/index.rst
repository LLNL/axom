
Lumberjack User Documentation
=============================

Lumberjack is a C++ library that provides scalable logging while reducing
the amount of messages written out the the screen or file system.


Introduction
------------

Lumberjack was created to provide scalable logging with a simple programming
model while allowing developers to customize its behavior. It
uses MPI and a scalable binary tree reduction scheme to combine duplicate
messages and limit output to only the root node.


**Contents:**

.. toctree::
   :maxdepth: 2

   classes


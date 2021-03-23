.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _inlet_writer_page_label:

#######
Writers
#######

Inlet provides a ``Writer`` interface for exporting metadata once the input file has been read in
and the schema has been defined.

Inlet also provides a few implementations of this interface:

------
Sphinx
------

The ``SphinxWriter`` generates a `reStructuredText <https://docutils.sourceforge.io/rst.html>`_ (.rst) file
compatible with the `Sphinx <https://www.sphinx-doc.org/en/master/>`_ documentation tool.

The ``SphinxWriter`` currently has two styles available:

- A table-based style that, for each ``Container``, generates a table for its child fields and functions
- A nested style that includes a "table of contents" for each ``Container`` with links to full descriptions
  of child fields and functions

For comparison, the following is produced from the ``documentation_generation.cpp`` example:

.. toctree::
  :maxdepth: 1

  example1_table_documentation
  example1_nested_documentation

...from the ``mfem_coefficient.cpp`` example:

.. toctree::
  :maxdepth: 1

  mfem_coefficient_table_documentation
  mfem_coefficient_nested_documentation

...and from the ``nested_structs.cpp`` example:

.. toctree::
  :maxdepth: 1

  nested_structs_table_documentation
  nested_structs_nested_documentation

-----------
JSON Schema
-----------

Inlet also provides a utility for generating a `JSON schema <https://json-schema.org/>`_ from your input file schema.
This allows for integration with text editors like Visual Studio Code, which allows you to associate a JSON schema
with an input file and subsequently provides autocompletion, linting, tooltips, and more.  VSCode and other editors
currently support verification of JSON and YAML input files with JSON schemas.

Using the same  ``documentation_generation.cpp`` example, the automatically generated schema can be used to assist
with input file writing:

.. image:: json_schema_example.gif

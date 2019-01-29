#------------------------------------------------------------------------------
# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

#
# Generate source
# quest/src/interface/CMakeLists.txt
#
set(YAML_INPUT_FILE quest_shroud.yaml)

add_shroud(
    YAML_INPUT_FILE      ${YAML_INPUT_FILE}
    C_FORTRAN_OUTPUT_DIR c_fortran
    PYTHON_OUTPUT_DIR    python
)

# Install setup.py file used to build quest Python module.
configure_file(python/setup.py.in ${CMAKE_CURRENT_BINARY_DIR}/python/setup.py)
configure_file(python/quest_test.py.in ${CMAKE_CURRENT_BINARY_DIR}/python/quest_test.py)

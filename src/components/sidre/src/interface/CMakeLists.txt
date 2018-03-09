#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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
# sidre/src/interface/CMakeLists.txt
#
set(YAML_INPUT_FILE sidre_shroud.yaml)

add_shroud(
    YAML_INPUT_FILE      ${YAML_INPUT_FILE}
    YAML_OUTPUT_DIR      yaml
    C_FORTRAN_OUTPUT_DIR c_fortran
    PYTHON_OUTPUT_DIR    python
    DEPENDS_SOURCE
        c_fortran/csidresplicer.c c_fortran/fsidresplicer.f
        python/pysidresplicer.c
    DEPENDS_BINARY genfsidresplicer.f
)

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.f
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/c_fortran/genfsidresplicer.py
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/c_fortran/genfsidresplicer.py fortran > 
            ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.f
)

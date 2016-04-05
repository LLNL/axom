if(ENABLE_SHROUD)
add_shroud( YAML_INPUT_FILE shroud_sidre.yaml
            C_FORTRAN_OUTPUT_DIR c_fortran
            PYTHON_OUTOUT_DIR python
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
endif()

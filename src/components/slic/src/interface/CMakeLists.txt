#
# Generate source
#
if(ENABLE_SHROUD)
  add_shroud( YAML_INPUT_FILE shroud_slic.yaml
              C_FORTRAN_OUTPUT_DIR c_fortran
  )
endif(ENABLE_SHROUD)


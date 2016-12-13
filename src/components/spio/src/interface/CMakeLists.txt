#
# Generate source
#
if(ENABLE_SHROUD)
  add_shroud( YAML_INPUT_FILE spio_shroud.yaml
              C_FORTRAN_OUTPUT_DIR c_fortran
  )
endif()

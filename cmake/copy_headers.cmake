separate_arguments(HEADERS UNIX_COMMAND ${LIBHEADERS})


## Copy headers to the includes folder within the build directory
foreach(hdr ${HEADERS})
   message(STATUS "copy ${hdr}...")
   file(COPY ${hdr} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()
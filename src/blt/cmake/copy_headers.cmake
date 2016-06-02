separate_arguments(HEADERS UNIX_COMMAND ${LIBHEADERS})


## Copy headers to the includes folder within the build directory
foreach(hdr ${HEADERS})
   #message(STATUS "copy ${hdr} to ${HEADER_INCLUDES_DIRECTORY}")
   if(IS_ABSOLUTE ${hdr})
       file(COPY ${hdr} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
   else()
       file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/src/${hdr} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
   endif()
endforeach()
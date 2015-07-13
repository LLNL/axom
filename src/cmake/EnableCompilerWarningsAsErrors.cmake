# - Add flags to treat compiler warnings as errors
#
#  enable_compiler_warnings_as_errors(<targetname>)
#  globally_enable_ompiler_warnings_as_errors() - to modify CMAKE_CXX_FLAGS, etc
#    to change for all targets declared after the command, instead of per-command
#
# Got idea from module for adding extra compiler warnings by Ryan Pavlik
# -- Aaron Black

if(__enable_compiler_warnings_as_errors)
	return()
endif()
set(__enable_compiler_warnings_as_errors YES)

macro(_enable_compiler_warnings_as_errors)
	set(_flags)
	if(MSVC)
		set(_flags "/WX")
	else()
		include(CheckCXXCompilerFlag)
		set(_flags)

		check_cxx_compiler_flag(-Werror SUPPORTS_WERROR_FLAG)
		if(SUPPORTS_WERROR_FLAG)
			set(_flags "${_flags} -Werror")
      else()
         message(WARNING "Compiler does not support -Werror flag, add support for this compiler in EnableCompilerWarningsAsErrors.cmake")
		endif()
	endif()
endmacro()

function(enable_compiler_warnings_as_errors _target)
	_enable_compiler_warnings_as_errors()
	get_target_property(_origflags ${_target} COMPILE_FLAGS)
	if(_origflags)
		set_property(TARGET
			${_target}
			PROPERTY
			COMPILE_FLAGS
			"${_flags} ${_origflags}")
	else()
		set_property(TARGET
			${_target}
			PROPERTY
			COMPILE_FLAGS
			"${_flags}")
	endif()

endfunction()

function(globally_enable_compiler_warnings_as_errors)
	_enable_compiler_warnings_as_errors()
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${_flags}" PARENT_SCOPE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_flags}" PARENT_SCOPE)
endfunction()

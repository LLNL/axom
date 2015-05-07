#!/bin/bash
###############################################################################
# file: config-build.sh
#
# description:
#  Configures an out-of-source debug build of the toolkit, 
#  optionally includes a initial cmake settings file.
#
###############################################################################

#------------------------------------------------------------------------------
#clean up existing build + install directories
#------------------------------------------------------------------------------
rm -rf build-debug install-debug

#------------------------------------------------------------------------------
# create a new build directory
#------------------------------------------------------------------------------
mkdir build-debug
mkdir install-debug
cd build-debug

#------------------------------------------------------------------------------
# setup desired basic cmake options
#------------------------------------------------------------------------------
export CMAKE_OPTS=" -DCMAKE_BUILD_TYPE=Debug"
export CMAKE_OPTS="$CMAKE_OPTS -DCMAKE_INSTALL_PREFIX=../install-debug"

#------------------------------------------------------------------------------
# include an initial cmake settings file if appropriate
#------------------------------------------------------------------------------

# first look for a specific config for this machine
export HOST_CONFIG=../host-configs/`hostname`.cmake
echo "Looking for host-config file: $HOST_CONFIG"
if [[ -e  "$HOST_CONFIG" ]]; then
    echo "FOUND: $HOST_CONFIG"
    export CMAKE_OPTS="$CMAKE_OPTS -C $HOST_CONFIG"
# then check for a sys-type based config
elif [[ "$SYS_TYPE" != "" ]]; then
    export HOST_CONFIG=../host-configs/$SYS_TYPE.cmake
    echo "Looking for SYS_TYPE based host-config file: $HOST_CONFIG"
    if [[ -e  "$HOST_CONFIG" ]]; then
        echo "FOUND: $HOST_CONFIG"
        export CMAKE_OPTS="$CMAKE_OPTS -C $HOST_CONFIG"
    fi
else 
    # fallback to simple a uname based config (Linux / Darwin / etc)
    export HOST_CONFIG=../host-configs/`uname`.cmake
    echo "Looking for uname based host-config file: $HOST_CONFIG"
    if [[ -e  "$HOST_CONFIG" ]]; then
        echo "FOUND: $HOST_CONFIG"
        export CMAKE_OPTS="$CMAKE_OPTS -C $HOST_CONFIG"
    fi
fi

#------------------------------------------------------------------------------
# parse other command line arguments
#------------------------------------------------------------------------------
for arg in "$@"
do
    arg="$1"
    echo "ARGUMENT: $arg"
    case $arg in
        --with-eclipse)
        export CMAKE_ECLIPSE_GEN="Eclipse CDT4 - Unix Makefiles"
        export CMAKE_OPTS="-G \"$CMAKE_ECLIPSE_GEN\" $CMAKE_OPTS"
        shift
        ;;
        --with-warnings)
        export CMAKE_OPTS="-DENABLE_WARNINGS:BOOL=ON $CMAKE_OPTS"
        shift
        ;;        
        # option to dump a compilation database used by the clang tools such as clang-modernize.  
        # This generates a file called compile_commands.json at the top of the build directory
        --with-clang-tooling)
        export CMAKE_OPTS="-DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON $CMAKE_OPTS"
        shift
        ;;        
        *)
        # unknown option, skip
        shift
        ;;
    esac

done


#------------------------------------------------------------------------------
# run cmake to configure
#------------------------------------------------------------------------------
echo "cmake $CMAKE_OPTS ../src"
eval "cmake  $CMAKE_OPTS ../src"

# return to the starting dir
cd ../

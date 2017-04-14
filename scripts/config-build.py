#!/bin/sh
"exec" "python" "-u" "-B" "$0" "$@"

# Python wrapper script for generating the correct cmake line with the options specified by the user.
#
# Please keep parser option names as close to possible as the names of the cmake options they are wrapping.

import sys
import os
import subprocess
import argparse
import platform
import shutil

def extract_cmake_location(file_path):
    # print "Extracting cmake entry from host config file ", file_path
    if os.path.exists(file_path):
        cmake_line_prefix = "# cmake executable path: "
        file_handle = open(file_path, "r")
        content = file_handle.readlines()
        for line in content:
            if line.startswith(cmake_line_prefix):
                return line.split(" ")[4].strip()
        print "Could not find a cmake entry in host config file."
    return None


def parse_arguments():
    parser = argparse.ArgumentParser(description="Configure cmake build.",
                                     epilog="Note: Additional or unrecognized parameters will be passed directly to cmake."
                                            " For example, append '-DENABLE_OPENMP=ON' to enable OpenMP."
                                            " See the toolkit QuickStart guide for a list of available CMake options."
                                            )

    parser.add_argument("-bp",
                        "--buildpath",
                        type=str,
                        default="",
                        help="specify path for build directory.  If not specified, will create in current directory.")

    parser.add_argument("-ip",
                        "--installpath", 
                        type=str, default="",
                        help="specify path for installation directory.  If not specified, will create in current directory.")

    parser.add_argument("-bt",
                        "--buildtype",
                        type=str,
                        choices=["Release", "Debug", "RelWithDebInfo", "MinSizeRel"],
                        default="Debug",
                        help="build type.")

    parser.add_argument("-e",
                        "--eclipse",
                        action='store_true',
                        help="create an eclipse project file.")

    parser.add_argument("-x",
                        "--xcode",
                        action='store_true',
                        help="create an xcode project.")

    parser.add_argument("-ecc",
                        "--exportcompilercommands",
                        action='store_true',
                        help="generate a compilation database.  Can be used by the clang tools such as clang-modernize.  Will create a 'compile_commands.json' file in build directory.")

    parser.add_argument("-hc",
                        "--hostconfig",
                        required=True,
                        type=str,
                        help="select a specific host-config file to initalize CMake's cache")

    args, unknown_args = parser.parse_known_args()
    if unknown_args:
        print "[config-build]: Passing the following arguments directly to cmake... %s" % unknown_args

    return args, unknown_args


########################
# Find CMake Cache File
########################
def find_host_config(args):
    hostconfigpath = os.path.abspath(args.hostconfig)
    assert os.path.exists( hostconfigpath ), "Could not find CMake host config file '%s'." % hostconfigpath
    print "Using host config file: '%s'." % hostconfigpath
    return hostconfigpath


########################
# Get Platform information from host config name
########################
def get_platform_info(hostconfigpath):
    platform_info = ""
    platform_info = os.path.split(hostconfigpath)[1]
    if platform_info.endswith(".cmake"):
        platform_info = platform_info[:-6]
    return platform_info
        

#####################
# Setup Build Dir
#####################
def setup_build_dir(args, platform_info):
    if args.buildpath != "":
        # use explicit build path
        buildpath = args.buildpath
    else:
        # use platform info & build type
        buildpath = "-".join(["build", platform_info, args.buildtype.lower()])

    buildpath = os.path.abspath(buildpath)

    if os.path.exists(buildpath):
        print "Build directory '%s' already exists.  Deleting..." % buildpath
        shutil.rmtree(buildpath)

    print "Creating build directory '%s'..." % buildpath
    os.makedirs(buildpath)
    return buildpath


#####################
# Setup Install Dir
#####################
def setup_install_dir(args, platform_info):    
    # For install directory, we will clean up old ones, but we don't need to create it, cmake will do that.
    if args.installpath != "":
        installpath = os.path.abspath(args.installpath)
    else:
        # use platform info & build type
        installpath = "-".join(["install", platform_info, args.buildtype.lower()])

    installpath = os.path.abspath(installpath)

    if os.path.exists(installpath):
        print "Install directory '%s' already exists, deleting..." % installpath
        shutil.rmtree(installpath)

    print "Creating install path '%s'..." % installpath
    os.makedirs(installpath)
    return installpath


############################
# Build CMake command line
############################
def create_cmake_command_line(args, unknown_args, buildpath, hostconfigpath, installpath):
    cmakeline = extract_cmake_location(hostconfigpath)
    assert cmakeline, "Host config file doesn't contain valid cmake location, value was %s" % cmakeline

    # Add cache file option
    cmakeline += " -C %s" % hostconfigpath
    # Add build type (opt or debug)
    cmakeline += " -DCMAKE_BUILD_TYPE=" + args.buildtype
    # Set install dir
    cmakeline += " -DCMAKE_INSTALL_PREFIX=%s" % installpath

    if args.exportcompilercommands:
        cmakeline += " -DCMAKE_EXPORT_COMPILE_COMMANDS=on"

    if args.eclipse:
        cmakeline += ' -G "Eclipse CDT4 - Unix Makefiles"'

    if args.xcode:
        cmakeline += ' -G "Xcode"'

    if unknown_args:
        cmakeline += " " + " ".join( unknown_args )

    scriptsdir = os.path.dirname( os.path.abspath(sys.argv[0]) )
    cmakeline += " %s/../src " % scriptsdir

    # Dump the cmake command to file for convenience
    with open("%s/cmake_cmd" % buildpath, "w") as cmdfile:
       cmdfile.write(cmakeline)

    import stat
    st = os.stat("%s/cmake_cmd" % buildpath)
    os.chmod("%s/cmake_cmd" % buildpath, st.st_mode | stat.S_IEXEC)
    return cmakeline


############################
# Run CMake
############################
def run_cmake(buildpath, cmakeline):
    print "Changing to build directory..."
    os.chdir(buildpath)
    print "Executing CMake line: '%s'" % cmakeline
    print 
    returncode = subprocess.call(cmakeline, shell=True)
    if not returncode == 0:
        print "Error: CMake command failed with return code: {0}".format(returncode)
        return False
    return True


############################
# Main
############################
def main():
    args, unknown_args = parse_arguments()

    hostconfigpath = find_host_config(args)
    platform_info = get_platform_info(hostconfigpath)
    buildpath = setup_build_dir(args, platform_info)
    installpath = setup_install_dir(args, platform_info)

    cmakeline = create_cmake_command_line(args, unknown_args, buildpath, hostconfigpath, installpath)
    return run_cmake(buildpath, cmakeline)

if __name__ == '__main__':
    exit(0 if main() else 1)

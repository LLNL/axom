#!/usr/bin/env python
# Python wrapper script for generating the correct cmake line with the options specified by the user.
#
# Please keep parser option names as close to possible as the names of the cmake options they are wrapping.

import sys
import os
import subprocess
import argparse
import platform
import shutil

# Helper function to get SYS_TYPE on LC systems.
def get_systype():
    import os
    lc_home_config_filename = "/etc/home.config"
    if os.path.exists(lc_home_config_filename):
        file_handle = open(lc_home_config_filename, "r")
        content = file_handle.readlines()
        for line in content:
            if line.startswith("SYS_TYPE"):
                return line.split(" ")[1].strip()
    return None

def extract_cmake_location(file_path):
    cmake_line_prefix = "# cmake executable path: "
    if os.path.exists(file_path):
        file_handle = open(file_path, "r")
        content = file_handle.readlines()
        for line in content:
            if line.startswith(cmake_line_prefix):
                return line.split(" ")[4].strip()
        print "Could not find a cmake entry in host config file."
        return None


default_compiler = "gcc@4.9.3"

parser = argparse.ArgumentParser(description="Configure cmake build.")

parser.add_argument("-bp",
                    "--buildpath",
                    type=str,
                    default="",
                    help="specify path for build directory.  If not specified, will create one under current directory.")

parser.add_argument("-ip",
                    "--installpath", 
                    type=str, default="",
                    help="specify path for installation directory.  If not specified, will create one under current directory.")

parser.add_argument("-c",
                    "--compiler",
                    type=str,
                    default=default_compiler,
                    help="compiler to use.")

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
	                help="generate a compilation database.  Can be used by the clang tools such as clang-modernize.  Will create a file called 'compile_commands.json' in your build directory.")

parser.add_argument("-co",
                    "--cmakeoption",
                    type=str,
                    help="specify additional cmake option to add to cmake line.  A '-D' will automatically be prepended.  Use caution, if you are doing something non-trivial, use ccmake or cmake-gui.")

parser.add_argument("-hc",
                    "--hostconfig",
                    default="",
                    type=str,
                    help="select a specific host-config file to initalize CMake's cache")

args = parser.parse_args()

########################
# Find CMake Cache File
########################
platform_info = ""
scriptsdir = os.path.dirname( os.path.abspath(sys.argv[0]) )
systype = get_systype()

if args.hostconfig != "":
    cachefile = os.path.abspath(args.hostconfig)
    platform_info = os.path.split(cachefile)[1]
    if platform_info.endswith(".cmake"):
        platform_info = platform_info[:-6]
    print "Using user specified host config file: '%s'." % cachefile
else:
    # If not specified, then check for a host-config file for this SYS_TYPE with default compiler.
    cachefile = scriptsdir.replace("scripts","host-configs")
    if systype:
        platform_info = systype.split("_")[0]
        import glob
        names = glob.glob(cachefile + "/*" + systype + "-" + args.compiler + ".cmake")
        assert len(names) <= 1, "Could not determine correct host-config file for SYS_TYPE %s, more than one file matched this SYS_TYPE and compiler." % systype
        assert len(names) > 0, "Could not find host-config file for SYS_TYPE %s and compiler %s" % (systype, args.compiler)
        cachefile = os.path.join( cachefile, names[0] )
        print "Found host config file for SYS_TYPE %s, compiler %s: %s" % (systype, args.compiler, cachefile)
    else:
        platform_info = platform.node()
        cachefile = os.path.join(cachefile, "other", "%s.cmake" % platform_info )
        print "No /etc/home.config file found, must not be a LC system.  Using hostname config file: '%s'" % ( cachefile )

assert os.path.exists( cachefile ), "Could not find cmake cache file '%s'." % cachefile


#####################
# Setup Build Dir
#####################
if args.buildpath != "":
    # use explicit build path
    buildpath = args.buildpath
elif args.hostconfig != "":
    # use host config name (via platform_info) as build dir base
    buildpath = "-".join(["build",platform_info,args.buildtype.lower()]).replace("@","-")    
elif args.buildpath == "":
    # Generate build directory name based on platform, buildtype, compiler
    buildpath = "-".join(["build",platform_info, args.compiler, args.buildtype.lower()]).replace("@","-")

buildpath = os.path.abspath(buildpath)

if os.path.exists(buildpath):
    print "Build directory '%s' already exists.  Deleting..." % buildpath
    shutil.rmtree(buildpath)

print "Creating build directory '%s'..." % buildpath
os.makedirs(buildpath)

#####################
# Setup Install Dir
#####################
# For install directory, we will clean up old ones, but we don't need to create it, cmake will do that.
if args.installpath != "":
    installpath = os.path.abspath(args.installpath)
elif args.hostconfig != "":
    # use host config name (via platform_info) as install dir base
    installpath = "-".join(["install",platform_info,args.buildtype.lower()]).replace("@","-")
else:
    # Generate install directory name based on platform, buildtype, compiler
    installpath = "-".join(["install",platform_info, args.compiler, args.buildtype.lower()]).replace("@","-")

installpath = os.path.abspath(installpath)

if os.path.exists(installpath):
    print "Install directory '%s' already exists, deleting..." % installpath
    import shutil
    shutil.rmtree(installpath)

print "Creating install path '%s'..." % installpath
os.makedirs(installpath)

############################
# Build CMake command line
############################

cmakeline = extract_cmake_location(cachefile)
assert cmakeline, "Host config file doesn't contain valid cmake location, value was %s" % cmakeline

# Add cache file option
cmakeline += " -C %s" % cachefile
# Add build type (opt or debug)
cmakeline += " -DCMAKE_BUILD_TYPE=" + args.buildtype
# Set install dir
cmakeline += " -DCMAKE_INSTALL_PREFIX=%s" % installpath

if args.exportcompilercommands:
    cmakeline += " -DCMAKE_EXPORT_COMPILE_COMMANDS=on"

if args.eclipse:
    cmakeline += ' -G "Eclipse CDT4 - Unix Makefiles"'

if args.xcode:
    cmakeline += ' -G "XCode"'

if args.cmakeoption:
    cmakeline += " -D" + args.cmakeoption

cmakeline += " %s/../src " % scriptsdir

# Dump the cmake command to file for convenience
cmdfile = open("%s/cmake_cmd" % buildpath, "w")
cmdfile.write(cmakeline)
cmdfile.close()
import stat
st = os.stat("%s/cmake_cmd" % buildpath)
os.chmod("%s/cmake_cmd" % buildpath, st.st_mode | stat.S_IEXEC)

############################
# Run CMake
############################
print "Changing to build directory..."
os.chdir(buildpath)
print "Executing cmake line: '%s'" % cmakeline
print 
subprocess.call(cmakeline,shell=True)


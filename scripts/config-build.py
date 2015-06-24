#!/usr/bin/env python
# Python wrapper script for generating the correct cmake line with the options specified by the user.
#
# Please keep parser option names as close to possible as the names of the cmake options they are wrapping.

import argparse

parser = argparse.ArgumentParser(description="Configure cmake build.")

parser.add_argument("-bp", "--buildpath", type=str, default="", help="specify path for build directory.  If not specified, will create one under current directory.")
parser.add_argument("-ip", "--installpath", type=str, default="", help="specify path for installation directory.  If not specified, will create one under current directory.")
parser.add_argument("-e", "--eclipse", action='store_true', help="create an eclipse project file.")
parser.add_argument("-c", "--compiler", type=str, choices=["intel","gnu","clang","xl"], default="gnu", help="compiler to use.")
parser.add_argument("-bt", "--buildtype", type=str, choices=["Release", "Debug", "RelWithDebInfo", "MinSizeRel"], default="Debug", help="build type.")
parser.add_argument("-ecc", "--exportcompilercommands", action='store_true',
	 help="generate a compilation database.  Can be used by the clang tools such as clang-modernize.  Will create a file called 'compile_commands.json' in your build directory.")
parser.add_argument("-co", "--cmakeoption", type=str, help="specify additional cmake options to add to cmake line.  Use with caution, if you are doing something non-trivial, use ccmake or cmake-gui.")

args = parser.parse_args()

import sys
import os
scriptsdir = os.path.dirname( os.path.join( os.getcwd(), sys.argv[0] ) )

cmakeline = "cmake"

# Add build type (opt or debug)
cmakeline += " -DCMAKE_BUILD_TYPE=" + args.buildtype

# Check if 'SYS_TYPE' exists, and look for cache file there.
cachefile = scriptsdir.replace("scripts","host-configs")

shortsystype = ""
if "SYS_TYPE" in os.environ:
    systype = os.environ["SYS_TYPE"]
    shortsystype = systype.split("_")[0]
    cachefile = os.path.join( cachefile, shortsystype, "%s.cmake" % args.compiler ) 
else:
    import platform
    shortsystype = platform.node()
    cachefile = os.path.join(cachefile, "other", "%s.cmake" % shortsystype )

assert os.path.exists( cachefile ), "Could not find cmake cache file '%s'." % cachefile
cmakeline += " -C %s" % cachefile

if args.buildpath == "":
    # Generate build directory name based on platform, buildtype, compiler
    buildpath = "-".join([shortsystype, args.compiler, args.buildtype])
else:
    buildpath = args.buildpath

if os.path.exists(buildpath):
    print "Build directory '%s' already exists.  Deleting..." % buildpath
    import shutil
    shutil.rmtree(buildpath)
print "Creating build directory '%s'..." % buildpath
os.makedirs(buildpath)

# For install directory, we will clean up old ones, but we don't need to create it, cmake will do that.
if args.installpath == "":
    # Generate install directory name based on platform, buildtype, compiler
    installpath = os.path.abspath("-".join([shortsystype, args.compiler, args.buildtype,"install"]))
else:
    installpath = os.path.abspath(args.installpath)

if os.path.exists(installpath):
    print "Install directory '%s' already exists, deleting..." % installpath
    import shutil
    shutil.rmtree(installpath)
cmakeline += " -DCMAKE_INSTALL_PREFIX=%s" % installpath
print "Creating install path '%s'..." % installpath
os.makedirs(installpath)

if args.exportcompilercommands:
    cmakeline += " -DCMAKE_EXPORT_COMPILE_COMMANDS=on"

if args.eclipse:
    cmakeline += ' -G "Eclipse CDT4 - Unix Makefiles"'

if args.cmakeoption:
    cmakeline += " " + parser.extraoptions

cmakeline += " %s/../src " % scriptsdir

print "Changing to build directory..."
os.chdir(buildpath)
print "Executing cmake line: '%s'..." % cmakeline
print 
os.system(cmakeline)

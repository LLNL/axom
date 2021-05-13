#!/usr/bin/env python

# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
# 
# SPDX-License-Identifier: (BSD-3-Clause)

# Python script to add Axom LLNL copyright notice at top of files in a directory
#
# The script takes a directory and checks all files in the directory.
# If the second line in the file does not match the copyright string, we add it.
#
# Modified from an initial script by P. Sinha

import os
import sys
import argparse

axom_copyright_str = """// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
\n\n"""

axom_copyright_begin_str = "Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and"

def checkAndAddCopyrightHeader(filename, testOnly=False):

  with open(filename, "r+") as f:
    first_line = f.readline()   # First line is only a c-style comment opener
    second_line = f.readline()  # So, we check the second line
    
    print "  Processing file:", filename,
    
    if not axom_copyright_begin_str in second_line:
        if testOnly:
            print "\t missing copyright statement."
        else:
            lines = f.readlines()
            f.seek(0)
            f.write(axom_copyright_str)
            f.write(first_line)
            f.write(second_line)
            f.writelines(lines)
            print "\t prepended copyright statement."
    else:
        print "\t already has copyright statement."


def fileNameGenerator(rootDir, validExtensions, isRecursive=False):
    """Generator function for file names whose extensions are in the validExtensions tuple.
       Files are rooted in rootDir, and process is recursive if isRecursive==True.
    """

    if isRecursive:
        for path,dirlist,filelist in os.walk(rootDir):
            for f in (f for f in filelist if f.lower().endswith(validExtensions) ):
                yield os.path.join(path, f)
    else:
        for f in os.listdir(rootDir):
            if f.lower().endswith( validExtensions):
                yield os.path.join(rootDir, f)

if __name__ == "__main__":

    ## Setup the argument parser, dir is required first argument
    parser = argparse.ArgumentParser(description="Append LLNL copyright message to files.")
    parser.add_argument("dir", type=str, help="specify directory containing files on which we want to operate.") # TODO -- should we accept multiple directories?

    # Option to recursively search for files
    parser.add_argument("-r", "--recursive", dest='isRecursive', action='store_true', help="add flag to recursively descend to subdirectories.")
    parser.add_argument("--no-recursive", dest='isRecursive', action='store_false')
    parser.set_defaults(isRecursive=False)

    # Test run to see which files might require a copyright notice.
    parser.add_argument("-t", "--test", action='store_true', help="add flag if we only want to see which files are missing copyright statements, but not to modify any files.") 

    # Additional possible featues to be implemented
    #
    # Specify a collection of file types.  
    #parser.add_argument("-f", "--filetypes", type=str, nargs="*", help="add file types on which to operate ")    # default should be *.hpp and *.cpp
    # For now, we hardcode to c/c++ header and source files
    valid_extensions = (".hpp", ".cpp", ".h", ".c", ".cxx", ".cc")
    
    args = parser.parse_args()

    ## Iterate through files, check for and add copyright notice
    print "Looking at directory {}".format( args.dir )   
    for fullFileName in fileNameGenerator(args.dir, valid_extensions, args.isRecursive):
        checkAndAddCopyrightHeader(fullFileName, args.test)


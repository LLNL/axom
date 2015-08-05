#!/usr/bin/env python
# Python script to add the ATK LLNL copyright notice to the top of files in a directory
#
# The script takes a directory and checks all files in the directory.
# If the second line in the file does not match the copyright string, we add it.
#
# Modified from an initial script by P. Sinha

import os
import sys
import argparse

atk_copyright_str = """/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */
\n\n"""

atk_copyright_begin_str = "Copyright (c) 2015, Lawrence Livermore National Security, LLC."

def checkAndAddCopyrightHeader(filename):

  with open(filename, "r+") as f:
    first_line = f.readline()   # First line is only a c-style comment opener
    second_line = f.readline()  # So, we check the second line
    
    print "  Processing file:", filename,
    
    if not atk_copyright_begin_str in second_line:
        lines = f.readlines()
        f.seek(0)
        f.write(atk_copyright_str)
        f.write(first_line)
        f.write(second_line)
        f.writelines(lines)
        print "\tadded copyright."
    else:
        print "\t<n/a>."


if __name__ == "__main__":

    # Setup the argument parser
    parser = argparse.ArgumentParser(description="Append LLNL copyright message to files.")
    parser.add_argument("dir", type=str, help="specify directory containing files on which we want to operate.") # TODO -- should we accept multiple directories?

    # Additional possible featues to be implemented
    #
    # Recursively search for files
    #parser.add_argument("-r", "--recursive", dest='isRecursive', action='store_true', help="add flag to recursively descend to subdirectories.")
    #parser.add_argument("--no-recursive", dest='isRecursive', action='store_false')
    #parser.set_defaults(isRecursive=False)
    #
    # Specify a collection of file types.  For now, we hardcode to hpp and cpp files.
    #parser.add_argument("-f", "--filetypes", type=str, nargs="*", help="add file types on which to operate ")    # default should be *.hpp and *.cpp
    #
    # Test run to see which files might require a copyright notice.
    #parser.add_argument("-t", "--testRun", action='store_true', help="add flag if we only want to see which files are missing copyright statements, but not to modify any files.") 
    
    args = parser.parse_args()

    valid_extensions = (".hpp", ".cpp", ".h", ".c", ".cxx")

    dir_of_interest = args.dir
    print "Looking at directory {}".format( dir_of_interest)   
    
    cpp_src_files = (f for f in os.listdir(dir_of_interest) if f.lower().endswith( valid_extensions) )
    for f in cpp_src_files:
        checkAndAddCopyrightHeader(os.path.join(dir_of_interest,f))
        

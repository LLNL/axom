#!/usr/bin/env python
# Python wrapper script for generating the correct cmake line with the options specified by the user.
#
import os

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

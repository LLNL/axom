#!/bin/sh
"exec" "python3" "-u" "-B" "$0" "$@"

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

"""
 file: build_src.py

 description:
  Builds all Axom with the host-configs for the current machine.

"""

from llnl_lc_build_tools import *

from argparse import ArgumentParser


def parse_args():
    "Parses args from command line"
    parser = ArgumentParser()
    # Location of source directory to build
    parser.add_argument("-d", "--directory",
                        dest="directory",
                        default="",
                        help="Directory of source to be built (Defaults to current)")
    # Whether to build a specific hostconfig
    parser.add_argument("--host-config",
                        dest="hostconfig",
                        default="",
                        help="Specific host-config file to build (Tries multiple known paths to locate given file)")
    # Build type for the configuration
    parser.add_argument("--build-type",
                        dest="buildtype",
                        default="Debug",
                        choices = ("Debug", "RelWithDebInfo", "Release", "MinSizeRel"),
                        help="The CMake build type to use")
    # Run unit tests serially (MPI Bug on El Capitan)
    parser.add_argument("--test-serial",
                        action="store_true",
                        dest="testserial",
                        default=False,
                        help="Run unit tests serially")
    # Extra cmake options to pass to config build
    parser.add_argument("--extra-cmake-options",
                        dest="extra_cmake_options",
                        default="",
                        help="Extra cmake options to add to the cmake configure line")
    parser.add_argument("--automation-mode",
                        action="store_true",
                        dest="automation",
                        default=False,
                        help="Toggle automation mode which uses env $HOST_CONFIG then $SYS_TYPE/$COMPILER if found")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        default=False,
                        help="Output logs to screen as well as to files")
    ###############
    # parse args
    ###############
    opts = parser.parse_args()
    # we want a dict b/c the values could
    # be passed without using argparse
    opts = vars(opts)

    # Ensure correctness
    if opts["automation"] and opts["hostconfig"] != "":
        print("[ERROR: automation and host-config modes are mutually exclusive]")
        sys.exit(1)

    return opts


def main():
    opts = parse_args()

    # Determine source directory to be built
    if os.environ.get("UBERENV_PREFIX") != None:
        repo_dir = os.environ["UBERENV_PREFIX"]
        if not os.path.isdir(repo_dir):
            print("[ERROR: Given environment variable 'UBERENV_PREFIX' is not a valid directory]")
            print("[    'UBERENV_PREFIX' = %s]" % repo_dir)
            return 1
    if opts["directory"] != "":
        repo_dir = opts["directory"]
        if not os.path.isdir(repo_dir):
            print("[ERROR: Given command line variable '--directory' is not a valid directory]")
            print("[    '--directory' = %s]" % repo_dir)
            return 1
    else:
        repo_dir = get_repo_dir()

    try:
        original_wd = os.getcwd()
        os.chdir(repo_dir)
        timestamp = get_timestamp()

        # Default to build all SYS_TYPE's host-configs in host-config/
        build_all = not opts["hostconfig"] and not opts["automation"]
        if build_all:
            res = build_and_test_host_configs(repo_dir, timestamp, False,
                                              report_to_stdout = opts["verbose"],
                                              extra_cmake_options = opts["extra_cmake_options"],
                                              build_type = opts["buildtype"],
                                              test_serial = opts["testserial"])
        # Otherwise try to build a specific host-config
        else:
            # Command-line arg has highest priority
            if opts["hostconfig"]:
                hostconfig = opts["hostconfig"]

            # Otherwise try to reconstruct host-config path from SYS_TYPE and COMPILER
            elif opts["automation"]:
                if not "SYS_TYPE" in os.environ:
                    print("[ERROR: Automation mode required 'SYS_TYPE' environment variable]")
                    return 1
                if not "COMPILER" in os.environ:
                    print("[ERROR: Automation mode required 'COMPILER' environment variable]")
                    return 1
                import socket
                hostname = socket.gethostname()
                # Remove any numbers after the end
                hostname = hostname.rstrip('0123456789')
                sys_type = os.environ["SYS_TYPE"]
                # Remove everything including and after the last hyphen
                sys_type = sys_type.rsplit('-', 1)[0]
                compiler = os.environ["COMPILER"]
                compiler = compiler.rsplit('-', 1)[0]
                hostconfig = "%s-%s-%s.cmake" % (hostname, sys_type, compiler)

            # First try with where uberenv generates host-configs.
            hostconfig_path = os.path.join(repo_dir, hostconfig)
            if not os.path.isfile(hostconfig_path):
                print("[INFO: Looking for hostconfig at %s]" % hostconfig_path)
                print("[WARNING: Spack generated host-config not found, trying with predefined]")

                # Then look into project predefined host-configs.
                hostconfig_path = os.path.join(repo_dir, "host-configs", hostconfig)
                if not os.path.isfile(hostconfig_path):
                    print("[INFO: Looking for hostconfig at %s]" % hostconfig_path)
                    print("[WARNING: Predefined host-config not found, trying with Docker]")

                    # Otherwise look into project predefined Docker host-configs.
                    hostconfig_path = os.path.join(repo_dir, "host-configs", "docker", hostconfig)
                    if not os.path.isfile(hostconfig_path):
                        print("[INFO: Looking for hostconfig at %s]" % hostconfig_path)
                        print("[WARNING: Predefined Docker host-config not found]")
                        print("[ERROR: Could not find any host-configs in any known path. Try giving fully qualified path.]")
                        return 1

            print("[build info]")
            binfo_str = json.dumps(build_info(),indent=2)
            print(binfo_str)

            test_root = get_build_and_test_root(repo_dir, timestamp)
            test_root = "{0}_{1}".format(test_root, hostconfig.replace(".cmake", "").replace("@","_"))
            os.mkdir(test_root)
            res = build_and_test_host_config(test_root, hostconfig_path,
                                             report_to_stdout = opts["verbose"],
                                             extra_cmake_options = opts["extra_cmake_options"],
                                             build_type = opts["buildtype"],
                                             test_serial = opts["testserial"])

    finally:
        os.chdir(original_wd)

    return res

if __name__ == "__main__":
    sys.exit(main())

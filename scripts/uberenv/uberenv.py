#!/bin/sh
"exec" "python" "-u" "-B" "$0" "$@"
###############################################################################
# Copyright (c) 2014-2019, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-666778
#
# All rights reserved.
#
# This file is part of Conduit.
#
# For details, see https://lc.llnl.gov/conduit/.
#
# Please also read conduit/LICENSE
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the disclaimer below.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the name of the LLNS/LLNL nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
# LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
###############################################################################

"""
 file: uberenv.py

 description: automates using spack to install a project.

"""

import os
import sys
import subprocess
import shutil
import socket
import platform
import json
import datetime
import glob

from optparse import OptionParser

from os import environ as env
from os.path import join as pjoin


def sexe(cmd,ret_output=False,echo = False):
    """ Helper for executing shell commands. """
    if echo:
        print("[exe: {}]".format(cmd))
        sys.stdout.flush()
    if ret_output:
        p = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        res = p.communicate()[0]
        res = res.decode('utf8')
        return p.returncode,res
    else:
        return subprocess.call(cmd,shell=True)


def parse_args():
    "Parses args from command line"
    parser = OptionParser()
    parser.add_option("--install",
                      action="store_true",
                      dest="install",
                      default=False,
                      help="Install `package_name` instead of `uberenv_package_name`.")
    # where to install
    parser.add_option("--prefix",
                      dest="prefix",
                      default="uberenv_libs",
                      help="destination directory")
    # what compiler to use
    parser.add_option("--spec",
                      dest="spec",
                      default=None,
                      help="spack compiler spec")
    # for vcpkg, what architecture to target
    parser.add_option("--triplet",
                      dest="triplet",
                      default=None,
                      help="vcpkg architecture triplet")
    # optional location of spack mirror
    parser.add_option("--mirror",
                      dest="mirror",
                      default=None,
                      help="spack mirror directory")
    # flag to create mirror
    parser.add_option("--create-mirror",
                      action="store_true",
                      dest="create_mirror",
                      default=False,
                      help="Create spack mirror")
    # this option allows a user to explicitly to select a
    # group of spack settings files (compilers.yaml , packages.yaml)
    parser.add_option("--spack-config-dir",
                      dest="spack_config_dir",
                      default=None,
                      help="dir with spack settings files (compilers.yaml, packages.yaml, etc)")

    # a file that holds settings for a specific project
    # using uberenv.py
    parser.add_option("--project-json",
                      dest="project_json",
                      default=pjoin(uberenv_script_dir(),"project.json"),
                      help="uberenv project settings json file")

    # flag to use insecure curl + git
    parser.add_option("-k",
                      action="store_true",
                      dest="ignore_ssl_errors",
                      default=False,
                      help="Ignore SSL Errors")

    # option to force a spack pull
    parser.add_option("--pull",
                      action="store_true",
                      dest="repo_pull",
                      default=False,
                      help="Pull if spack repo already exists")

    # option to force for clean of packages specified to
    # be cleaned in the project.json
    parser.add_option("--clean",
                      action="store_true",
                      dest="spack_clean",
                      default=False,
                      help="Force uninstall of packages specified in project.json")

    # option to tell spack to run tests
    parser.add_option("--run_tests",
                      action="store_true",
                      dest="run_tests",
                      default=False,
                      help="Invoke build tests during spack install")

    # option to init osx sdk env flags
    parser.add_option("--macos-sdk-env-setup",
                      action="store_true",
                      dest="macos_sdk_env_setup",
                      default=False,
                      help="Set several env vars to select OSX SDK settings."
                           "This was necessary for older versions of macOS "
                           " but can cause issues with macOS versions >= 10.13. "
                           " so it is disabled by default.")


    ###############
    # parse args
    ###############
    opts, extras = parser.parse_args()
    # we want a dict b/c the values could
    # be passed without using optparse
    opts = vars(opts)
    if not opts["spack_config_dir"] is None:
        opts["spack_config_dir"] = os.path.abspath(opts["spack_config_dir"])
        if not os.path.isdir(opts["spack_config_dir"]):
            print("[ERROR: invalid spack config dir: {} ]".format(opts["spack_config_dir"]))
            sys.exit(-1)
    return opts, extras


def uberenv_script_dir():
    # returns the directory of the uberenv.py script
    return os.path.dirname(os.path.abspath(__file__))

def load_json_file(json_file):
    # reads json file
    return json.load(open(json_file))

def is_darwin():
    return "darwin" in platform.system().lower()

def is_windows():
    return "windows" in platform.system().lower()

class UberEnv():
    """ Base class for package manager """

    def __init__(self, opts, extra_opts):
        self.opts = opts
        self.extra_opts = extra_opts

        # load project settings
        self.project_opts = load_json_file(opts["project_json"])
        print("[uberenv project settings: {}]".format(str(self.project_opts)))
        print("[uberenv options: {}]".format(str(self.opts)))

        # setup main package name
        if opts["install"]:
            self.pkg_name = self.project_opts["package_name"]
        else:
            self.pkg_name = self.project_opts["uberenv_package_name"]

    def setup_paths_and_dirs(self):
        self.uberenv_path = os.path.split(os.path.abspath(__file__))[0]

        self.dest_dir = os.path.abspath(self.opts["prefix"])
        print("[installing to: {0}]".format(self.dest_dir))

        # print a warning if the dest path already exists
        if not os.path.isdir(self.dest_dir):
            os.mkdir(self.dest_dir)
        else:
            print("[info: destination '{}' already exists]".format(self.dest_dir))

    def detect_platform(self):
        # find supported sets of compilers.yaml, packages,yaml
        res = None
        if is_darwin():
            res = "darwin"
        elif "SYS_TYPE" in os.environ.keys():
            sys_type = os.environ["SYS_TYPE"].lower()
            res = sys_type
        return res


class VcpkgEnv(UberEnv):
    """ Helper to clone vcpkg and install libraries on Windows """

    def __init__(self, opts, extra_opts):
        UberEnv.__init__(self,opts,extra_opts)

        # setup architecture triplet
        self.triplet = opts["triplet"]
        if self.triplet is None:
           self.triplet = os.getenv("VCPKG_DEFAULT_TRIPLET", "x86-windows")

    def setup_paths_and_dirs(self):
        # get the current working path, and the glob used to identify the
        # package files we want to hot-copy to vcpkg

        UberEnv.setup_paths_and_dirs(self)

        self.ports = pjoin(self.uberenv_path, "vcpkg_ports","*")

        # setup path for vcpkg repo
        self.dest_vcpkg = pjoin(self.dest_dir,"vcpkg")

        if os.path.isdir(self.dest_vcpkg):
            print("[info: destination '{}' already exists]".format(self.dest_vcpkg))

    def clone_repo(self):
        if not os.path.isdir(self.dest_vcpkg):
            # compose clone command for the dest path, vcpkg url and branch
            vcpkg_branch = self.project_opts.get("vcpkg_branch", "master")
            vcpkg_url = self.project_opts.get("vcpkg_url", "https://github.com/microsoft/vcpkg")

            print("[info: cloning vcpkg '{}' branch from {} into {}]"
                .format(vcpkg_branch,vcpkg_url, self.dest_vcpkg))

            os.chdir(self.dest_dir)

            clone_opts = ("-c http.sslVerify=false " 
                          if self.opts["ignore_ssl_errors"] else "")

            clone_cmd =  "git {} clone -b {} {}".format(clone_opts, vcpkg_branch,vcpkg_url)
            sexe(clone_cmd, echo=True)

            # optionally, check out a specific commit
            if "vcpkg_commit" in self.project_opts:
                sha1 = self.project_opts["vcpkg_commit"]
                print("[info: using vcpkg commit {}]".format(sha1))
                os.chdir(self.dest_vcpkg)
                sexe("git checkout {}".format(sha1),echo=True)
                
        if self.opts["repo_pull"]:
            # do a pull to make sure we have the latest
            os.chdir(self.dest_vcpkg)
            sexe("git stash", echo=True)
            sexe("git pull", echo=True)

        # Bootstrap vcpkg
        os.chdir(self.dest_vcpkg)
        sexe("bootstrap-vcpkg.bat -disableMetrics", echo=True)
            

    def patch(self):
        """ hot-copy our ports into vcpkg """
        
        import distutils
        from distutils import dir_util

        src_vcpkg_ports = pjoin(self.uberenv_path, "vcpkg_ports")
        dest_vcpkg_ports = pjoin(self.dest_vcpkg,"ports")

        print("[info: copying from {} to {}]".format(src_vcpkg_ports,dest_vcpkg_ports))
        distutils.dir_util.copy_tree(src_vcpkg_ports,dest_vcpkg_ports)


    def clean_build(self):
        pass

    def show_info(self):
        os.chdir(self.dest_vcpkg)
        print("[info: Details for package '{}']".format(self.pkg_name))
        sexe("vcpkg.exe search " + self.pkg_name, echo=True)

        print("[info: Dependencies for package '{}']".format(self.pkg_name))
        sexe("vcpkg.exe depend-info " + self.pkg_name, echo=True)

    def create_mirror(self):
        pass

    def use_mirror(self):
        pass

    def install(self):
        
        os.chdir(self.dest_vcpkg)
        install_cmd = "vcpkg.exe "
        install_cmd += "install {0}:{1}".format(self.pkg_name, self.triplet)

        res = sexe(install_cmd, echo=True)

        # Running the install_cmd eventually generates the host config file.
        # Copy this to the target directory.
        src_hc = pjoin(self.dest_vcpkg, "installed", self.triplet, "include", self.pkg_name, "hc.cmake")
        hcfg_fname = pjoin(self.dest_dir, "{0}.{1}.cmake".format(platform.uname()[1], self.triplet))
        print("[info: copying host config file to {}]".format(hcfg_fname))
        shutil.copy(os.path.abspath(src_hc), hcfg_fname)
        print("")
        print("[install complete!]")
        return res


class SpackEnv(UberEnv):
    """ Helper to clone spack and install libraries on MacOS an Linux """

    def __init__(self, opts, extra_opts):
        UberEnv.__init__(self,opts,extra_opts)

        # Some additional setup for macos
        if is_darwin():
            if opts["macos_sdk_env_setup"]:
                # setup osx deployment target and sdk settings
                setup_osx_sdk_env_vars()
            else:
                print("[skipping MACOSX env var setup]")

        # setup default spec
        if opts["spec"] is None:
            if is_darwin():
                opts["spec"] = "%clang"
            else:
                opts["spec"] = "%gcc"
        print("[spack spec: {}]".format(opts["spec"]))

    def setup_paths_and_dirs(self):
        # get the current working path, and the glob used to identify the
        # package files we want to hot-copy to spack

        UberEnv.setup_paths_and_dirs(self)

        self.pkgs = pjoin(self.uberenv_path, "packages","*")

        # setup path for spack 
        self.dest_spack = pjoin(self.dest_dir,"spack")

        if os.path.isdir(self.dest_spack):
            print("[info: destination '{}' already exists]".format(self.dest_spack))


    def find_spack_pkg_path(self,pkg_name):
        r,rout = sexe("spack/bin/spack find -p " + pkg_name,ret_output = True)
        for l in rout.split("\n"):
            if l.startswith(" "):
                return {"name": pkg_name, "path": l.split()[-1]}

    def read_spack_full_spec(self,pkg_name,spec):
        rv, res = sexe("spack/bin/spack spec " + pkg_name + " " + spec, ret_output=True)
        for l in res.split("\n"):
            if l.startswith(pkg_name) and l.count("@") > 0 and l.count("arch=") > 0:
                return l.strip()

    def clone_repo(self):
        if not os.path.isdir(self.dest_spack):
            # compose clone command for the dest path, spack url and branch
            spack_branch = self.project_opts.get("spack_branch", "develop")
            spack_url = self.project_opts.get("spack_url", "https://github.com/spack/spack.git")

            print("[info: cloning spack '{}' branch from {} into {}]"
                .format(spack_branch, spack_url, self.dest_spack))

            os.chdir(self.dest_dir)

            clone_opts = ("-c http.sslVerify=false " 
                          if self.opts["ignore_ssl_errors"] else "")


            clone_cmd =  "git {} clone -b {} {}".format(clone_opts, spack_branch,spack_url)
            sexe(clone_cmd, echo=True)

            # optionally, check out a specific commit
            if "spack_commit" in self.project_opts:
                sha1 = self.project_opts["spack_commit"]
                print("[info: using spack commit {}]".format(sha1))
                os.chdir(self.dest_spack)
                sexe("git checkout {}".format(sha1),echo=True)
                
        if self.opts["repo_pull"]:
            # do a pull to make sure we have the latest
            os.chdir(self.dest_spack)
            sexe("git stash", echo=True)
            sexe("git pull", echo=True)

    def config_dir(self):
        """ path to compilers.yaml, which we will use for spack's compiler setup"""
        spack_config_dir = self.opts["spack_config_dir"]
        if spack_config_dir is None:
            uberenv_plat = self.detect_platform()
            if not uberenv_plat is None:
                spack_config_dir = os.path.abspath(pjoin(self.uberenv_path,"spack_configs",uberenv_plat))
        return spack_config_dir


    def disable_spack_config_scopes(self,spack_dir):
        # disables all config scopes except "default", which we will
        # force our settings into
        spack_lib_config = pjoin(spack_dir,"lib","spack","spack","config.py")
        print("[disabling config scope (except default) in: {}]".format(spack_lib_config))
        cfg_script = open(spack_lib_config).read()
        for cfg_scope_stmt in ["('system', os.path.join(spack.paths.system_etc_path, 'spack')),",
                            "('site', os.path.join(spack.paths.etc_path, 'spack')),",
                            "('user', spack.paths.user_config_path)"]:
            cfg_script = cfg_script.replace(cfg_scope_stmt,
                                            "#DISABLED BY UBERENV: " + cfg_scope_stmt)
        open(spack_lib_config,"w").write(cfg_script)



    def patch(self):

        cfg_dir = self.config_dir()
        spack_dir = self.dest_spack

        # force spack to use only default config scope
        self.disable_spack_config_scopes(spack_dir)
        spack_etc_defaults_dir = pjoin(spack_dir,"etc","spack","defaults")

        # copy in default config.yaml
        config_yaml = os.path.abspath(pjoin(self.uberenv_path,"spack_configs","config.yaml"))
        sexe("cp {} {}/".format(config_yaml, spack_etc_defaults_dir ), echo=True)

        # copy in other settings per platform
        if not cfg_dir is None:
            print("[copying uberenv compiler and packages settings from {0}]".format(cfg_dir))

            config_yaml    = pjoin(cfg_dir,"config.yaml")
            compilers_yaml = pjoin(cfg_dir,"compilers.yaml")
            packages_yaml  = pjoin(cfg_dir,"packages.yaml")

            if os.path.isfile(config_yaml):
                sexe("cp {} {}/".format(config_yaml , spack_etc_defaults_dir ), echo=True)

            if os.path.isfile(compilers_yaml):
                sexe("cp {} {}/".format(compilers_yaml, spack_etc_defaults_dir ), echo=True)

            if os.path.isfile(packages_yaml):
                sexe("cp {} {}/".format(packages_yaml, spack_etc_defaults_dir ), echo=True)
        else:
            # let spack try to auto find compilers
            sexe("spack/bin/spack compiler find", echo=True)

        dest_spack_pkgs = pjoin(spack_dir,"var","spack","repos","builtin","packages")
        # hot-copy our packages into spack
        sexe("cp -Rf {} {}".format(self.pkgs,dest_spack_pkgs))


    def clean_build(self):
        # clean out any temporary spack build stages
        cln_cmd = "spack/bin/spack clean "
        res = sexe(cln_cmd, echo=True)

        # clean out any spack cached downloads
        cln_cmd = "spack/bin/spack clean -d"
        res = sexe(cln_cmd, echo=True)

        # check if we need to force uninstall of selected packages
        if self.opts["spack_clean"]:
            if self.project_opts.has_key("spack_clean_packages"):
                for cln_pkg in self.project_opts["spack_clean_packages"]:
                    if not self.find_spack_pkg_path(cln_pkg) is None:
                        unist_cmd = "spack/bin/spack uninstall -f -y --all --dependents " + cln_pkg
                        res = sexe(unist_cmd, echo=True)

    def show_info(self):
        spec_cmd = "spack/bin/spack spec " + self.pkg_name + self.opts["spec"]
        return sexe(spec_cmd, echo=True)

    def install(self):
        # use the uberenv package to trigger the right builds
        # and build an host-config.cmake file
        install_cmd = "spack/bin/spack "
        if self.opts["ignore_ssl_errors"]:
            install_cmd += "-k "
        install_cmd += "install "
        if self.opts["run_tests"]:
            install_cmd += "--test=root "
        install_cmd += self.pkg_name + self.opts["spec"]
        res = sexe(install_cmd, echo=True)
        if res != 0:
            return res

        if "spack_activate" in self.project_opts:
            print("[activating dependent packages]")
            # get the full spack spec for our project
            full_spec = self.read_spack_full_spec(self.pkg_name,self.opts["spec"])
            pkg_names = self.project_opts["spack_activate"].keys()
            for pkg_name in pkg_names:
                pkg_spec_requirements = self.project_opts["spack_activate"][pkg_name]
                activate=True
                for req in pkg_spec_requirements:
                    if req not in full_spec:
                        activate=False
                        break
                if activate:
                    activate_cmd = "spack/bin/spack activate " + pkg_name
                    sexe(activate_cmd, echo=True)
        # note: this assumes package extends python when +python
        # this may fail general cases
        if self.opts["install"] and "+python" in full_spec:
            activate_cmd = "spack/bin/spack activate " + self.pkg_name
            sexe(activate_cmd, echo=True)
        # if user opt'd for an install, we want to symlink the final ascent
        # install to an easy place:
        if self.opts["install"]:
            pkg_path = self.find_spack_pkg_path(self.pkg_name)
            if self.pkg_name != pkg_path["name"]:
                print("[ERROR: Could not find install of {}]".format(self.pkg_name))
                return -1
            else:
                pkg_lnk_dir = "{}-install".format(self.pkg_name)
                if os.path.islink(pkg_lnk_dir):
                    os.unlink(pkg_lnk_dir)
                print("")
                print("[symlinking install to {}]").format(pjoin(self.dest_dir,pkg_lnk_dir))
                os.symlink(pkg_path["path"],os.path.abspath(pkg_lnk_dir))
                hcfg_glob = glob.glob(pjoin(pkg_lnk_dir,"*.cmake"))
                if len(hcfg_glob) > 0:
                    hcfg_path  = hcfg_glob[0]
                    hcfg_fname = os.path.split(hcfg_path)[1]
                    if os.path.islink(hcfg_fname):
                        os.unlink(hcfg_fname)
                    print("[symlinking host config file to {}]".format(pjoin(self.dest_dir,hcfg_fname)))
                    os.symlink(hcfg_path,hcfg_fname)
                print("")
                print("[install complete!]")

    def get_mirror_path(self):
        mirror_path = self.opts["mirror"]
        if not mirror_path:
            print("[--create-mirror requires a mirror directory]")
            sys.exit(-1)
        return os.path.abspath(mirror_path)

    def create_mirror(self):
        """
        Creates a spack mirror for pkg_name at mirror_path.
        """

        mirror_path = self.get_mirror_path()

        mirror_cmd = "spack/bin/spack "
        if self.opts["ignore_ssl_errors"]:
            mirror_cmd += "-k "
        mirror_cmd += "mirror create -d {} --dependencies {}".format(mirror_path,
                                                                    self.pkg_name)
        return sexe(mirror_cmd, echo=True)

    def find_spack_mirror(self, spack_dir, mirror_name):
        """
        Returns the path of a site scoped spack mirror with the
        given name, or None if no mirror exists.
        """
        rv, res = sexe("spack/bin/spack mirror list", ret_output=True)
        mirror_path = None
        for mirror in res.split('\n'):
            if mirror:
                parts = mirror.split()
                if parts[0] == mirror_name:
                    mirror_path = parts[1]
        return mirror_path

    def use_mirror(self):
        """
        Configures spack to use mirror at a given path.
        """
        spack_dir = self.dest_spack
        mirror_name = self.pkg_name
        mirror_path = self.get_mirror_path()
        existing_mirror_path = self.find_spack_mirror(spack_dir, mirror_name)

        if existing_mirror_path and mirror_path != existing_mirror_path:
            # Existing mirror has different URL, error out
            print("[removing existing spack mirror `{}` @ {}]".format(mirror_name,
                                                                    existing_mirror_path))
            #
            # Note: In this case, spack says it removes the mirror, but we still
            # get errors when we try to add a new one, sounds like a bug
            #
            sexe("spack/bin/spack mirror remove --scope=site {} ".format(mirror_name),
                echo=True)
            existing_mirror_path = None
        if not existing_mirror_path:
            # Add if not already there
            sexe("spack/bin/spack mirror add --scope=site {} {}".format(
                    mirror_name, mirror_path), echo=True)
            print("[using mirror {}]".format(mirror_path))


def find_osx_sdks():
    """
    Finds installed osx sdks, returns dict mapping version to file system path
    """
    res = {}
    sdks = glob.glob("/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX*.sdk")
    for sdk in sdks:
        sdk_base = os.path.split(sdk)[1]
        ver = sdk_base[len("MacOSX"):sdk_base.rfind(".")]
        res[ver] = sdk
    return res

def setup_osx_sdk_env_vars():
    """
    Finds installed osx sdks, returns dict mapping version to file system path
    """
    # find current osx version (10.11.6)
    dep_tgt = platform.mac_ver()[0]
    # sdk file names use short version (ex: 10.11)
    dep_tgt_short = dep_tgt[:dep_tgt.rfind(".")]
    # find installed sdks, ideally we want the sdk that matches the current os
    sdk_root = None
    sdks = find_osx_sdks()
    if dep_tgt_short in sdks.keys():
        # matches our osx, use this one
        sdk_root = sdks[dep_tgt_short]
    elif len(sdks) > 0:
        # for now, choose first one:
        dep_tgt  = sdks.keys()[0]
        sdk_root = sdks[dep_tgt]
    else:
        # no valid sdks, error out
        print("[ERROR: Could not find OSX SDK @ /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/]")
        sys.exit(-1)

    env["MACOSX_DEPLOYMENT_TARGET"] = dep_tgt
    env["SDKROOT"] = sdk_root
    print("[setting MACOSX_DEPLOYMENT_TARGET to {}]".format(env["MACOSX_DEPLOYMENT_TARGET"]))
    print("[setting SDKROOT to {}]".format(env[ "SDKROOT"]))



def main():
    """
    Clones and runs a package manager to setup third_party libs.
    Also creates a host-config.cmake file that can be used by our project.
    """

    # parse args from command line
    opts, extra_opts = parse_args()

    # Initialize the environment -- use vcpkg on windows, spack otherwise
    env = SpackEnv(opts, extra_opts) if not is_windows() else VcpkgEnv(opts, extra_opts)

    # Setup the necessary paths and directories
    env.setup_paths_and_dirs()

    # Clone the package manager
    env.clone_repo()

    os.chdir(env.dest_dir)

    # Patch the package manager, as necessary
    env.patch()

    # Clean the build
    env.clean_build()

    # Show the spec for what will be built
    env.show_info()


    ##########################################################
    # we now have an instance of spack configured how we
    # need it to build our tpls at this point there are two
    # possible next steps:
    #
    # *) create a mirror of the packages
    #   OR
    # *) build
    #
    ##########################################################
    if opts["create_mirror"]:
        return env.create_mirror()
    else:
        if not opts["mirror"] is None:
            env.use_mirror()
        res = env.install()

        return res

if __name__ == "__main__":
    sys.exit(main())



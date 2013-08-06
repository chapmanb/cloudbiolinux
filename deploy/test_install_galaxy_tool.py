"""
Test script for building Python API for installing Galaxy tools using
CBL without any dependencies (i.e. it clones down CBL and utilizes it
like bcbio-nextgen's installer).

Goal is to ultimately fold something like this to Galaxy tool shed
client code to provide high-level support for easy CloudBioLinux based
tool installations as @chapmanb described at the 2013 BOSC Codefest.

  <action type="cloudbiolinux_install"
          [cbl_revision="<cbl_git_changeset(default=master)>"]
          [cbl_url="<cbl_repo_url(default=https://github.com/chapmanb/cloudbiolinux)>"]
          [tool_name="<tool_name(default=use dependency package name)>"]
          [tool_version="<tool_version(default=use dependency package version)>"]
          />

"""

import os
from subprocess import check_call
from tempfile import mkdtemp
from getpass import getuser


DEFAULT_CBL_URL = "https://github.com/chapmanb/cloudbiolinux.git"


def __clone_cloudbiolinux(cbl_config):
    """Clone CloudBioLinux to a temporary directory.

    TODO: Support particular revision.
    """
    cbl_url = cbl_config.get("repository", DEFAULT_CBL_URL)
    cbl_dir = mkdtemp(suffix="cbl")
    check_call(["git", "clone", cbl_url, cbl_dir])

    revision = cbl_config.get("revision", None)
    if revision:
        git_dir = os.path.join(cbl_dir, ".git")
        check_call(["git", "--work-tree", cbl_dir, "--git-dir", git_dir, "checkout", revision])
    return cbl_dir


def install_cbl_tool(tool_name, tool_version, install_dir, cbl_config={}):
    cbl_dir = __clone_cloudbiolinux(cbl_config)
    cbl_install_command = [os.path.join(cbl_dir, "deploy", "deploy.sh"), "--action", "install_galaxy_tool"]
    deployer_args = {"vm_provider": "novm",
                     "galaxy_tool_name": tool_name,
                     "galaxy_tool_version": tool_version,
                     "galaxy_tool_dir": install_dir,
                     "settings": "__none__"}
    for prop, val in deployer_args.iteritems():
        cbl_install_command.append("--%s" % prop)
        cbl_install_command.append(val)

    fabric_properties = {"use_sudo": "False",
                         "galaxy_user": getuser()}
    for prop, val in fabric_properties.iteritems():
        cbl_install_command.append("--fabric_property")
        cbl_install_command.append(prop)
        cbl_install_command.append("--fabric_value")
        cbl_install_command.append(val)
    check_call(cbl_install_command)

cbl_config = {"repository": "https://github.com/jmchilton/cloudbiolinux.git"}
install_cbl_tool("tint_proteomics_scripts", "1.19.20", os.path.abspath("test_tool_dir"), cbl_config)

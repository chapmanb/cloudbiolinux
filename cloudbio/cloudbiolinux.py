"""CloudBioLinux specific scripts
"""
import os
from fabric.api import *
from fabric.contrib.files import *

from cloudbio.custom import shared

def _freenx_scripts(env):
    """Provide graphical access to clients via FreeNX.
    """
    home_dir = env.safe_run_output("echo $HOME")
    setup_script = "setupnx.sh"
    bin_dir = shared._get_bin_dir(env)
    install_file_dir = os.path.join(env.config_dir, os.pardir, "installed_files")
    if not env.safe_exists(os.path.join(bin_dir, setup_script)):
        env.safe_put(os.path.join(install_file_dir, setup_script),
                     os.path.join(home_dir, setup_script))
        env.safe_run("chmod 0777 %s" % os.path.join(home_dir, setup_script))
        env.safe_sudo("mv %s %s" % (os.path.join(home_dir, setup_script), bin_dir))
    remote_login = "configure_freenx.sh"
    if not env.safe_exists(os.path.join(home_dir, remote_login)):
        env.safe_put(os.path.join(install_file_dir, 'bash_login'), os.path.join(home_dir, remote_login))
        env.safe_run("chmod 0777 %s" % os.path.join(home_dir, remote_login))
    _configure_gnome(env)

def _cleanup_space(env):
    """Cleanup to recover space from builds and packages.
    """
    env.logger.info("Cleaning up space from package builds")
    with settings(warn_only=True):
        env.safe_sudo("rm -rf .cpanm")
        env.safe_sudo("rm -f /var/crash/*")
        env.safe_run("rm -f ~/*.dot")
        env.safe_run("rm -f ~/*.log")

def _configure_gnome(env):
    """Configure NX server to use classic GNOME.

    http://askubuntu.com/questions/50503/why-do-i-get-unity-instead-of-classic-when-using-nx
    http://notepad2.blogspot.com/2012/04/install-freenx-server-on-ubuntu-1110.html
    """
    add = 'COMMAND_START_GNOME="gnome-session --session gnome-fallback"'
    fname = "/etc/nxserver/node.conf"
    if env.safe_exists("/etc/nxserver/"):
        env.safe_append(fname, add, use_sudo=True)

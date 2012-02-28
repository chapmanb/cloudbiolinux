"""CloudBioLinux specific scripts
"""
import os
from fabric.api import *
from fabric.contrib.files import *

def _freenx_scripts(env):
    """Provide graphical access to clients via FreeNX.
    """
    setup_script = "setupnx.sh"
    remote_setup = "%s/bin/%s" % (env.system_install, setup_script)
    if not exists(os.path.dirname(remote_setup)):
        sudo('mkdir -p {0}'.format(os.path.dirname(remote_setup)))
    install_file_dir = os.path.join(env.config_dir, os.pardir, "installed_files")
    if not exists(remote_setup):
        put(os.path.join(install_file_dir, setup_script), setup_script,
                mode=0777)
        env.safe_sudo("mv %s %s" % (setup_script, remote_setup))
    remote_login = "configure_freenx.sh"
    if not exists(remote_login):
        put(os.path.join(install_file_dir, 'bash_login'), remote_login,
                mode=0777)
    _configure_gnome(env)

def _cleanup_space(env):
    """Cleanup to recover space from builds and packages.
    """
    env.logger.info("Cleaning up space from package builds")
    env.safe_sudo("rm -rf .cpanm")
    env.safe_sudo("rm -f /var/crash/*")

def _configure_gnome(env):
    """Configure NX server to use classic GNOME.

    http://askubuntu.com/questions/50503/why-do-i-get-unity-instead-of-classic-when-using-nx
    """
    add = 'COMMAND_START_GNOME="/usr/bin/nx-session-launcher-suid gnome-session --session=gnome-classic"'
    fname = "/etc/nxserver/node.conf"
    append(fname, add, use_sudo=True)

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
    if env.edition.short_name not in ["minimal"]:
        env.logger.info("Cleaning up space from package builds")
        env.safe_sudo("rm -rf .cpanm")
        env.safe_sudo("rm -f /var/crash/*")
        run("rm -f ~/*.dot")
        run("rm -f ~/*.log")

def _configure_gnome(env):
    """Configure NX server to use classic GNOME.

    http://askubuntu.com/questions/50503/why-do-i-get-unity-instead-of-classic-when-using-nx
    http://notepad2.blogspot.com/2012/04/install-freenx-server-on-ubuntu-1110.html
    """
    add = 'COMMAND_START_GNOME="gnome-session --session gnome-fallback"'
    fname = "/etc/nxserver/node.conf"
    if exists("/etc/nxserver/"):
        append(fname, add, use_sudo=True)

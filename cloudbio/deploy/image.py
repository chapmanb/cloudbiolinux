# Based almost entirely on version from Dr. Enis Afgan at
# (https://bitbucket.org/afgane/mi-deployment)
import os
import os.path

from fabric.api import sudo
from fabric.contrib.files import exists, contains, append

from util import ensure_can_sudo_into, start_service

from cloudbio.galaxy import _setup_users, _setup_xvfb, _install_nginx_standalone, _setup_postgresql
from cloudbio.galaxy.utils import _chown_galaxy
from cloudbio.package import _configure_and_install_native_packages


def configure_MI(env):
    # Clean this next line up.
    _configure_and_install_native_packages(env, ["minimal", "cloudman", "galaxy"])
    # _update_system()
    _setup_users(env)
    _setup_xvfb(env)
    _required_programs(env)


# == required programs
def _required_programs(env):
    """ Install required programs """
    if not exists(env.install_dir):
        sudo("mkdir -p %s" % env.install_dir)
        sudo("chown %s %s" % (env.user, env.install_dir))

    # Setup global environment for all users
    install_dir = os.path.split(env.install_dir)[0]
    exports = ["export PATH=%s/bin:%s/sbin:$PATH" % (install_dir, install_dir),
               "export LD_LIBRARY_PATH=%s/lib" % install_dir]
    for e in exports:
        _ensure_export(e)
    # Install required programs
    _install_nginx_standalone(env)
    _start_nginx(env)
    _deploy_setup_postgresql(env)

    # Verify this is not needed.
    # _install_samtools()


def _ensure_export(command):
    if not contains('/etc/bash.bashrc', command):
        append('/etc/bash.bashrc', command, use_sudo=True)


def _start_nginx(env):
    galaxy_data = env.galaxy_data_mount
    env.safe_sudo("mkdir -p '%s'" % env.galaxy_data)
    _chown_galaxy(env, galaxy_data)
    start_service("nginx")


def _deploy_setup_postgresql(env):
    ensure_can_sudo_into("postgres")
    _setup_postgresql(env)

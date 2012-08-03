"""
"""

from cloudbio.galaxy import _setup_users
from cloudbio.galaxy import _setup_galaxy_env_defaults
from cloudbio.galaxy import _install_galaxy
from cloudbio.galaxy import _configure_galaxy_options
from cloudbio.galaxy.tools import _install_tools


def install_galaxy_webapp(env):
    _prep_galaxy(env)
    _install_galaxy(env)
    _configure_galaxy_options(env)


def install_galaxy_tools(env):
    _prep_galaxy(env)
    _install_tools(env)


def _prep_galaxy(env):
    _setup_users(env)
    _setup_galaxy_env_defaults(env)

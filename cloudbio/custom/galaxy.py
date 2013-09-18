"""
Install any components that fall under 'galaxy' directive in main.yaml
"""
from cloudbio.galaxy import _setup_users
from cloudbio.galaxy import _setup_galaxy_env_defaults
from cloudbio.galaxy import _install_galaxy
from cloudbio.galaxy import _configure_galaxy_options


def install_galaxy_webapp(env):
    _prep_galaxy(env)
    _install_galaxy(env)
    _configure_galaxy_options(env)


def _prep_galaxy(env):
    _setup_users(env)
    _setup_galaxy_env_defaults(env)

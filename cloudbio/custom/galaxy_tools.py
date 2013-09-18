"""
Install any components that fall under 'galaxy_tools' directive in main.yaml
"""
from cloudbio.galaxy.tools import _install_tools
from cloudbio.custom.galaxy import _prep_galaxy


def install_cbl_galaxy_tools(env):
    _prep_galaxy(env)
    _install_tools(env)

"""
Millstone flavor.
"""

from fabric.api import *
from fabric.contrib.files import *

from cloudbio.custom import shared
from cloudbio.flavor import Flavor

from installer import install_millstone


class MillstoneFlavor(Flavor):
    def __init__(self, env):
        Flavor.__init__(self, env)
        self.name = "Millstone Flavor"

    def rewrite_config_items(self, name, packages):
        return packages

    def post_install(self):
        # Manually install logrotate, postgresql-9.3, rabbitmq-server.
        # postgresql-9.3 and rabbitmq-server require incompatible versions of
        # logrotate by default, so we force install the working version.
        env.safe_sudo("apt-get -y --force-yes install logrotate=3.7.8-6ubuntu5")
        env.safe_sudo("apt-get -y --force-yes install postgresql-9.3")
        env.safe_sudo("apt-get -y --force-yes install libpq-dev")
        env.safe_sudo("apt-get -y --force-yes install pgdg-keyring")
        env.safe_sudo("apt-get -y --force-yes install rabbitmq-server")

        # Now install psycogpg2, which requires postgresql-9.3 and related
        # to have been installed first.
        env.safe_sudo("{0} install --upgrade {1}".format(shared._pip_cmd(env), 'psycopg2'))

        # Finally, install Millstone.
        install_millstone(self.env)

env.flavor = MillstoneFlavor(env)

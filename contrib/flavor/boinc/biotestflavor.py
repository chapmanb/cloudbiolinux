from fabric.api import *
from fabric.contrib.files import *

from cloudbio.flavor import Flavor

from cloudbio.custom.shared import (_fetch_and_unpack)

class BoincFlavor(Flavor):
    """A VM flavor for running Boinc
    """
    def __init__(self, env):
        Flavor.__init__(self,env)
        self.name = "Boinc client"

    def rewrite_config_items(self, name, items):
        return []

    def post_install(self):
        env.logger.info("Starting post-install")

env.flavor = BoincFlavor(env)

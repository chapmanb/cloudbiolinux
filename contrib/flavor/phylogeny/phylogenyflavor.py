from fabric.api import *
from fabric.contrib.files import *

from cloudbio.flavor import Flavor

from cloudbio.custom.shared import (_fetch_and_unpack)

class PhylogenyFlavor(Flavor):
    """A VM flavor for running Phylogeny
    """
    def __init__(self, env):
        Flavor.__init__(self,env)
        self.name = "Phylogeny Flavor"

    def rewrite_config_items(self, name, packages):
        if name == 'packages':
          packages += [ 'openssh-server', 'unzip', 'tar', 'sudo', 'openjdk-6-jre']
          packages += [ 'openmpi-bin' ]  # required for MrBayes-MPI
          # if 'bio-linux-mrbayes-multi' in packages:
          #   (Debian version is still not OK)
          #   packages.remove('bio-linux-mrbayes-multi')

        for package in packages:
          env.logger.info("Selected: "+name+" "+package)
        return packages

    def post_install(self):
        env.logger.info("Starting post-install")
        pass

env.flavor = PhylogenyFlavor(env)

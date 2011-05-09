from fabric.api import *

from cloudbio.edition import Edition
import logging

class BioNode(Edition):
    """BioNode specialization of BioLinux
    """
    def __init__(self, env):
        self.name = "BioNode Edition"
        self.env = env
        self.include_oracle_virtualbox = False
        self.include_freenx = False
        self.logger = logging.getLogger("cloudbiolinux")

    def check_packages_source(self):
        # Bionode removes sources, just to be sure
        self.logger.debug("Clearing %s" % self.env.sources_file)
        sudo("cat /dev/null > %s" % self.env.sources_file)


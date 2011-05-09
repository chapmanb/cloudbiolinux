from fabric.api import *

from cloudbio.edition import Edition

class BioNode(Edition):
    """BioNode specialization of BioLinux
    """
    def __init__(self, env):
        self.name = "BioNode Edition"
        self.env = env
        self.include_oracle_virtualbox = False
        self.include_freenx = False

    def check_packages_source(self):
        # Bionode removes sources, just to be sure
        self.env.logger.debug("Clearing %s" % self.env.sources_file)
        sudo("cat /dev/null > %s" % self.env.sources_file)


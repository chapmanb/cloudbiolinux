from fabric.api import *

from cloudbio.edition.base import Edition

class BioNode(Edition):
    """BioNode specialization of BioLinux
    """
    def __init__(self, env):
        Edition.__init__(self,env)
        self.name = "BioNode Edition"
        self.short_name = "bionode"
        self.env = env

    def check_packages_source(self):
        # Bionode removes sources, just to be sure
        self.env.logger.debug("Clearing %s" % self.env.sources_file)
        sudo("cat /dev/null > %s" % self.env.sources_file)


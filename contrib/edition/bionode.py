
from fabric.api import *

from edition import *
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
# create formatter
formatter = logging.Formatter('%(name)s %(levelname)s: %(message)s')
# add formatter to ch
ch.setFormatter(formatter)
logger.addHandler(ch)

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
        logger.debug("Clearing %s" % self.env.sources_file)
        sudo("cat /dev/null > %s" % self.env.sources_file)


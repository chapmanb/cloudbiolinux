from fabric.api import *

from cloudbio.flavor import Flavor

class BioTestFlavor(Flavor):
    """A Flavor for cross Bio* tests
    """
    def __init__(self, env):
        Flavor.__init__(self,env)
        self.name = "Bio* cross-lang flavor"

env.flavor = BioTestFlavor(env)

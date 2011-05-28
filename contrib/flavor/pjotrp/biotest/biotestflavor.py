from fabric.api import *

from cloudbio.flavor import Flavor

class BioTestFlavor(Flavor):
    """A Flavor for cross Bio* tests
    """
    def __init__(self, env):
        Flavor.__init__(self,env)
        self.name = "Bio* cross-lang flavor"

    def rewrite_packages_list(self, list):
        # list.remove('screen')
        # list.append('test')
        return list

env.flavor = BioTestFlavor(env)

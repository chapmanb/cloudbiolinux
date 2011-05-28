from fabric.api import *
from fabric.contrib.files import *

from cloudbio.flavor import Flavor

from cloudbio.custom.shared import (_fetch_and_unpack)

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

    def post_install(self):
        env.logger.info("Starting post-install")
        if exists('Scalability'):
            with cd('Scalability'):
               run('git pull')
        else:
           _fetch_and_unpack("git clone git://github.com/pjotrp/Scalability.git")
        # Now run a post installation routine
        run('./Scalability/scripts/hello.sh')


env.flavor = BioTestFlavor(env)

from fabric.api import *
from fabric.contrib.files import *

from cloudbio.flavor import Flavor

# This flavour installs several neuroinformatics software from NeuroDebian
#
#   http://neuro.debian.net/
#
# Author:  Roman Valls Guimera <roman@incf.org>

class NeuroFlavor(Flavor):
    """ A flavour for installing NeuroDebian. A debian neuroinformatics repository
        of software and datasets.
    """

    def __init__(self, env):
        Flavor.__init__(self, env)
        self.name = "Neuroinformatics Flavor"

    def rewrite_config_items(self, name, packages):
        if name == 'packages':
            packages.extend([
                "git"
            ])
        return packages


    def rewrite_apt_sources_list(self, name, sources):
        sources = [
            'deb http://neurodeb.pirsquared.org data main contrib non-free',
            '#deb-src http://neurodeb.pirsquared.org data main contrib non-free',
            'deb http://neurodeb.pirsquared.org saucy main contrib non-free',
            '#deb-src http://neurodeb.pirsquared.org saucy main contrib non-free'
        ]

#wget -O- http://neuro.debian.net/lists/saucy.us-ca.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
#sudo apt-key adv --recv-keys --keyserver pgp.mit.edu 2649A5A9

env.flavor = NeuroFlavor(env)
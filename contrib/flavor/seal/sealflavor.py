from fabric.api import *
from fabric.contrib.files import *

from cloudbio.flavor import Flavor

from cloudbio.custom.shared import (_fetch_and_unpack)

import sys

# This flavour installs the Seal toolkit for processing high-throughput
# sequencing data on Hadoop.
#   http://biodoop-seal.sf.net/
#
# It pulls in quite a few dependencies, including Hadoop itself and
# Pydoop (http://pydoop.sf.net/).
#
# The dependencies it pulls into Cloudbiolinux are structured as follows:
#
# contrib/flavor/seal/main.yaml
#   sealdist
#   customsealdist
#
# config/packages-yum.yaml
#   sealdist (metapackage)
# config/custom.yaml
#   customsealdist (metapackage)
#     - pydoop
#     - seal
#
# The components of the customsealdist metapackage are installed through
# the functions in cloudbio/custom/customsealdist.py
#
#
# This flavour has only been installed on Scientific Linux and has not
# yet been well tested.
#
# To try installing it run the following:
#   cd <your cloudbiolinux directory>
#   fab -f ./fabfile.py -H root@<your host> -c ./contrib/flavor/seal/fabricrc_sl.txt  install_biolinux:packagelist=contrib/flavor/seal/main.yaml
#
# Authors:  Roman Valls Guimera <roman.valls.guimera@scilifelab.se>
#           Luca Pireddu <luca.pireddu@crs4.it>

class SealFlavor(Flavor):
	"""A flavour for installing Seal
	"""
	def __init__(self, env):
		Flavor.__init__(self,env)
		self.name = "Seal Flavor"

	def rewrite_config_items(self, name, packages):
		if name == 'packages':
			if sys.version_info < (2,7):
				# for versions of Python prior to 2.7 we need to add importlib
				# and argparse
				packages.extend([ 
					"python-importlib",
					"python-argparse"
				])
		return packages


	def post_install(self):
		env.logger.info("Starting post-install")
		pass

env.flavor = SealFlavor(env)
